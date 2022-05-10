#' Check data has the correct variables
#' @noRd
assert_particle_data_df <- function(data){
  pars_not_in <- setdiff(c("deaths", "date_start"), names(data))
  if(length(pars_not_in) > 1){
    stop(paste0("data is missing variables: ", paste0(pars_not_in, collapse = ", ")))
  }
  squire:::assert_int(data$deaths)
  squire:::assert_date(data$date_start)
  squire:::assert_date(data$date_end)
  #check for overlaps and inconsistencies
  if(any(data$date_start >= data$date_end)){
    stop("date_start greater than or equal to date_end in data")
  }
  if(any(diff(data$date_start) < as.numeric(data$date_start - data$date_end)[-length(data$date_start)])){
    stop("A time period covered in data overlaps another time period")
  }
}
#' check distribution is in the right format
#' @noRd
assert_distribution <- function(distribution){
  #each element should have the same set of names
  par_names <- purrr::map(distribution, ~names(.x))
  if(!purrr::every(par_names, ~all(.x %in% par_names[[1]]))){
    stop("Every element of distribution must have the same set of named parameters")
  }
}
#' Function factory for setting up the model output generating function for the model
#' should work for all types, with the caveat that there must be a psedo method in
#' squire:::beta_est that works with the model.
#' @noRd
generate_model_function <- function(squire_model, parameters){
  squire_parameters <- do.call(squire_model$parameter_func, parameters)
  #define a function that gets the deaths for a requested period
  function(Rt, t_start, t_end, initial_state = NULL, tt_Rt = NULL){
    force(squire_parameters)
    force(squire_model)
    #setup the initial state for the odin model
    if(!is.null(initial_state)){
      #this might need to be model specific!
      #set all in initial state so that t = 0 corresponds to t_start
      pars_to_change <- names(initial_state)[stringr::str_detect(names(initial_state), "tt_")]
      for(par in pars_to_change){
        if(length(initial_state[[par]]) > 1){
          #if the parameter is not constant
          initial_state[[par]] <- initial_state[[par]] - t_start
        }
      }
      squire_parameters <- initial_state
    }
    #calculate beta value for Rt
    squire_parameters$beta_set <-
      squire:::beta_est(squire_model, squire_parameters, Rt)
    if(!is.null(tt_Rt)){
      squire_parameters$tt_beta <- tt_Rt
    }
    #create odin model with these parameters
    odin_model <- squire_model$odin_model$new(
      user = squire_parameters,
      unused_user_action = "ignore"
    )
    #run the model over the requested timeperiod
    model_output <- odin_model$run(c(0, seq(t_start, t_end, by = 1)))[-1,]
  }
}
#' Function factory for setting up the death generating function for the model
#' should for all types, with the caveat that there must be a psedo method in
#' squire:::beta_est that works with the model.
#' @noRd
generate_deaths_function <- function(model_func, squire_model, parameters){
  N_age <- do.call(squire_model$parameter_func, parameters)$N_age
  #define a function that gets the deaths for a requested period
  get_deaths <- function(Rt, t_start, t_end, initial_state = NULL){
    force(model_func)
    force(N_age)
    model_output <- model_func(Rt, t_start, t_end, initial_state)
    model_output[, paste0("D[", seq_len(N_age), "]")] %>%
      rowSums()
  }
}

#' Spread initial infections out over a model
#' @noRd
assign_infections <- function(initial_state, initial_infections){
  #update E1_0, might need to be model specific!
  #spread evenly across adult population
  initial_state$E1_0[4:14] <- initial_infections * initial_state$population[4:14]/sum(initial_state$population[4:14])
  initial_state
}

#' Update initial state with model output
#' @noRd
update_initial_state <- function(initial_state, model_output){
  model_output <- utils::tail(model_output, 1)
  #get_values to update
  pars <- stringr::str_replace(names(initial_state)[stringr::str_detect(names(initial_state), "_0")], "_0", "")
  names(pars) <- pars
  new_values <- purrr::map(
    pars,
    function(par){
      par_names <- paste0(par, "[", 1:initial_state$N_age, "]")
      as.numeric(model_output[, par_names])
    }
  )
  #update the values
  initial_state[paste0(names(new_values), "_0")] <- new_values
  initial_state
}

#' Get the delay in the Rt effects on death
#' A change in Rt has a delay until it effects deaths
#' For now we will set it to be the shortest possible time
#' @noRd
get_delay <- function(squire_model, parameters){
  pars <- do.call(squire_model$parameter_func, parameters)
  vals <- pars$dur_ICase + pars$dur_E + unlist(pars[c("dur_get_mv_die", "dur_get_ox_die",
                                                      "dur_not_get_mv_die", "dur_not_get_ox_die")])
  list(min = round(min(vals)), max = round(max(vals)))#max(vals)))
}
#' Get the time series of Rt trend changes and split data up
#' @noRd
get_time_series <- function(squire_model, parameters, data, rt_spacing){
  # Set up delay parameter
  rt_death_delay <- get_delay(squire_model, parameters)
  #first t occures before first death so we ignore 0s
  #calculate when rt should change, so that it aligns with the data
  rt_change_t <- seq(data$t_start[data$deaths > 0][1], utils::tail(data$t_end, 1), by = rt_spacing) - rt_death_delay$min
  #remove any trends that occures on or before 0, we add R0 later anyway
  rt_change_t <- rt_change_t[rt_change_t > 0]
  #ensure no changes estiamted in last max delay period
  to_remove <- which(rt_change_t > utils::tail(data$t_end, 1) - rt_death_delay$max)
  if(length(to_remove) > 0){
    rt_change_t <- c(0, rt_change_t[-to_remove])
  } else {
    rt_change_t <- c(0, rt_change_t)
  }

  rt_df <- data.frame(rt_change_t = rt_change_t)
  rt_df$initial_target <- c(rt_change_t[-1], utils::tail(data$t_end, 1))
  rt_df$death_start_t <- rt_df$rt_change_t + rt_death_delay$min
  rt_df$death_start_t[1] <- 0
  rt_df$death_end_t <- rt_df$initial_target + rt_death_delay$max
  rt_df$death_end_t[length(rt_df$death_end_t)] <- utils::tail(data$t_end, 1)

  split_data <- purrr::map(seq_len(nrow(rt_df)), function(x){
    data %>%
      dplyr::filter(
        .data$t_start >= rt_df$death_start_t[x] &
          .data$t_end <= rt_df$death_end_t[x]
      )
  })
  #Some kind of check on the data?

  list(rt_df = rt_df, split_data = split_data)
}
#' Find the particle that maximises the likelihood for R0 and initial conditions
#' @noRd
particle_optimise_R0 <- function(initial_r, initial_infections, n_particles, proposal_width, R0_likelihood){
  repeat_limit <- 5
  repeats <- 0
  searching <- TRUE
  while(searching){
    values <- data.frame(
      R0 = get_Rt_to_explore(initial_r, proposal_width, n_particles),
      initial_infections = seq(initial_infections * (1-proposal_width), initial_infections/(1-proposal_width), length.out = n_particles)
    ) %>%
      tidyr::expand(.data$R0, .data$initial_infections) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        ll = R0_likelihood(.data$R0, .data$initial_infections)
      ) %>%
      dplyr::ungroup()
    #get the maximum and return R0 and initial_infections
    mle <- values %>%
      dplyr::filter(.data$ll == max(.data$ll))
    if(nrow(mle) > 1){
      #catch if we have equal likelihoods
      mle <- mle[round(nrow(mle)/2),]
    }
    update_R0 <- mle$R0 %in% c(min(values$R0), max(values$R0))
    update_initial_infections <-  mle$initial_infections %in%
      c(min(values$initial_infections), max(values$initial_infections))
    if(update_R0){
      initial_r <- mle$R0
    }
    if(update_initial_infections){
      initial_infections <- mle$initial_infections
    }
    if(!update_R0 & !update_initial_infections){
      searching <- FALSE
    } else {
      repeats <- repeats + 1
    }
    if(repeats == repeat_limit){
      searching <- FALSE
    }
  }
  c(R0 = mle$R0, initial_infections = mle$initial_infections)
}
#' Find the particle that optimises the likelihood for Rt
#' @noRd
particle_optimise_Rt <- function(rt, initial_state, rt_index, n_particles, proposal_width, likelihood){
  repeat_limit <- 5
  repeats <- 0
  searching <- TRUE
  while(searching){
    rt_to_explore <- get_Rt_to_explore(rt, proposal_width, n_particles)
    ll <- purrr::map_dbl(rt_to_explore, ~likelihood(.x, rt_index, initial_state))
    best <- rt_to_explore[which.max(ll)]

    update_Rt <- best %in% c(min(rt_to_explore), max(rt_to_explore))

    if(update_Rt){
      rt <- best
      repeats <- repeats + 1
    } else{
      searching <- FALSE
    }
    if(repeats == repeat_limit){
      searching <- FALSE
    }
  }
  best
}
#' Get the rt values to explore
#' @noRd
get_Rt_to_explore <- function(rt, proposal_width, n_particles){
  rt_max <- 20
  upper <- rt/(1-proposal_width)
  if(upper > rt_max){
    upper <- rt_max
  }
  seq(rt * (1-proposal_width), upper, length.out = n_particles)
}
