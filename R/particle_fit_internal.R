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
#' squire:::beta_est that works with the model or a method for squire.page::beta_est.
#' @noRd
generate_model_function <- function(squire_model, parameters){
  #get parameters for the model, essentially fills out missing with defaults
  squire_parameters <- setup_parameters(squire_model, parameters)
  #define a function that gets the deaths for a requested period
  function(Rt, t_start, t_end, initial_state = NULL, tt_Rt = NULL, atol= 10^-6, rtol = 10^-6){
    #ensures the environment contains these objects
    force(squire_parameters)
    force(squire_model)
    #set the initial state for the odin model
    if (!is.null(initial_state)) {
      squire_parameters <- initial_state
    }
    #calculate beta value from Rt
    squire_parameters$beta_set <-
      beta_est(squire_model, squire_parameters, Rt)
    if(!is.null(tt_Rt)){
      squire_parameters$tt_beta <- tt_Rt
    }
    #create odin model with these parameters
    #with catch if nimue type
    if(class(squire_model$odin_model) == "function"){
      odin_model <- squire_model$odin_model(
        user = squire_parameters,
        unused_user_action = "ignore"
      )
    } else {
      odin_model <- squire_model$odin_model$new(
        user = squire_parameters,
        unused_user_action = "ignore"
      )
    }
    #run the model over the requested time-period
    model_output <- odin_model$run(seq(t_start, t_end, by = 1), atol = atol, rtol = rtol)
  }
}
#' Function factory for setting up the death generating function for the model
#' should for all types, with the caveat that there must be a psedo method in
#' squire:::beta_est that works with the model or a method for squire.page::beta_est.
#' @noRd
generate_deaths_function <- function(model_func){
  #determine with columns related to deaths
  temp_rn <- model_func(1, 0, 2, NULL)
  #save there column indexes
  death_cols <- stringr::str_detect(colnames(temp_rn), "D\\[")
  rm(temp_rn)
  #define a function that gets the deaths for a requested period
  get_deaths <- function(Rt, t_start, t_end, initial_state = NULL){
    force(model_func)
    force(death_cols)
    #generate output
    model_output <- model_func(Rt, t_start, t_end, initial_state)
    #reduce to death compartments and reduce
    model_output[, death_cols] %>%
      rowSums()
  }
}

#' Spread initial infections out over a model
#' @noRd
assign_infections <- function(initial_state, initial_infections){
  #Adds infections uniformly across adult population in the E1 compartment
  #this should work regardless of type due to how R handles [] calls to matrices
  initial_state$E1_0[4:14] <- initial_infections * initial_state$population[4:14]/sum(initial_state$population[4:14])
  initial_state
}

#' Update initial state with model output
#' @noRd
update_initial_state <- function(initial_state, model_output){
  #get_values to update, they will have _0 in their name but we remove this to get the final values
  pars <- stringr::str_replace(names(initial_state)[stringr::str_detect(names(initial_state), "_0")], "_0", "")
  names(pars) <- pars
  #get each value from the output
  new_values <- purrr::map(
    pars,
    function(par) {
      cols <- stringr::str_detect(colnames(model_output), paste0(par, "\\["))
      as.numeric(model_output[, cols])
    }
  )
  if("N_vaccine" %in% names(initial_state)){
    #convert to matrices if vaccine model
    new_values <- purrr::map(
      new_values, ~matrix(.x, nrow = initial_state$N_age, ncol = initial_state$N_vaccine)
    )
  }
  #update the values
  initial_state[paste0(names(new_values), "_0")] <- new_values
  initial_state
}

#' Get the delay in the Rt effects on death
#' A change in Rt has a delay until it effects deaths
#' For now we will set it to be the shortest possible time
#' @noRd
get_delay <- function(squire_model, parameters){
  #convert parameters into those fed into odin model
  pars <- setup_parameters(squire_model, parameters)
  #get the range of rate from infection to death (if dieing)
  vals <- (2/pars$gamma_ICase) + (2/pars$gamma_E) + (2/unlist(pars[c("gamma_get_mv_die", "gamma_get_ox_die",
                                                      "gamma_not_get_mv_die", "gamma_not_get_ox_die")]))
  #return the minimum and the maximum
  list(min = floor(min(vals)), max = ceiling(max(vals)))
}
#' Get the time series of Rt trend changes and split data up
#' @noRd
get_time_series <- function(squire_model, parameters, data, rt_spacing){
  # Set up delay parameter
  rt_death_delay <- get_delay(squire_model, parameters)
  #calculate when rt should change, so that it aligns with the data
  rt_change_t <- seq(data$t_start[data$deaths > 0][1], utils::tail(data$t_end, 1), by = rt_spacing) - rt_death_delay$min
  #remove any trends that occures on or before 0, we add R0 later anyway
  rt_change_t <- rt_change_t[rt_change_t > 0]
  #ensure no changes estimated in last max delay period, (i.e. don't estimate an Rt in the last ~20 days)
  to_remove <- which(rt_change_t > utils::tail(data$t_end, 1) - rt_death_delay$max)
  if(length(to_remove) > 0){
    rt_change_t <- c(0, rt_change_t[-to_remove])
  } else {
    rt_change_t <- c(0, rt_change_t)
  }
  #setup data frame
  rt_df <- data.frame(rt_change_t = rt_change_t)
  #add parameter for when the next trend is started (used for simulating up to the next trend)
  rt_df$initial_target <- c(rt_change_t[-1], utils::tail(data$t_end, 1))
  #parameter for when the fitting period starts
  rt_df$death_start_t <- rt_df$rt_change_t + rt_death_delay$min
  rt_df$death_start_t[1] <- 0
  #parameter for when the fitting period ends
  rt_df$death_end_t <- rt_df$initial_target + rt_death_delay$max
  rt_df$death_end_t[length(rt_df$death_end_t)] <- utils::tail(data$t_end, 1)
  #for convience we also split the data up into seperate data frames based on
  #which rt_trend they relate too, note these groups likely overlap and that is
  #intentional
  split_data <- purrr::map(seq_len(nrow(rt_df)), function(x){
    data %>%
      dplyr::filter(
        .data$t_start >= rt_df$death_start_t[x] &
          .data$t_end <= rt_df$death_end_t[x]
      ) %>%
      #format so relative to rt_change_t
      dplyr::mutate(
        index_start = .data$t_start - rt_df$rt_change_t[x] + 1,
        index_end = .data$t_end - rt_df$rt_change_t[x] + 1
      )
  })
  #potentially add a check in case Rt trends have no data to fit too

  list(rt_df = rt_df, split_data = split_data)
}
#' Find the particle that maximises the likelihood for R0 and initial conditions
#' @noRd
particle_optimise_R0 <- function(initial_r, initial_infections, n_particles, proposal_width, R0_likelihood, rt_upper_bound_max = 20, rt_upper_bound_min = 5){
  #this is identical to particle_optimise_Rt except it also explores initial infections
  repeat_limit <- 5
  repeats <- 0
  searching <- TRUE
  while(searching){
    values <- data.frame(
      #get the range of Rt values
      R0 = get_Rt_to_explore(initial_r, proposal_width, n_particles, rt_upper_bound_max, rt_upper_bound_min),
      #get the range of initial infections to explore
      initial_infections = seq(initial_infections * (1-proposal_width), initial_infections/(1-proposal_width), length.out = n_particles)
    ) %>%
      #get all possible compations
      tidyr::expand(.data$R0, .data$initial_infections) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        #calculate their likelihood, with catch if it fails
        ll = tryCatch(R0_likelihood(.data$R0, .data$initial_infections),
                 error = function(e){-Inf} #negative infinity so it is not picked
        )
      ) %>%
      dplyr::ungroup()
    #get the maximum and return R0 and initial_infections
    mle <- values %>%
      dplyr::filter(.data$ll == max(.data$ll))
    if(nrow(mle) > 1){
      #catch if we have equal likelihoods (shouldn't happen, unless data is 0)
      mle <- mle[round(nrow(mle)/2),]
    }
    #check if we need to repeat this (if R0 or initial_infections are on boundaries)
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
particle_optimise_Rt <- function(rt, initial_state, rt_index, n_particles, proposal_width, likelihood, rt_upper_bound_max = 20, rt_upper_bound_min = 5){
  #we allow this search to repeat a maximum of 5 times, if the maximum likelihood
  #is on the boundary of the searched values, then we repeat with it recentered on
  #that value
  repeat_limit <- 5
  repeats <- 0
  searching <- TRUE
  while(searching){
    #get values to explore
    rt_to_explore <- get_Rt_to_explore(rt, proposal_width, n_particles, rt_upper_bound_max, rt_upper_bound_min)
    #calculate R0
    ll <- purrr::map_dbl(rt_to_explore, function(this_rt){
      #add a catch if this fails
      tryCatch(likelihood(this_rt, rt_index, initial_state),
               error = function(e){-Inf} #rerun negative infinity so it is not picked
      )
    })
    #pick the highest
    best <- rt_to_explore[which.max(ll)]
    #check if on a boundary
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
get_Rt_to_explore <- function(rt, proposal_width, n_particles, rt_upper_bound_max = 20, rt_upper_bound_min = 5){
  #set bounds on maximums to avoid ridiculous values or being stuck below 1
  upper <- rt/(1-proposal_width)
  if(upper > rt_upper_bound_max){
    upper <- rt_upper_bound_max
  } else if(upper < rt_upper_bound_min){
    upper <- rt_upper_bound_min
  }
  #generate values
  seq(rt * (1-proposal_width), upper, length.out = n_particles)
}
#' Function to calculate ll with the negative binomial
#'
#' Separate function so it can be used in trimming too
#'
#' @noRd
ll_negative_binomial <- function(model_deaths, data_deaths, k){
  model_deaths <- dplyr::if_else(model_deaths < 0, 10^-10, model_deaths)
  #calculate likelihood for each time-period then sum
  squire:::ll_nbinom(
    data = data_deaths,
    model = model_deaths,
    phi = 1,
    k = k,
    exp_noise = Inf #no noise
  ) %>%
    sum()
}

