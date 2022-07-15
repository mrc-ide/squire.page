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
    if(inherits(squire_model$odin_model, "function")){
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
    #run the model over the requested time-period, with decreasing tolerance if it fails
    ts <- seq(t_start, t_end, by = 1)
    model_output <- tryCatch(
      odin_model$run(ts, atol = atol, rtol = rtol),
      error = function(e){
        tryCatch(
          odin_model$run(ts, atol = atol*0.1, rtol = rtol*0.1),
          error = function(e){
            odin_model$run(ts, atol = atol*0.01, rtol = rtol*0.01)
          }
        )
      }
    )
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
  #set S1 back to the full population
  initial_state$S_0 <- initial_state$S_0 + initial_state$E1_0
  #Adds infections uniformly across adult population in the E1 compartment
  #this should work regardless of type due to how R handles [] calls to matrices
  initial_state$E1_0[4:14] <- initial_infections * initial_state$population[4:14]/sum(initial_state$population[4:14])
  #remove these from s1
  initial_state$S_0 <- initial_state$S_0 - initial_state$E1_0
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
  to_remove <- which(rt_change_t > utils::tail(data$t_end, 1) - (rt_death_delay$max))
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
#' Find the particle that optimises the likelihood for Rt
#' @noRd
optimise_Rt <- function(initial_state, rt_index, likelihood, rt_interval){
  stats::optimize(
    function(x){tryCatch(likelihood(x, rt_index, initial_state),
                    error = function(e){-Inf} #rerun negative infinity so it is not picked
             )},
    upper = rt_interval[2], lower = rt_interval[1], maximum = TRUE
    )$maximum
}
#' Find the particle that maximises the likelihood for R0 and initial conditions
#' @noRd
particle_optimise_R0 <- function(initial_infections_interval, n_particles, likelihood, rt_interval, initial_state){
  #we optimise R0 as with Rt values but pick initial infections as a particle
  #we can avoid multidimensional optimisation

  #get the range of initial infections to explore
  initial_infections <- seq(initial_infections_interval[1], initial_infections_interval[2], length.out = n_particles)
  lls <- purrr::map(
    initial_infections,
    ~stats::optimize(
      function(x){tryCatch(likelihood(x , 1, assign_infections(initial_state, .x)),
                           error = function(e){-Inf} #rerun negative infinity so it is not picked
      )},
      upper = rt_interval[2], lower = rt_interval[1], maximum = TRUE
    )
  )
  mle_index <- which.max(purrr::map_dbl(lls, ~.x$objective))
  c(R0 = lls[[mle_index]]$maximum, initial_infections = initial_infections[mle_index])
}
#' Function to calculate ll with the negative binomial
#'
#' Separate function so it can be used in trimming too
#'
#' @noRd
ll_negative_binomial <- function(model_deaths, data_deaths, k){
  #nb cannot handle negative means nor 0 means so we set a positive floor for the modelled deaths
  floor_value <- 10^-6
  model_deaths <- dplyr::if_else(model_deaths < floor_value, floor_value, as.numeric(model_deaths))
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

