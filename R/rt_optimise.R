#' Fit nimue to a given set of death data via the particle and optimisation exploration approach.
#'
#' For a given number of samples from given parameter uncertainty or distribution
#' function, this approach iteratively fits the Rt trend to the provided death
#' data and returns a nimue_simulation object for future usage in scenario
#' modelling.
#'
#' NOTE: For death curves with periods of 0's rt_interval's lower bound must be greater than 0
#' else it will likely fail to overcome a low infective population and Rt will
#' tend to some unrealistically large number.
#'
#' This function is progressr enabled, so progressr::handlers(global = TRUE)
#' can be used to view progress through the samples.
#'
#' @param data A data frame of deaths occuring over a given time frame. Given in the
#' format: deaths(integer), date_start(date), date_end(date). Must have at least
#' one death period in each set of Rt trend changes (i.e. a 14 day period by
#' default) and no period can overlap these changes.
#' @param distribution A list of samples(a list) with names specifying parameters.
#' @param squire_model A model object of the desired type to use, i.e. squire, nimue.
#' @param parameters A list of parameters to keep the same across all samples.
#' @param start_date Date when the epidemic begins, parameters and distribution
#' should be formatted relative to this date, where necessary.
#'
#' @param parallel Run each sample concurrently, uses the future::plan set by the user. Default = FALSE.
#' @param rt_spacing Number of days between each Rt trend, default = 14 days.
#' @param k Control the dispersion on the negative binomial likelihood, default = 2.
#' @param n_particles How many particles to explore uniformly the interval of initial infections, default = 7.
#' @param initial_infections_interval The range of initial number of infections to explore, default = c(5, 500).
#' @param rt_interval The range of values that Rt can take, default = c(0.5, 20).
#' @param dt If the passed squire_model has a difference model attached, this is the step size we shall use.
#' The difference model is only used if dt is non-null and squire_model$odin_difference_model is non-null.
#' Defaults to NULL.
#'
#' @return An object of type rt_optimised, (model type).
#'
#' @export
rt_optimise <- function(data, distribution, squire_model, parameters,
                         start_date,
                         parallel = FALSE,
                         rt_spacing = 14,
                         k = 2,
                         n_particles = 14,
                         initial_infections_interval = c(5, 500),
                         rt_interval = c(0.5, 20),
                         dt = NULL
                        ) {
  ##Initial Checks and Housekeeping

  #checks
  squire:::assert_dataframe(data)
  #sort data by start date
  data <- dplyr::arrange(data, .data$date_start) %>%
    dplyr::filter(.data$date_start > start_date) #must have occured after start date
  assert_particle_data_df(data)
  if(!inherits(squire_model, "squire_model")){
    stop("squire_model must be an object inheriting from class squire_model")
  }
  assert_distribution(distribution)
  #no overlap in distribution and parameters
  if(length(intersect(names(distribution[[1]]), names(parameters))) > 1){
    stop("There should be no common parameters in parameters and distribution")
  }
  squire:::assert_list(parameters)
  squire:::assert_date(start_date)
  squire:::assert_int(rt_spacing)
  squire:::assert_int(n_particles)

  #calculate convience parameters
  data_start_date <- min(data$date_start)
  data_end_date <- max(data$date_end)
  #set dates so that start_date = 0
  data$t_start <- as.numeric(data$date_start - start_date)
  data$t_end <- as.numeric(data$date_end - start_date)

  #check that distribution and model are compatible?

  #check that rt_spacing and data are compatible
  if(!purrr::every(seq(data_start_date, data_end_date, by = rt_spacing), function(rt_date){
    #check a death period occurs in its coverage period
    (data %>%
     dplyr::filter(.data$date_start >= rt_date, .data$date_start == min(.data$date_start)) %>%
     dplyr::pull(.data$date_end) %>%
     `<=`(rt_date + rt_spacing)) &
      #check that there is no overlap on the low end
      (data %>%
       dplyr::filter(.data$date_end > rt_date, .data$date_end == min(.data$date_end)) %>%
       dplyr::pull(.data$date_start) %>%
       `>=`(rt_date)) &
      #check that there is no overlap on the high end
      (data %>%
       dplyr::filter(.data$date_start < rt_date + rt_spacing, .data$date_start == max(.data$date_start)) %>%
       dplyr::pull(.data$date_end) %>%
       `<=`(rt_date + rt_spacing))
  })){
    stop("data is not compatible with rt_spacing, every period of change in Rt must
         have a death period occuring within it and there can be no death periods
         occuring on the overlap of Rt change periods.")
  }

  #set up difference
  if(!is.null(dt) & !is.null(squire_model$odin_difference_model)){
    if(dt > 1){
      stop("dt (step size) must be smaller than 1")
    }
    use_difference <- TRUE
  } else {
    use_difference <- FALSE
    dt <- 1
  }

  map_func <- function(parameters, data, initial_infections_interval, n_particles, k, squire_model, rt_interval, p, use_difference, dt){
    #Determine the deaths used to fit each Rt
    temp <-  get_time_series(squire_model, parameters, data, rt_spacing)
    rt_df <- temp$rt_df
    split_data <- temp$split_data
    rm(temp)

    if(use_difference){
      #adjust parameter and data timings
      dt_multi <- 1/dt
      #parameters <- adjust_params_for_difference(squire_model, parameters, dt_multi)
      rt_df <- rt_df * dt_multi
      split_data <- purrr::map(split_data, function(x){
        x[, c("t_start", "t_end")] <-
          x[, c("t_start", "t_end")] * dt_multi
        x
      })
    }

    ##series of model specific functions are generated
    #generate function that just returns model output
    model <- generate_model_function(squire_model, parameters, use_difference, dt)
    #generate function that returns the deaths from this output
    deaths_function <- generate_deaths_function(model, dt)
    #basic likelihood function
    likelihood <- function(Rt, rt_index, initial_state){
      this_data <- split_data[[rt_index]]
      #get deaths for entire rt period until end of death period, this way index in this_data is the correct position in this time series
      cumulative_deaths <- deaths_function(Rt, rt_df$rt_change_t[rt_index], rt_df$death_end_t[rt_index], initial_state)
      deaths <- cumulative_deaths[this_data$index_end] - cumulative_deaths[this_data$index_start]
      #call generic likelihood function, add penalty for distance from current Rt value
      ll_negative_binomial(deaths, this_data$deaths, k)
    }
    #function to get the initial state value from a model
    get_initial_state <- function(Rt, rt_index, initial_state){
      #generate odin output from current Rt trend time to when the next Rt trend begins
      model_output <- model(Rt, rt_df$rt_change_t[rt_index], rt_df$initial_target[rt_index], initial_state) %>%
        #only keep the final value
        utils::tail(1)
      #convert this output into initial state values
      update_initial_state(initial_state, model_output)
    }
    R0_initial_state <- function(R0, initial_infections){
      #get the initial state as per country parameters
      initial_state <- setup_parameters(squire_model, parameters)
      get_initial_state(Rt = R0, rt_index = 1, initial_state =
                          assign_infections(initial_state, initial_infections)
      )
    }
    #calculate R0 + initial number of infections
    res <- particle_optimise_R0(initial_infections_interval, n_particles, likelihood, rt_interval, initial_state = setup_parameters(squire_model, parameters))
    Rt_trend <- res[1]
    initial_infections <- res[2]
    names(Rt_trend) <- "R0"
    #calculate the model state when the next rt trend starts
    initial_state <- R0_initial_state(res[1], res[2])

    #calculate each Rt trend onwards
    rt_trend <- purrr::reduce(
      #skip the first index as that is R0
      seq_along(rt_df$rt_change_t)[-1],
      function(state_trend, rt_index){
        #using state_trend$initial state, i.e. the state for the last Rt trend
        #we search over the range Rt values based on the previous value
        #(state_trend$rt_trend) and pick the highest likelihood
        res <- optimise_Rt(state_trend$initial_state, rt_index,
                           likelihood, rt_interval)
        #append our choice to the Rt_trend
        Rt_trend <- c(state_trend$Rt_trend, res)
        names(Rt_trend)[rt_index] <- paste0("Rt_", rt_index - 1)
        list(
          #get the new initial state
          initial_state = get_initial_state(res, rt_index, state_trend$initial_state),
          Rt_trend = Rt_trend
        )
        #and repeat process
      },
      .init = list(initial_state = initial_state, Rt_trend = Rt_trend)
    )$Rt_trend

    #calculate final model output for the whole time period covered by the data
    model_output <- model(rt_trend, t_start = 0,
                                   t_end = utils::tail(utils::tail(split_data, 1)[[1]]$t_end, 1),
                                   initial_state = assign_infections(setup_parameters(squire_model, parameters), initial_infections),
                                   tt_Rt = rt_df$rt_change_t,
                                   #run with higher tolerance, the odin model should never fail in the fitting though it can here
                                   atol = 10^-6, rtol = 10^-6)

    #diagnostics are null for now
    diagnostics <- NULL

    #update progressr so it knows that this sample is finished
    p()

    if(use_difference){
      #correct if using difference model
      rt_df$rt_change_t <- rt_df$rt_change_t * dt
      model_output[, "t"] <- model_output[, "t"] * dt
      model_output[, "time"] <- model_output[, "time"] * dt
    }

    #output values
    list(
      #our fitted Rt trend
      rt_trend = rt_trend,
      #the times relative to the start date that Rt changes
      rt_change_t = rt_df$rt_change_t,
      #our fitted initial number of infections
      initial_infections = initial_infections,
      model_output = model_output,
      diagnostics = diagnostics
    )
  }

  #setup progressr for user controlled feedback
  p <- progressr::progressor(steps = length(distribution), enable = TRUE)

  #Run the Fitting Process
  if(parallel & Sys.getenv("SQUIRE_PARALLEL_DEBUG") != "TRUE"){
    #set seed = NULL to suppress warnings due to squires (unused) random initial infections
    particle_output <- furrr::future_map(
      .x = purrr::map(distribution, ~append(.x, parameters)),
      .f = map_func,
      data = data,
      initial_infections_interval = initial_infections_interval,
      n_particles = n_particles, k = k,
      squire_model = squire_model, rt_interval = rt_interval,
      p = p, use_difference = use_difference, dt = dt,
      .options = furrr::furrr_options(seed = NULL)
    )
  } else {
    #otherwise use purrr for easier debugging
    particle_output <- purrr::map(
      .x = purrr::map(distribution, ~append(.x, parameters)),
      .f = map_func,
      data = data,
      initial_infections_interval = initial_infections_interval,
      n_particles = n_particles, k = k,
      squire_model = squire_model, rt_interval = rt_interval,
      p = p, use_difference = use_difference, dt = dt
    )
  }

  #remove null model objects where the solver has failed
  failed <- purrr::map_lgl(particle_output, ~is.null(.x$model_output))
  if(sum(failed) > 0){
    message(paste0(sum(failed), " draws failed to solve, removing from output, see $diagnostics$failed for the removed distribution entries"))
    particle_output <- particle_output[-which(failed)]
    failed_distributions <- distribution[which(failed)]
    distribution <- distribution[-which(failed)]
  } else {
    failed_distributions <- list()
  }
  #return model object
  #we must calculate some default user parameters with given initial infections
  #else these would be random
  default_user <- assign_infections(setup_parameters(squire_model, parameters), mean(initial_infections_interval))
  #get the model itself with basic parameters with a catch for different formats this takes
  if(inherits(squire_model$odin_model, "function")){
    odin_model <- squire_model$odin_model(user = default_user,
                                        unused_user_action = "ignore")
  } else {
    odin_model <- squire_model$odin_model$new(user = default_user,
                                        unused_user_action = "ignore")
  }
  return_obj <- list(
    parameters = parameters, #parameters fixed across each sample
    model = odin_model, #redundancy to keep compatible with squire/nimue functions
    squire_model = squire_model, #model object it self, kept for access to parameter function
    samples = purrr::map(seq_along(distribution), ~append( #add fitted Rt trends to distribution
      distribution[[.x]], list(
        R0 = particle_output[[.x]]$rt_trend,
        tt_R0 = particle_output[[.x]]$rt_change_t,
        initial_infections = particle_output[[.x]]$initial_infections
      )
    )),
    output = #merge simulation outputs into one
      abind::abind(purrr::map(particle_output, ~.x$model_output), along = 3, new.names = list(
        as.character(start_date + seq_len(nrow(particle_output[[1]]$model_output)) - 1),
        colnames(particle_output[[1]]$model_output),
        NULL
      )),
    inputs = list(
      #start_date and diagnostics
      start_date = start_date,
      data = data,
      rt_spacing = rt_spacing,
      k = k,
      n_particles = n_particles,
      initial_infections_interval = initial_infections_interval,
      rt_interval = rt_interval
    ),
    diagnostics = list(
      failed = failed_distributions
    )
  )
  #add a tag if vaccine model
  if(any(stringr::str_detect(class(squire_model), "nimue"))){
    class(return_obj) <- c("rt_optimised", "nimue_simulation")
  } else {
    class(return_obj) <- "rt_optimised"
  }
  return_obj
}
