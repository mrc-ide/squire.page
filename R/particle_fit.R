#' Fit nimue to a given set of death data via the particle exploration approach.
#'
#' For a given number of samples from given parameter uncertainty or distribution
#' function, this approach iteratively fits the Rt trend to the provided death
#' data and returns a nimue_simulation object for future usage in scenario
#' modelling.
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
#' @param rt_spacing Number of days between each Rt trend, default = 14 days.
#' @param initial_r The value of R0 to centre on initially, default = 3.
#' @param initial_infections The initial number of infections to centre on initially, default = 100.
#' @param proposal_width What proportion of the previous Rt trend to explore, default = 50%.
#' @param n_particles How many particles to explore uniformly across that range, default = 7.
#' @param k Control the dispersion on the negative binomial likelihood, default = 2.
#'
#' @return An object of type particle_fit, (model type).
#'
#' @export
particle_fit <- function(data, distribution, squire_model, parameters,
                         start_date,
                         rt_spacing = 14,
                         initial_r = 3,
                         initial_infections = 100,
                         proposal_width = 0.5,
                         n_particles = 14,
                         k = 2) {
  ##Initial Checks and Housekeeping

  #checks
  squire:::assert_dataframe(data)
  #sort data by start date
  data <- dplyr::arrange(data, .data$date_start) %>%
    dplyr::filter(.data$date_start > start_date) #must have occured after start date
  assert_particle_data_df(data)
  if(!("squire_model" %in% class(squire_model))){
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
  squire:::assert_numeric(proposal_width)
  squire:::assert_numeric(initial_r)

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

  #assume no deaths in this start period
  # data <- rbind(
  #   data.frame(deaths = 0, date_start = start_date, date_end = data$date_start[1],
  #              t_start = 0, t_end = data$t_start[1]),
  #   data
  # )


  ##Run the Fitting Process
  particle_output <- purrr::map(
    purrr::map(distribution, ~append(.x, parameters)),
    function(parameters, data, initial_r, initial_infections, proposal_width, n_particles, k, squire_model){
      temp <-  get_time_series(squire_model, parameters, data, rt_spacing)
      rt_df <- temp$rt_df
      split_data <- temp$split_data
      rm(temp)

      #generate function that just returns model output
      model <- generate_model_function(squire_model, parameters)
      #generate function that returns the deaths from this output
      deaths_function <- generate_deaths_function(model, squire_model, parameters)
      #basic likelihood function
      likelihood <- function(Rt, rt_index, initial_state){
        #format times so Rt_start_date is 0 and data
        this_data <- split_data[[rt_index]]
        this_data$t_end <-  this_data$t_end - rt_df$rt_change_t[rt_index]
        this_data$t_start <-  this_data$t_start - rt_df$rt_change_t[rt_index]
        death_start_t <- rt_df$death_start_t[rt_index] - rt_df$rt_change_t[rt_index]
        death_end_t <- rt_df$death_end_t[rt_index] - rt_df$rt_change_t[rt_index]

        cumulative_deaths <- deaths_function(Rt, death_start_t, death_end_t, initial_state)
        #calculate likelihood for each timeperiod then sum
        squire:::ll_nbinom(
          data = this_data$deaths,
          model = cumulative_deaths[this_data$t_end - death_start_t + 1] - cumulative_deaths[this_data$t_start - death_start_t + 1],
          phi = 1,
          k = k,
          exp_noise = Inf #no noise
        ) %>%
          sum()
      }
      #function to get the initial state value from a model
      get_initial_state <- function(Rt, rt_index, initial_state){
        #will probably need a model specific method
        #format dates so that rt_change_t = 0
        next_rt_t <- rt_df$initial_target[rt_index] - rt_df$rt_change_t[rt_index]
        model_output <- model(Rt, 0, next_rt_t, initial_state)
        update_initial_state(initial_state, model_output)
      }
      #R0 specific likelihood function
      R0_likelihood <- function(R0, initial_infections){
        #get the initial state as per country parameters
        initial_state <- do.call(squire_model$parameter_func, parameters)
        likelihood(Rt = R0, rt_index = 1, initial_state =
                     assign_infections(initial_state, initial_infections)
        )
      }
      # #R0 specific likelihood function
      # R0_likelihood <- function(R0, initial_infections){
      #   #get the initial state as per country parameters
      #   initial_state <- do.call(squire_model$parameter_func, parameters)
      #   cumulative_deaths <- deaths_function(R0, rt_df$death_start_t[1], rt_df$death_end_t[1], assign_infections(initial_state, initial_infections))
      #   #calculate likelihood for each timeperiod then sum
      #   this_data <- split_data[[1]]
      #   this_data$t_end <-  this_data$t_end - this_data$t_start[1] + 1
      #   this_data$t_start <-  this_data$t_start - this_data$t_start[1] + 1
      #   this_data %>%
      #     mutate(deaths = deaths/(t_end - t_start),
      #            t = t_start)
      # }
      R0_initial_state <- function(R0, initial_infections){
        #get the initial state as per country parameters
        initial_state <- do.call(squire_model$parameter_func, parameters)
        get_initial_state(Rt = R0, rt_index = 1, initial_state =
                     assign_infections(initial_state, initial_infections)
        )
      }
      #calculate R0 + initial number of infections, just optim for now
      res <- particle_optimise_R0(initial_r, initial_infections, n_particles, proposal_width, R0_likelihood)
      Rt_trend <- res[1]
      initial_infections <- res[2]
      names(Rt_trend) <- "R0"
      initial_state <- R0_initial_state(res[1], res[2])

      #calculate each Rt trend onwards
      rt_trend <- purrr::reduce(
        seq_along(rt_df$rt_change_t)[-1],
        function(state_trend, rt_index){
          res <- particle_optimise_Rt(rt = utils::tail(state_trend$Rt_trend, 1),
                               state_trend$initial_state, rt_index, n_particles,
                               proposal_width, likelihood)
          Rt_trend <- c(state_trend$Rt_trend, res)
          names(Rt_trend)[rt_index] <- paste0("Rt_", rt_index - 1)
          list(
            #get the initial state
            initial_state = get_initial_state(res, rt_index, state_trend$initial_state),
            Rt_trend = Rt_trend
          )
        },
        .init = list(initial_state = initial_state, Rt_trend = Rt_trend)
      )$Rt_trend

      #calculate final model output
      model_output <- model(rt_trend, t_start = 0,
                            t_end = utils::tail(utils::tail(split_data, 1)[[1]]$t_end, 1),
                            initial_state = assign_infections(do.call(squire_model$parameter_func, parameters), initial_infections),
                            tt_Rt = rt_df$rt_change_t)

      #diagnostics are null for now
      diagnostics <- NULL

      #output values
      list(
        rt_trend = rt_trend,
        rt_change_t = rt_df$rt_change_t,
        initial_infections = initial_infections,
        model_output = model_output,
        diagnostics = diagnostics
      )
    },
    data, initial_r, initial_infections, proposal_width, n_particles, k, squire_model
  )
  #return model object
  return_obj <- list(
    parameters = parameters,
    model = squire_model$odin_model$new(user = do.call(squire_model$parameter_func, parameters),
                                        unused_user_action = "ignore"), #redundancy to keep compatible with squire/nimue functions
    squire_model = squire_model,
    samples = purrr::map(seq_along(distribution), ~append( #add fitted Rt trends to distribution
      distribution[[.x]], list(
        R0 = particle_output[[.x]]$rt_trend,
        tt_R0 = particle_output[[.x]]$rt_change_t,
        initial_infections = particle_output[[.x]]$initial_infections
      )
    )),
    output = #merge simulation outputs into one
      purrr::reduce(particle_output, function(final_array, output){
        if(is.null(final_array)){
          n_arrays <- 1
        } else {
          n_arrays <- dim(final_array)[3] + 1
        }
        array(c(final_array, output$model_output), dim = c(dim(output$model_output), n_arrays),
              dimnames = list(
                as.character(start_date + seq_len(nrow(output$model_output)) - 1),
                colnames(output$model_output), NULL))
      }, .init = NULL),
    inputs = list(
      #start_date and diagnostics
      start_date = start_date,
      data = data
    )
  )
  if(any(stringr::str_detect(class(squire_model), "nimue"))){
    class(return_obj) <- c("particle_fit", "nimue_simulation")
  } else {
    class(return_obj) <- "particle_fit"
  }
  return_obj
}
