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
      t_end <- t_end - t_start
      t_start <- 0
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
generate_deaths_function <- function(model_func, parameters){
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
update_initial_state <- function(initial_state, model_output){
  model_output <- tail(model_output, 1)
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
