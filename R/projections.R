#' Provide projections from calibrated simulations by changing RO, contact
#' matrices or bed availability.
#'
#' This extends previous \code{projections} as you can pass in lists of each argument
#' that then get passed to each simulation replicate.
#'
#' @details The user can specify changes to R0, contact matrices and bed
#' provision, which will come into effect from the current day in the calibration.
#' These changes can either set these to be specific values or change them
#' relative to their values in the original simulation. If no change is
#' requested, the simulation will use parameters chosen for the calibration run.
#' This extends previous versions of \code{projections} as you can now pass in
#' lists of each argument that then get passed to each simulation replicate.
#'
#' @param r Calibrated \code{{squire_simulation}} or \code{{rt_optimised}} object.
#' @param time_period How many days is the projection. Default = NULL, which will
#'   carry the projection forward from t = 0 in the calibration (i.e. the number
#'   of days set in calibrate using forecast)
#' @param R0 Numeric vector for R0 from t = 0 in the calibration.
#'   E.g. \code{R0 = c(2, 1)}. Default = NULL, which will use \code{R0_change}
#'   to alter R0 if provided.
#' @param R0_change Numeric vector for relative changes in R0 relative to the
#'   final R0 used in the calibration (i.e. at t = 0 in the calibration)
#'   E.g. \code{R0 = c(0.8, 0.5)}. Default = NULL, which will use \code{R0} to
#'   parameterise changes in R0 if provided.
#' @param tt_R0 Change time points for R0
#'
#' @param contact_matrix_set Contact matrices used in simulation. Default =
#'   NULL, which will use \code{contact_matrix_set_change} to alter the contact
#'   matrix if provided.
#' @param contact_matrix_set_change Numeric vector for relative changes in the
#'   contact matrix realtive to the final contact matrix used in the calibration
#'   (i.e. at t = 0 in the calibration).
#'   E.g. \code{contact_matrix_set_change = c(0.8, 0.5)}. Default = NULL, which
#'   will use \code{contact_matrix_set} to parameterise changes in contact
#'   matrices if if provided.
#' @param tt_contact_matrix Time change points for matrix change. Default = 0
#'
#' @param hosp_bed_capacity Numeric vector for hospital bed capacity
#'   from t = 0 in the calibration. Default = NULL, which will use
#'   \code{hosp_bed_capacity_change} to alter hosp_bed_capacity if provided.
#' @param hosp_bed_capacity_change Numeric vector for relative changes in
#'   hospital bed capacity relative to the final hospital bed capacity used in the
#'   calibration (i.e. at t = 0 in the calibration).
#'   E.g. \code{hosp_bed_capacity_change = c(0.8, 0.5)}. Default = NULL, which
#'   will use \code{hosp_bed_capacity} to parameterise changes in hospital bed capacity
#'   if provided.
#' @param tt_hosp_beds Change time points for hosp_bed_capacity
#'
#' @param ICU_bed_capacity Numeric vector for ICU bed capacity
#'   from t = 0 in the calibration. Default = NULL, which will use
#'   \code{ICU_bed_capacity_change} to alter ICU_bed_capacity if provided.
#' @param ICU_bed_capacity_change Numeric vector for relative changes in
#'   ICU bed capacity relative to the final ICU bed capacity used in the
#'   calibration (i.e. at t = 0 in the calibration).
#'   E.g. \code{ICU_bed_capacity_change = c(0.8, 0.5)}. Default = NULL, which
#'   will use \code{ICU_bed_capacity} to parameterise changes in ICU bed capacity
#'   if provided.
#' @param tt_ICU_beds Change time points for ICU_bed_capacity
#' @param to_be_run List of logicals for whether each replicate should be run.
#'   Default = TRUE, which causes all replictes to be run.
#' @param model_user_args List of other parameters to be passed to the model to
#'   be run. Default = NULL. An example would be:
#'
#'   \code{list(
#'   list(
#'   "prob_severe" = runif(17),
#'   "tt_dur_get_ox_survive" = c(0, 10),
#'   "gamma_get_ox_survive" = 0.2),
#'   list(
#'   "prob_severe" = runif(17),
#'   "tt_dur_get_mv_survive" = c(0, 5),
#'   "gamma_get_mv_survive" = 0.1)
#'   )}
#'
#'   The list should be the same length as the number of replicates in the
#'   simulations. Each list element should then be a list with elements
#'   named to match the arguments expected by the odin model with \code{r}.
#'   Above would be suitable to set the model parameters for a simulation with
#'   2 replicates. You do not have to have the same arguments in each list.
#' @param ... any other parameters, only used to pass these parameters to the methods
#'
#' @export
projections <- function(r, ...){
  UseMethod("projections")
}
#' @rdname projections
#' @export
projections.default <- function(r, ...){
  squire::projections
}
#' @rdname projections
#' @export
projections.rt_optimised <- function(
  r,
  time_period,
  R0 = NULL,
  R0_change = NULL,
  tt_R0 = NULL,
  contact_matrix_set = NULL,
  contact_matrix_set_change = NULL,
  tt_contact_matrix = NULL,
  hosp_bed_capacity = NULL,
  hosp_bed_capacity_change = NULL,
  tt_hosp_beds = NULL,
  ICU_bed_capacity = NULL,
  ICU_bed_capacity_change = NULL,
  tt_ICU_beds = NULL,
  to_be_run = TRUE,
  model_user_args = NULL,
  ...
){
  #get parameters for each sample
  parameters <- purrr::map(r$samples, function(sample){
    setup_parameters(r$squire_model, append(r$parameters, sample))
  })
  #convenience parameters
  final_t <- max(r$inputs$data$t_end)
  n_samples <- length(parameters)
  #format the variables into lists for each sample
  R0 <- convert_change_to_value(R0, R0_change, "R0", parameters, r$samples, final_t) %>%
    format_projection_variable(tt_R0, n_samples)
  contact_matrix_set <- convert_change_to_value(contact_matrix_set, contact_matrix_set_change, "contact_matrix_set", parameters, r$samples, final_t) %>%
    format_projection_variable(tt_contact_matrix, n_samples)
  hosp_bed_capacity <- convert_change_to_value(hosp_bed_capacity, hosp_bed_capacity_change, "hosp_beds", parameters, r$samples, final_t) %>%
    format_projection_variable(tt_hosp_beds, n_samples)
  ICU_bed_capacity <- convert_change_to_value(ICU_bed_capacity, ICU_bed_capacity_change, "ICU_beds", parameters, r$samples, final_t) %>%
    format_projection_variable(tt_ICU_beds, n_samples)
  #format to be run
  if(!(identical(to_be_run, TRUE) | identical(to_be_run, FALSE))){
    warning("For Rt Optimise outputs only generating all or no values is supported, proceding to generate all values")
    to_be_run <- TRUE
  }
  #check if model_args is long enough
  if(is.null(model_user_args)){
    model_user_args <- purrr::map(seq_along(parameters), ~list())
  } else if(length(model_user_args) < n_samples){
    stop("Ensure model_user_args has an entry for each sample")
  }
  #now setup these changes as changes to $samples
  r_output <- r
  r_output$samples <- purrr::map(seq_along(r$samples), function(sample_index){
    #gather together new parameters
    new_parameters <- append(model_user_args[[sample_index]],
                           list(
                              R0 = R0$par[[sample_index]],
                              tt_R0 = R0$tt[[sample_index]]
                           )) %>%
      append(
        list(
          contact_matrix_set = contact_matrix_set$par[[sample_index]],
          tt_contact_matrix = contact_matrix_set$tt[[sample_index]]
        )
      ) %>%
      append(
        list(
          hosp_bed_capacity = hosp_bed_capacity$par[[sample_index]],
          tt_hosp_beds = hosp_bed_capacity$tt[[sample_index]]
        )
      ) %>%
      append(
        list(
          ICU_bed_capacity = ICU_bed_capacity$par[[sample_index]],
          tt_ICU_beds = ICU_bed_capacity$tt[[sample_index]]
        )
      )
    #remove all null values
    new_parameters_clean <- list()
    for(par_index in seq_along(new_parameters)) {
      if(!is.null(new_parameters[[par_index]])){
        new_parameters_clean[[names(new_parameters)[par_index]]] <- new_parameters[[par_index]]
      }
    }
    rm(par_index)
    new_parameters <- new_parameters_clean
    rm(new_parameters_clean)
    #correct timings so that these changes occur on the final day
    tt_pars <- names(new_parameters)[stringr::str_detect(names(new_parameters), "tt_")]
    for(tt_par in tt_pars){
      new_parameters[[tt_par]] <- new_parameters[[tt_par]] + final_t
    }
    rm(tt_par)
    #append to samples
    pars_in_both <- intersect(names(new_parameters), names(r$samples[[sample_index]]))
    pars_in_samples <- setdiff(names(r$samples[[sample_index]]), pars_in_both)
    pars_in_new <- setdiff(names(new_parameters), pars_in_both)
    all_pars <- c(pars_in_both, pars_in_samples, pars_in_new)
    names(all_pars) <- all_pars
    purrr::map(
      all_pars, function(par){
        if(par %in% pars_in_both){
          #if its a time parameter, just add it on
          if(par %in% tt_pars | par == "R0"){
            c(r$samples[[sample_index]][[par]], new_parameters[[par]])
          } else {
            #append new value , making both a list if not already
            c(make_list(r$samples[[sample_index]][[par]]), make_list(new_parameters[[par]]))
          }
        } else if(par %in% pars_in_samples){
          r$samples[[sample_index]][[par]]
        } else if(par %in% pars_in_new){
          #if its a time par we add 0 for the baseline
          if(par %in% tt_pars){
            c(0, new_parameters[[par]])
          } else {
            new_parameters[[par]]
          }
        }
      }
    )
  })
  #now if asked generate outputs for each sample
  if(to_be_run){
    generate_draws(r_output, final_t + time_period, project_forwards = !is.null(r_output$output))
  } else {
    r_output
  }
}

#'@noRd
format_projection_variable <- function(x, tt_x, n_samples){
  if(is.null(x)){
    #we won't change it
    if(length(tt_x) > 1 | tt_x[1] > 0){
      message("Non NULL tt given for a NULL variable, assuming tt is NULL")
    }
    list(par = NULL, tt = NULL)
  } else {
    x_list <- is.list(x)
    tt_list <- is.list(tt_x)
    if(!x_list){
      x <-purrr::map(seq_len(n_samples), ~x)
    } else {
      if(length(x) != n_samples){
        stop("Ensure all variables have an entry for each sample")
      }
    }
    if(!tt_list){
      tt_x <-purrr::map(seq_len(n_samples), ~tt_x)
    } else {
      if(length(tt_x) != n_samples){
        stop("Ensure all tt have an entry for each sample")
      }
    }
    if(length(tt_x) != length(x)){
      stop("Ensure all tt and variable are the same length")
    }
    list(par = x, tt = tt_x)
  }
}

#'@noRd
convert_change_to_value <- function(main_par, change_par, par_name, parameters, samples, final_t){
  if(is.null(change_par)){
    if(!is.list(main_par)){
      main_par <- purrr::map(seq_along(samples), ~main_par)
    }
    return(main_par)
  } else {
    if(!is.list(change_par)){
      change_par <- purrr::map(seq_along(samples), ~change_par)
    }
    #special rules for some
    if(par_name == "R0"){
      return(purrr::map(seq_along(samples), function(sample_index){
        final_R0 <- block_interpolate(final_t, samples[[sample_index]]$R0,
                                      samples[[sample_index]]$tt_R0)
        #apply change
        as.numeric(final_R0 * change_par[[sample_index]])
      }))
    } else {
      return(purrr::map(seq_along(parameters), function(sample_index){
        final_par <- block_interpolate(final_t, parameters[[sample_index]][[par_name]],
                                      parameters[[sample_index]][[paste0("tt_", par_name)]])
        #apply change
        as.numeric(final_par * change_par[[sample_index]])
      }))
    }
  }
}

#' This can be made OOP, could handle vector t
#'@noRd
block_interpolate <- function(t, x, tt_x){
  indexes <- stats::approx(x = tt_x, y=seq_along(tt_x), xout = t, method = "constant", rule = 2)$y
  if(is.vector(x)){
    x[indexes]
  } else {
    if(is.list(x)){
      val <- x[indexes]
    } else if(is.array(x)){
      val <- x[indexes,,]
    } else if(is.matrix(x)){
      val <- x[indexes,]
    }
    if(length(indexes) == 1){
      val[[1]]
    } else {
      val
    }
  }
}

#'@noRd
make_list <- function(x){
  if(!is.list(x)){
    list(x)
  }
}
