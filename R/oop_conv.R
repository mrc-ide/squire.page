#' S3 Generic to get final date in the data from a model
#'
#' Smooths over the difference in data save locations in squire/nimue/squire.page
#' models.
#'
#' @noRd
get_data_end_date <- function(model_out){
  UseMethod("get_data_end_date")
}
#' S3 Method to get final date in the data from a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_data_end_date.nimue_simulation <- function(model_out){
  max(model_out$pmcmc_results$inputs$data$date)
}
#' S3 Method to get final date in the data from a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_data_end_date.excess_nimue_simulation <- function(model_out){
  max(model_out$pmcmc_results$inputs$data$week_end)
}
#' S3 Method to get final date in the data from a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_data_end_date.rt_optimised <- function(model_out){
  max(model_out$inputs$data$date_end)
}

#' S3 Generic to get start date of the final time-period in the data from a model
#'
#' Smooths over the difference in data save locations in squire/nimue/squire.page
#' models.
#'
#' @noRd
get_data_end_date_inner <- function(model_out){
  UseMethod("get_data_end_date_inner")
}
#' S3 Method to get start date of the final time-period in the data from a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_data_end_date_inner.nimue_simulation <- function(model_out){
  max(model_out$pmcmc_results$inputs$data$date)
}
#' S3 Method to get start date of the final time-period in the data from a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_data_end_date_inner.excess_nimue_simulation <- function(model_out){
  max(model_out$pmcmc_results$inputs$data$week_start)
}
#' S3 Method to get start date of the final time-period in the data from a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_data_end_date_inner.rt_optimised <- function(model_out){
  max(model_out$inputs$data$date_start)
}

#' S3 Generic to get the start date of the data in a model
#'
#' Smooths over the difference in data save locations in squire/nimue/squire.page
#' models.
#'
#' @noRd
get_data_start_date <- function(model_out){
  UseMethod("get_data_start_date")
}
#' S3 Method to get the start date of the data in a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_data_start_date.nimue_simulation <- function(model_out){
  min(model_out$pmcmc_results$inputs$data$date)
}
#' S3 Method to get the start date of the data in a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_data_start_date.excess_nimue_simulation <- function(model_out){
  min(model_out$pmcmc_results$inputs$data$week_start)
}
#' S3 Method to get the start date of the data in a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_data_start_date.rt_optimised <- function(model_out){
  min(model_out$inputs$data$date_start)
}

#' S3 Generic to get all the start dates of time-periods in the data of a model
#'
#' Smooths over the difference in data save locations in squire/nimue/squire.page
#' models.
#'
#' @noRd
get_dates <- function(model_out){
  UseMethod("get_dates")
}
#' S3 Method to get all the start dates of time-periods in the data of a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_dates.nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data$date
}
#' S3 Method to get all the start dates of time-periods in the data of a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_dates.excess_nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data$week_start
}
#' S3 Method to get all the start dates of time-periods in the data of a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_dates.rt_optimised <- function(model_out){
  model_out$inputs$data$date_start
}

#' S3 Generic to get all the end dates of time-periods in the data of a model
#'
#' Smooths over the difference in data save locations in squire/nimue/squire.page
#' models.
#'
#' @noRd
get_dates_greater <- function(model_out){
  UseMethod("get_dates_greater")
}
#' S3 Method to get all the end dates of time-periods in the data of a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_dates_greater.nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data$date
}
#' S3 Method to get all the end dates of time-periods in the data of a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_dates_greater.excess_nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data$week_end
}
#' S3 Method to get all the end dates of time-periods in the data of a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_dates_greater.rt_optimised <- function(model_out){
  model_out$inputs$data$date_end
}

#' S3 Generic to get all the data of a model
#'
#' Smooths over the difference in data save locations in squire/nimue/squire.page
#' models.
#'
#' @noRd
get_data <- function(model_out){
  UseMethod("get_data")
}
#' S3 Method to get all the data of a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_data.nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data
}
#' S3 Method to get all the data of a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_data.excess_nimue_simulation <- function(model_out){
  data <- model_out$pmcmc_results$inputs$data
  #make consistent with particle fits, not worth changing in functions
  data$date_start <- data$week_start
  data$date_end <- data$week_end
  data$week_end <- data$week_start <- NULL
  data
}
#' S3 Method to get all the data of a model
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_data.rt_optimised <- function(model_out){
  model_out$inputs$data
}

#' Add some adjustments for excess data in plotting nimue simulations.
#'
#' @param x An excess_nimue_simulation object
#' @param replicates Plot replicates
#' @param summarise Logical, add summary line
#' @param ci logical add confidence interval ribbon
#' @param q Quantiles for upper and lower of interval ribbon
#' @param var_select Vector of variable names to plot (default is all)
#' @param summary_f Function to summarise each compartment
#'   passed to the \code{fun} argument of \code{\link[ggplot2]{stat_summary}}
#' @param x_var X variable to use for plotting (default is \code{"t"},
#'   but can be set to, \code{"date"}, if \code{date_0} provided), which will
#'   cause the date to be plotted rather than time.
#' @param particle_fit If the squire_simulation provided is the result of
#'   running the particle filter, do we want to just plot the fit. Default =
#'   FALSE
#' @param date_0 Date of time 0 (e.g. "2020-03-01"), if specified a date column
#'   will be added
#' @param ... additional arguments affecting the plot produced.
#'
#' @export
plot.excess_nimue_simulation <- function(x, var_select = NULL, replicates = FALSE, summarise = TRUE,
                                         ci = TRUE, q = c(0.025, 0.975), summary_f = mean, x_var = "t",
                                         date_0 = NULL, particle_fit = FALSE, ...) {
  #set the dates correctly
  x$pmcmc_results$inputs$data$date <- x$pmcmc_results$inputs$data$week_start
  #call nimue function
  nimue:::plot.nimue_simulation(
    x = x, var_select = var_select, replicates = replicates,
    summarise = summarise, ci = ci, q = q, summary_f = summary_f, x_var = x_var,
    date_0 = date_0, particle_fit = particle_fit, ...
  )
}

#' Essentially just does squire::plot_pmcmc_sample. Not intended for use otherwise.
#'
#' @param x rt_optimised output object
#' @param replicates Plot replicates
#' @param summarise Logical, add summary line
#' @param ci logical add confidence interval ribbon
#' @param q Quantiles for upper and lower of interval ribbon
#' @param particle_fit For compatibility, if TRUE deaths are cumulative, else deaths are daily
#' @param ... placeholder for compatibility does nothing.
#'
#' @export
plot.rt_optimised <- function(x, q = c(0.025, 0.975), replicates = TRUE, summarise = FALSE, ci = TRUE, particle_fit = FALSE, ...){
  # #set up our parameters
  # x$parameters$day_return <- TRUE
  # x$parameters$replicates <- length(x$samples)
  # x$parameters$dt <- 1
  # x$pmcmc_results <- list(
  #   inputs = list(
  #     data = x$inputs$data
  #   )
  # )
  # x$pmcmc_results$inputs$data$date <- get_dates(x)
  # if("nimue_simulation" %in% class(x)){
  #   class(x) <- "nimue_simulation"
  #   nimue:::plot.nimue_simulation(x, ...)
  # } else {
  #   #
  #   class(x$squire_model)
  #   squire:::plot
  #   class(x) <- "squire_simulation"
  #   squire:::plot.squire_simulation(x, particle_fit = TRUE)
  # }
  #keep this simple and just return the deaths (squire doesn't seem capable of extracting deaths for some reason)
  df <- nimue_format(x, "D", date_0 = x$inputs$start_date) %>%
    dplyr::group_by(.data$replicate)
  if(particle_fit){
    df <- dplyr::mutate(df, y = diff(c(0, .data$y)))
  }
  p <- ggplot2::ggplot()

  if(ci & particle_fit){
    p <- p +
      ggplot2::geom_line(ggplot2::aes(y=.data$ymin, x=.data$date), data = df %>%
                                 dplyr::group_by(.data$date, .data$compartment) %>%
                                 dplyr::summarise(ymin = stats::quantile(.data$y, q[1]),
                                                  ymax = stats::quantile(.data$y, q[2]),
                                                  .groups = "keep") , linetype="dashed") +
      ggplot2::geom_line(ggplot2::aes(y=.data$ymax, x=.data$date), data = df %>%
                           dplyr::group_by(.data$date, .data$compartment) %>%
                           dplyr::summarise(ymin = stats::quantile(.data$y, q[1]),
                                            ymax = stats::quantile(.data$y, q[2]),
                                            .groups = "keep"), linetype="dashed")
  }

  if(ci & !particle_fit){
    p <- p +
      ggplot2::geom_ribbon(
        data = df %>%
          dplyr::group_by(.data$date, .data$compartment) %>%
          dplyr::summarise(ymin = stats::quantile(.data$y, q[1]),
                           ymax =  stats::quantile(.data$y, q[2]),
                           .groups = "keep"),
        ggplot2::aes(x = .data$date, ymin = .data$ymin, ymax = .data$ymax,
                     fill = .data$compartment),
        alpha = 0.25, col = NA
      )
  }

  if(replicates){
    p <- p +
      ggplot2::geom_line(data = df,
                         ggplot2::aes(x = .data$date,
                                      y = .data$y,
                                      col = .data$compartment,
                                      group = interaction(.data$compartment, .data$replicate)),
                         alpha = max(0.2, 1 / length(unique(df$replicate)))
                         )
  }
  if(summarise){
    p <- p +
      ggplot2::geom_line(data = df %>%
                           dplyr::group_by(.data$compartment, .data$date) %>%
                           dplyr::summarise(y = stats::median(.data$y), .groups = "keep"),
                         ggplot2::aes(x = .data$date,
                                      y = .data$y,
                                      col = .data$compartment)
                         )
  }
  p <- p +
    ggplot2::scale_color_discrete(name = "") +
    ggplot2::scale_fill_discrete(guide = "none") +
    ggplot2::xlab("Date") +
    ggplot2::ylab("Deaths") +
    ggplot2::theme_bw() +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 2))
}
#' Calls plot.rt_optimised and adds excluded fits . Not intended for use otherwise.
#'
#' @param x rt_optimised output object
#' @param replicates Plot replicates
#' @param summarise Logical, add summary line
#' @param ci logical add confidence interval ribbon
#' @param q Quantiles for upper and lower of interval ribbon
#' @param particle_fit For compatibility, if TRUE deaths are cumulative, else deaths are daily
#' @param ... placeholder for compatibility does nothing.
#'
#' @export
plot.rt_optimised_trimmed <- function(x, q = c(0.025, 0.975), replicates = TRUE, summarise = FALSE, ci = TRUE, particle_fit = FALSE, ...){
  gg <- plot.rt_optimised(x = x, q = q, replicates = replicates, summarise = summarise, ci = ci, partcile_fit = particle_fit, ...)
  #get the excluded trajcetories
  x$output <- x$excluded$output
  df <- nimue_format(x, "D", date_0 = x$inputs$start_date) %>%
    dplyr::group_by(.data$replicate)
  if(particle_fit){
    df <- dplyr::mutate(df, y = diff(c(0, .data$y)))
  }
  gg +
    ggplot2::geom_line(df, ggplot2::aes(x = .data$date, y = .data$y), colour = "yellow", line.type = "dashed")
}
#' S3 Generic to get total susceptible population
#' @noRd
get_parameters <- function(model_out){
  UseMethod("get_parameters")
}
#' S3 Method to get total susceptible population
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_parameters.default <- function(model_out){
  model_out$pmcmc_results$inputs$model_params
}
#' S3 Method to get total susceptible population
#' @param model_out A nimue/squire mcmc or particle fit output
#' @export
get_parameters.rt_optimised <- function(model_out){
  if("population" %in% names(model_out$samples[[1]])){
    stop("population modified per sample, either calculate AR manually or request this feature be added")
  }
  #ensure these values are not here
  model_out$parameters$day_return <- NULL
  model_out$parameters$replicates <- NULL
  setup_parameters(model_out$squire_model, model_out$parameters)
}

#' An S3 Generic to get the parameters of a given model
#' @noRd
setup_parameters <- function(model_obj, parameters){
  UseMethod("setup_parameters")
}
#' If no particular method we default to calling the parameter function
#' attached to the model
#' @param model_obj A nimue/squire model object
#' @param parameters The parameters to the pass to the parameter function
#' @export
setup_parameters.default <- function(model_obj, parameters){
  do.call(model_obj$parameter_func, parameters)
}
#' If no particular method we default to calling the parameter function
#' attached to the model
#' @param model_obj A nimue/squire model object
#' @param parameters The parameters to the pass to the parameter function
#' @export
setup_parameters.nimue_model <- function(model_obj, parameters){
  #fill population if missing
  if(!"population" %in% names(parameters)){
    parameters$population <- squire::get_population(parameters$country)$n
  }
  #fill contact matrix if missing
  if(!"contact_matrix_set" %in% names(parameters)){
    parameters$contact_matrix_set <-
      squire::get_mixing_matrix(parameters$country)
  }
  #fill tt_contact_matrix if missing
  if(!"tt_contact_matrix" %in% names(parameters)){
    if(length(dim(parameters$contact_matrix_set)) == 2){
      parameters$tt_contact_matrix <- 0
    } else {
      parameters$tt_contact_matrix <- dim(parameters$contact_matrix_set)[1]
    }
  }
  #fill tt_contact_matrix if missing
  if(!"tt_contact_matrix" %in% names(parameters)){
    if(length(dim(parameters$contact_matrix_set)) == 2){
      parameters$tt_contact_matrix <- 0
    } else {
      parameters$tt_contact_matrix <- dim(parameters$contact_matrix_set)[1]
    }
  }
  #fill hosp_bed_capacity if missing
  if(!"hosp_bed_capacity" %in% names(parameters)){
    parameters$hosp_bed_capacity <- squire::get_healthcare_capacity(parameters$country)$hosp_beds * sum(parameters$population)/1000
  }
  #fill tt_hosp_beds if missing
  if(!"tt_hosp_beds" %in% names(parameters)){
    parameters$tt_hosp_beds <- seq_along(parameters$hosp_bed_capacity) - 1
  }
  #fill ICU beds if missing
  if(!"ICU_bed_capacity" %in% names(parameters)){
    parameters$ICU_bed_capacity <- squire::get_healthcare_capacity(parameters$country)$ICU_beds * sum(parameters$population)/1000
  }
  #fill tt_ICU_beds if missing
  if(!"tt_ICU_beds" %in% names(parameters)){
    parameters$tt_ICU_beds <- seq_along(parameters$ICU_bed_capacity) - 1
  }
  #fill dt if missing
  if(!"dt" %in% names(parameters)){
    parameters$dt <- 1
  }

  do.call(model_obj$parameter_func, parameters)
}
#' An S3 generic for estimating Beta
#' @noRd
beta_est <- function(squire_model, model_params, R0) {
  UseMethod("beta_est")
}
#' An S3 method for getting Beta for a model
#' @param squire_model A model object
#' @param model_params Parameters for the model
#' @param R0 R0/Rt value to translate into beta
#' @export
beta_est.default <- function(squire_model, model_params, R0) {
  squire:::beta_est(squire_model, model_params, R0)
}
#' An S3 method for getting Beta for a model
#' @param squire_model A model object
#' @param model_params Parameters for the model
#' @param R0 R0/Rt value to translate into beta
#' @export
beta_est.booster_model <- function(squire_model, model_params, R0) {
  #treat this as nimue
  class(squire_model) <- c("nimue_model", "squire_model")
  squire:::beta_est(squire_model, model_params, R0)
}
