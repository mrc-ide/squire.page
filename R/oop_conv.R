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
#' @export
get_data_end_date.nimue_simulation <- function(model_out){
  max(model_out$pmcmc_results$inputs$data$date)
}
#' S3 Method to get final date in the data from a model
#' @export
get_data_end_date.excess_nimue_simulation <- function(model_out){
  max(model_out$pmcmc_results$inputs$data$week_end)
}
#' S3 Method to get final date in the data from a model
#' @export
get_data_end_date.particle_fit <- function(model_out){
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
#' @export
get_data_end_date_inner.nimue_simulation <- function(model_out){
  max(model_out$pmcmc_results$inputs$data$date)
}
#' S3 Method to get start date of the final time-period in the data from a model
#' @export
get_data_end_date_inner.excess_nimue_simulation <- function(model_out){
  max(model_out$pmcmc_results$inputs$data$week_start)
}
#' S3 Method to get start date of the final time-period in the data from a model
#' @export
get_data_end_date_inner.particle_fit <- function(model_out){
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
#' @export
get_data_start_date.nimue_simulation <- function(model_out){
  min(model_out$pmcmc_results$inputs$data$date)
}
#' S3 Method to get the start date of the data in a model
#' @export
get_data_start_date.excess_nimue_simulation <- function(model_out){
  min(model_out$pmcmc_results$inputs$data$week_start)
}
#' S3 Method to get the start date of the data in a model
#' @export
get_data_start_date.particle_fit <- function(model_out){
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
#' @export
get_dates.nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data$date
}
#' S3 Method to get all the start dates of time-periods in the data of a model
#' @export
get_dates.excess_nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data$week_start
}
#' S3 Method to get all the start dates of time-periods in the data of a model
#' @export
get_dates.particle_fit <- function(model_out){
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
#' @export
get_dates_greater.nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data$date
}
#' S3 Method to get all the end dates of time-periods in the data of a model
#' @export
get_dates_greater.excess_nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data$week_end
}
#' S3 Method to get all the end dates of time-periods in the data of a model
#' @export
get_dates_greater.particle_fit <- function(model_out){
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
#' @export
get_data.nimue_simulation <- function(model_out){
  model_out$pmcmc_results$inputs$data
}
#' S3 Method to get all the data of a model
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
#' @export
get_data.particle_fit <- function(model_out){
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
#' @param x particle_fit output object
#' @param replicates Plot replicates
#' @param summarise Logical, add summary line
#' @param ci logical add confidence interval ribbon
#' @param q Quantiles for upper and lower of interval ribbon
#' @param particle_fit For compatibility, if TRUE deaths are cumulative, else deaths are daily
#' @param ... placeholder for compatibility does nothing.
#'
#' @export
plot.particle_fit <- function(x, q = c(0.025, 0.975), replicates = TRUE, summarise = FALSE, ci = TRUE, particle_fit = FALSE, ...){
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

  if(ci){
    # p <- p +
    #   ggplot2::geom_ribbon(
    #     data = df %>%
    #       dplyr::group_by(.data$date, .data$compartment) %>%
    #       dplyr::summarise(ymin = quantile(.data$y, q[1]),
    #                        ymax = quantile(.data$y, q[2]),
    #                        .groups = "keep"),
    #     ggplot2::aes(x = .data$date, ymin = .data$ymin, ymax = .data$ymax,
    #                  fill = .data$compartment),
    #     alpha = 0.25, col = NA
    #   )
    p <- p +
      ggplot2::geom_line(ggplot2::aes(y=.data$ymin, x=.data$date), data = df %>%
                                 dplyr::group_by(.data$date, .data$compartment) %>%
                                 dplyr::summarise(ymin = quantile(.data$y, q[1]),
                                                  ymax = quantile(.data$y, q[2]),
                                                  .groups = "keep") , linetype="dashed") +
      ggplot2::geom_line(ggplot2::aes(y=.data$ymax, x=.data$date), data = df %>%
                           dplyr::group_by(.data$date, .data$compartment) %>%
                           dplyr::summarise(ymin = quantile(.data$y, q[1]),
                                            ymax = quantile(.data$y, q[2]),
                                            .groups = "keep"), linetype="dashed")
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
                           dplyr::summarise(y = median(.data$y), .groups = "keep"),
                         ggplot2::aes(x = .data$date,
                                      y = .data$y,
                                      col = .data$compartment)
                         )
  }
  p <- p +
    ggplot2::scale_color_discrete(name = "") +
    ggplot2::scale_fill_discrete(guide = "none") +
    ggplot2::xlab("Deaths") +
    ggplot2::ylab("Date") +
    ggplot2::theme_bw() +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 2))
}
#' S3 Generic to get total susceptible population
#' @noRd
get_total_s <- function(model_out){
  UseMethod("get_total_s")
}
#' S3 Method to get total susceptible population
#' @export
get_total_s.default <- function(model_out){
  sum(res$pmcmc_results$inputs$model_params$population)
}
#' S3 Method to get total susceptible population
#' @export
get_total_s.particle_fit <- function(model_out){
  if("population" %in% names(model_out$samples[[1]])){
    stop("population modified per sample, either calculate AR manuall or request this feature be added")
  }
  sum(do.call(model_out$squire_model$parameter_func, model_out$parameters)$population)
}
