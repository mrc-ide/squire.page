#' Methods to calculate Rt over time from each model type
#' @param model_out Model object
#'@export
get_Rt <- function(model_out){
  UseMethod("get_Rt")
}
#' @rdname get_Rt
#'@export
get_Rt.nimue_simulation <- function(model_out){
  date_0 <- max(as.Date(model_out$pmcmc_results$inputs$data$date))
  #get iso3c to attach
  iso3c <- squire::get_population(model_out$parameters$country)$iso3c[1]
  return(
    do.call(
      rbind,
      lapply(seq_along(model_out$replicate_parameters$R0), function(y) {
        tt <- squire:::intervention_dates_for_odin(dates = model_out$interventions$date_R0_change,
                                                   change = model_out$interventions$R0_change,
                                                   start_date = model_out$replicate_parameters$start_date[y],
                                                   steps_per_day = 1/model_out$parameters$dt)

        if(!("pmcmc_results" %in% names(model_out))) {
          Rt <- c(model_out$replicate_parameters$R0[y],
                  vapply(tt$change, model_out[["scan_results"]]$inputs$Rt_func, numeric(1),
                         R0 = model_out$replicate_parameters$R0[y], Meff = model_out$replicate_parameters$Meff[y]))
        } else {
          Rt <- squire:::evaluate_Rt_pmcmc(
            R0_change = tt$change,
            date_R0_change = tt$dates,
            R0 = model_out$replicate_parameters$R0[y],
            pars = as.list(model_out$replicate_parameters[y,]),
            Rt_args = model_out$pmcmc_results$inputs$Rt_args)
        }

        df <- data.frame(
          Rt = Rt,
          date = tt$dates
        ) %>%
          tidyr::complete(date = seq(min(.data$date), date_0, by = "days")) %>%
          dplyr::mutate(
            t = as.numeric(.data$date - min(.data$date))
          ) %>%
          tidyr::fill(.data$Rt) %>%
          dplyr::mutate(
            iso3c = iso3c,
            rep = y
          )
        return(df)
      })
    )
  )
}
#'
#' @rdname get_Rt
#'@export
get_Rt.excess_nimue_simulation <- function(model_out){
  date_0 <- max(as.Date(model_out$pmcmc_results$inputs$data$week_start))
  #get iso3c to attach
  iso3c <- squire::get_population(model_out$parameters$country)$iso3c[1]
  return(
    do.call(
      rbind,
      lapply(seq_along(model_out$replicate_parameters$R0), function(y) {

        #get the Rt values from R0 and the Rt_change values
        Rt <- evaluate_Rt_pmcmc_custom(R0 = model_out$replicate_parameters$R0[y],
                                       pars = as.list(model_out$replicate_parameters[y,]),
                                       Rt_args = model_out$pmcmc_results$inputs$Rt_args)
        #get the dates in t and the corresponding Rt indexes
        tt <- squire:::intervention_dates_for_odin(dates = model_out$interventions$date_Rt_change,
                                                   change = seq(2, length(Rt)),
                                                   start_date = model_out$replicate_parameters$start_date[y],
                                                   steps_per_day = 1/model_out$parameters$dt,
                                                   starting_change = 1)
        #reduce Rt to the values needed
        Rt <- Rt[tt$change]

        df <- data.frame(
          Rt = Rt,
          date = tt$dates
        ) %>%
          tidyr::complete(date = seq(min(.data$date), date_0, by = "days")) %>%
          dplyr::mutate(
            t = as.numeric(.data$date - min(.data$date))
          ) %>%
          tidyr::fill(.data$Rt) %>%
          dplyr::mutate(
            iso3c = iso3c,
            rep = y
          )
        return(df)
      })
    )
  )
}
#'
#' @rdname get_Rt
#'@export
get_Rt.vacc_durR_nimue_simulation <- function(model_out){
  date_0 <- max(as.Date(model_out$pmcmc_results$inputs$data$week_start))
  #get iso3c to attach
  iso3c <- squire::get_population(model_out$parameters$country)$iso3c[1]
  return(
    do.call(
      rbind,
      lapply(seq_along(model_out$replicate_parameters$R0), function(y) {

        #get the Rt values from R0 and the Rt_change values
        Rt <- evaluate_Rt_pmcmc_simple(R0 = model_out$replicate_parameters$R0[y],
                                       pars = as.list(model_out$replicate_parameters[y,]))
        #get the dates in t and the corresponding Rt indexes
        tt <- squire:::intervention_dates_for_odin(dates = model_out$interventions$date_Rt_change,
                                                   change = seq(2, length(Rt)),
                                                   start_date = model_out$replicate_parameters$start_date[y],
                                                   steps_per_day = 1/model_out$parameters$dt,
                                                   starting_change = 1)
        #reduce Rt to the values needed
        Rt <- Rt[tt$change]

        df <- data.frame(
          Rt = Rt,
          date = tt$dates
        ) %>%
          tidyr::complete(date = seq(min(.data$date), date_0, by = "days")) %>%
          dplyr::mutate(
            t = as.numeric(.data$date - min(.data$date))
          ) %>%
          tidyr::fill(.data$Rt) %>%
          dplyr::mutate(
            iso3c = iso3c,
            rep = y
          )
        return(df)
      })
    )
  )
}
#'
#' @rdname get_Rt
#'@export
get_Rt.rt_optimised <- function(model_out){
  date_0 <- model_out$inputs$start_date
  #get iso3c to attach
  iso3c <- squire::get_population(model_out$parameters$country)$iso3c[1]
  return(
    do.call(
      rbind,
      lapply(seq_along(model_out$samples), function(y) {

        #get the Rt values from R0 and the Rt_change values
        Rt <- model_out$samples[[y]]$R0
        #get the dates in t and the corresponding Rt indexes
        tt <- list(
          change = seq_along(Rt),
          dates = date_0 + model_out$samples[[y]]$tt_R0
        )
        #reduce Rt to the values needed

        df <- data.frame(
          Rt = Rt,
          date = date_0 + model_out$samples[[y]]$tt_R0
        ) %>%
          tidyr::complete(date = seq(date_0, max(.data$date), by = "days")) %>%
          dplyr::mutate(
            t = as.numeric(.data$date - min(.data$date))
          ) %>%
          tidyr::fill(.data$Rt) %>%
          dplyr::mutate(
            iso3c = iso3c,
            rep = y
          )
        return(df)
      })
    )
  )
}
