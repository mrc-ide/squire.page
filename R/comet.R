#' Convert a model output into a data frame for Comet/GAVI Analysis
#'
#' Given a Nimue or Squire simulation object with MCMC fits, this function
#' produces a data frame of Rt over time and the equivalent Beta values, and if
#' an Nimue object, vaccine related parameters. It also extends this data frame
#' for the next 240 days for use in forecasting.
#'
#' Currently can only take fits to reported deaths and not excess mortality.
#'
#' @param out Nimue or squire model output, must have MCMC results in
#' \code{pmcmc_results} slot.
#' @param ox_interventions Data frame of government interventions, used to add
#' graphical elements, if \code{NULL} then element is not added to the data frame.
#' @return A \code{tibble} containing the necessary data to run Comet or other
#' analyses.
#' @examples
#' #apply this process to the standard Afghanistan fit
#' df <- prepare_input_json_df(
#'    afg_fit,
#'    ox_interventions = NULL
#'    )
#' #can be written to a json as needed
#' \dontrun{
#'    writeLines(
#'       jsonlite::toJSON(df, pretty = TRUE),
#'       "input_params.json"
#'       )
#' }
#' @export
prepare_input_json_df <- function(out, ox_interventions = NULL){
  if(any(!c("pmcmc_results", "replicate_parameters", "parameters") %in%
         names(out))){
    stop("out does not contain the needed data, pleasure insure this is the output of MCMC fit.")
  }
  iso3c <- countrycode::countrycode(out$parameters$country, origin = "country.name",
                                    destination = "iso3c")
  ## and save the info for the interface
  all_chains <- do.call(rbind,lapply(out$pmcmc_results$chains, "[[", "results"))
  if(is.null(all_chains)) {
    all_chains <- out$pmcmc_results$results
  }

  date_0 <- out$pmcmc_results$inputs$data$date[1]

  #extract parameters with highest posterior density
  best <- as.data.frame(all_chains[which.max(all_chains$log_posterior), ])
  #convert date to date
  best$start_date <- squire:::offset_to_start_date(date_0 ,round(best$start_date))
  #set the model's replicate to this so we can use the get_Rt function
  old_replicate_parameters <- out$replicate_parameters
  out$replicate_parameters <- best
  #get the Rt values
  df <- get_Rt(out) %>%
    dplyr::select(.data$date, .data$Rt, .data$t) %>%
    dplyr::mutate(#get beta values
      beta_set =
        squire:::beta_est(squire_model = out$pmcmc_results$inputs$squire_model,
                          model_params = out$pmcmc_results$inputs$model_params,
                          R0 = .data$Rt),
      grey_bar_start = FALSE
    ) %>%
    dplyr::rename(tt_beta = .data$t)


  ## -----------------------------------------------------------------------------

  #use replicates to calculate the uncertainty
  out$replicate_parameters <- old_replicate_parameters
  df <- dplyr::left_join(
    df,
    get_Rt(out) %>%
      dplyr::group_by(.data$date) %>%
      dplyr::summarise(
        Rt_min = stats::quantile(.data$Rt, 0.025,na.rm=TRUE),
        Rt_max = stats::quantile(.data$Rt, 0.975,na.rm=TRUE)
      ),
    by = "date") %>%
    tidyr::fill(tidyselect::all_of(c("Rt_min", "Rt_max")), .direction = "downup") %>%
    dplyr::mutate(#calculate beta values
      beta_set_min = squire:::beta_est(squire_model = out$pmcmc_results$inputs$squire_model,
                                       model_params = out$pmcmc_results$inputs$model_params,
                                       R0 = .data$Rt_min),
      beta_set_max = squire:::beta_est(squire_model = out$pmcmc_results$inputs$squire_model,
                                       model_params = out$pmcmc_results$inputs$model_params,
                                       R0 = .data$Rt_max)
    )

  ## -----------------------------------------------------------------------------

  # add in grey bar start for interface
  if(!is.null(ox_interventions)){
    ox_interventions_unique <- squire:::interventions_unique(ox_interventions[[iso3c]], "C")
    df$grey_bar_start[which.min(abs(as.numeric(df$date - ox_interventions_unique$dates_change[1])))] <- TRUE
  }

  # mark as unlikely fit because of death issues:
  df$recent_deaths <- TRUE
  if(sum(utils::tail(out$pmcmc_results$inputs$data$deaths, 20)) == 0) {
    df$recent_deaths <- FALSE
  }

  # add in the deaths to the json fits themselves
  df$deaths <- out$pmcmc_results$inputs$data$deaths[match(df$date, out$pmcmc_results$inputs$data$date)]
  #extend the df for a given period of time
  df <- extend_df_for_covidsim(df = df, out = out, ext = 240)
  df$iso3c <- iso3c

  #make adjustments for the vaccines if needed
  if(class(out) == "nimue_simulation"){
    if(is.null(out$interventions$vaccine_strategy)){
      stop("Vaccine strategy is not attached, if this is output of regular Squire functions manually add to out$interventions$vaccine_strategy as a string.")
    }
    df <- ammend_df_covidsim_for_vaccs(df,
                                       out,
                                       strategy =
                                         out$interventions$vaccine_strategy,
                                       iso3c = iso3c)
  }

  return(df)
}
#' Internal function to add vaccine data to covid-sim output
#' @noRd
ammend_df_covidsim_for_vaccs <- function(df, out, strategy, iso3c) {
  # add in vaccine pars
  df$max_vaccine <- 0
  df$vaccine_efficacy_infection <- out$interventions$vaccine_efficacy_infection[[1]][1]
  df$vaccine_efficacy_disease <- out$interventions$vaccine_efficacy_disease[[1]][1]
  vacc_pos <- which(as.Date(df$date) %in% as.Date(out$interventions$date_vaccine_change))
  df$max_vaccine[vacc_pos] <- as.integer(out$interventions$max_vaccine[-1][seq_along(vacc_pos)])
  df$max_vaccine[(max(vacc_pos)+1):length(df$max_vaccine)] <- as.integer(mean(utils::tail(df$max_vaccine[vacc_pos], 7)))

  df$vaccine_efficacy_infection[vacc_pos] <- vapply(
    seq_along(out$interventions$vaccine_efficacy_infection),
    function(x){out$interventions$vaccine_efficacy_infection[[x]][1]},
    numeric(1)
  )[-1]

  df$vaccine_efficacy_disease[vacc_pos] <- vapply(
    seq_along(out$interventions$vaccine_efficacy_disease),
    function(x){ out$interventions$vaccine_efficacy_disease[[x]][1] },
    numeric(1)
  )[-1]

  # and adjust to be the reported efficacy rather than the breakthrough impact on disease after infection blocking

  df$vaccine_efficacy_disease <- (df$vaccine_efficacy_disease * (1 - df$vaccine_efficacy_infection)) + df$vaccine_efficacy_infection
  df$vaccine_strategy <- strategy
  df$vaccine_coverage <- max(out$pmcmc_results$inputs$model_params$vaccine_coverage_mat)

  #set available_doses_proportion according to covax or no

  if(iso3c %in% get_covax_iso3c()){
    df$vaccines_available <- 0.2
  } else {
    df$vaccines_available <- 0.95
  }
  return(df)
}
#' Internal function for extending the covid-sim data frame over a certain
#' number of days
#'@noRd
extend_df_for_covidsim <- function(df, out, ext = 240) {

  betas <- df$beta_set
  tt_R0 <- df$tt_beta
  dates <- df$date

  # get an initial with the same seeds as before
  population <- squire::get_population(country = out$parameters$country, simple_SEIR = FALSE)
  init <- squire:::init_check_explicit(NULL, population$n, seeding_cases = 5)
  init$S <- init$S + init$E1
  init$E1 <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
  init$S <- init$S - c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0)

  # run model
  det_out <- squire::run_deterministic_SEIR_model(
    country = out$parameters$country,
    beta_set = betas,
    walker_params = FALSE,
    init = init,
    day_return = TRUE,
    tt_R0 = tt_R0+1,
    R0 = tt_R0,
    time_period = length(dates) + ext)

  # summarise the changes on the susceptibles
  index <- squire:::odin_index(det_out$model)

  # get the ratios
  mixing_matrix <- squire:::process_contact_matrix_scaled_age(
    out$pmcmc_results$inputs$model_params$contact_matrix_set[[1]],
    out$pmcmc_results$inputs$model_params$population
  )
  dur_ICase <- out$parameters$dur_ICase
  dur_IMild <- out$parameters$dur_IMild
  prob_hosp <- out$parameters$prob_hosp
  pop <- out$parameters$population

  prop_susc <- t(t(det_out$output[, index$S, 1])/pop)
  relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild

  adjusted_eigens <-  unlist(lapply(seq_len(nrow(prop_susc)), function(y) {
    if(any(is.na(prop_susc[y,]))) {
      return(NA)
    } else {
      Re(eigen(mixing_matrix*prop_susc[y,]*relative_R0_by_age)$values[1])
    }
  }))

  betas <- squire:::beta_est(squire_model = out$pmcmc_results$inputs$squire_model,
                             model_params = out$pmcmc_results$inputs$model_params,
                             R0 = df$Rt[1])

  ratios <- (betas * adjusted_eigens) / df$Rt[1]

  # extend df
  df2 <- df %>%
    tidyr::complete(tt_beta = seq(0, max(df$tt_beta) + ext, 1)) %>%
    dplyr::mutate(date = seq.Date(min(df$date), max(df$date)+ext,1)) %>%
    tidyr::fill(c("beta_set",4:10), .direction = "down") %>%
    dplyr::mutate(Reff = ratios*.data$Rt)

  return(df2)

}
