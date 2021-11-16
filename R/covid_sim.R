#'
#'@export
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
    function(x){ out$interventions$vaccine_efficacy_infection[[x]][1] },
    numeric(1)
  )[-1]

  df$vaccine_efficacy_disease[vacc_pos] <- vapply(
    seq_along(out$interventions$vaccine_efficacy_disease),
    function(x){ out$interventions$vaccine_efficacy_disease[[x]][1] },
    numeric(1)
  )[-1]

  # and adjsut to be the reported efficacy rather than the breakthrough impact on disease after infection blocking
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
#'
#'@export
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
