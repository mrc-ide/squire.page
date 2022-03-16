#'
#'@export
excess_log_likelihood_vaccine <- function(pars, data, squire_model, model_params, pars_obs, n_particles,
                                  forecast_days = 0, return = "ll", Rt_args, interventions, ...) {
  squire:::assert_in(c("R0", "start_date"), names(pars), message = "Must specify R0, start date to infer")
  R0 <- pars[["R0"]]
  start_date <- pars[["start_date"]]
  squire:::assert_pos(R0)
  squire:::assert_date(start_date)
  date_Rt_change <- interventions$date_Rt_change
  date_contact_matrix_set_change <- interventions$date_contact_matrix_set_change
  date_ICU_bed_capacity_change <- interventions$date_ICU_bed_capacity_change
  date_hosp_bed_capacity_change <- interventions$date_hosp_bed_capacity_change
  date_vaccine_change <- interventions$date_vaccine_change
  date_vaccine_efficacy <- interventions$date_vaccine_efficacy

  #get efficacy from pars
  vaccine_efficacy <- vaccine_efficacies_to_efficacy_ts(
    pars[["ves"]], interventions$vaccine_efficacies, interventions$dose_ratio,
    date_vaccine_efficacy, pars_obs$delta_start_date,
    pars_obs$shift_duration, model_params$prob_hosp
  )

  if (is.null(date_Rt_change)) {
    tt_beta <- 0
  } else {
    #get the Rt values from R0 and the Rt_change values
    Rt <- evaluate_Rt_pmcmc_simple(R0 = R0, pars = pars)
    #get the dates in t and the corresponding Rt indexes
    tt_list <- squire:::intervention_dates_for_odin(dates = date_Rt_change,
                                                    change = seq(2, length(Rt)), start_date = start_date, steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)
    model_params$tt_beta <- tt_list$tt
    #reduce Rt to the values needed
    Rt <- Rt[tt_list$change]
  }
  if (is.null(date_contact_matrix_set_change)) {
    tt_contact_matrix <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_contact_matrix_set_change,
                                                    change = seq_along(interventions$contact_matrix_set)[-1],
                                                    start_date = start_date, steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)
    model_params$tt_matrix <- tt_list$tt
    model_params$mix_mat_set <- model_params$mix_mat_set[tt_list$change,, ]
  }
  if (is.null(date_ICU_bed_capacity_change)) {
    tt_ICU_beds <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_ICU_bed_capacity_change,
                                                    change = interventions$ICU_bed_capacity[-1], start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt), starting_change = interventions$ICU_bed_capacity[1])
    model_params$tt_ICU_beds <- tt_list$tt
    model_params$ICU_beds <- tt_list$change
  }
  if (is.null(date_hosp_bed_capacity_change)) {
    tt_hosp_beds <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_hosp_bed_capacity_change,
                                                    change = interventions$hosp_bed_capacity[-1], start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt), starting_change = interventions$hosp_bed_capacity[1])
    model_params$tt_hosp_beds <- tt_list$tt
    model_params$hosp_beds <- tt_list$change
  }
  if (is.null(date_vaccine_change)) {
    tt_vaccine <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_change,
                                                    change = interventions$max_vaccine[-1], start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt), starting_change = interventions$max_vaccine[1])
    model_params$tt_vaccine <- tt_list$tt
    model_params$max_vaccine <- tt_list$change
  }

  #vaccine taken from pars, estimated
  if(is.null(date_vaccine_efficacy)){
    #if not time varying
    model_params$vaccine_efficacy_infection <- vaccine_efficacy$vaccine_efficacy_infection
    model_params$tt_vaccine_efficacy_infection <- 0
    model_params$tt_vaccine_efficacy_disease <- 0
    model_params$prob_hosp <- vaccine_efficacy$vaccine_efficacy_disease
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy,
                                                    change = seq_along(date_vaccine_efficacy),
                                                    start_date = start_date, steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)

    model_params$tt_vaccine_efficacy_infection <- tt_list$tt
    model_params$vaccine_efficacy_infection <- vaccine_efficacy$vaccine_efficacy_infection[tt_list$change,
                                                                                       , ]
    model_params$tt_vaccine_efficacy_disease <- tt_list$tt
    model_params$prob_hosp <- vaccine_efficacy$vaccine_efficacy_disease[tt_list$change,
                                                     , ]
  }
  #calculate Beta from Rt
  beta_set <- squire:::beta_est(squire_model = squire_model, model_params = model_params,
                                R0 = Rt)
  model_params$beta_set <- beta_set
  #delta adjustments
  #get the dates in the shift as t
  shift_start <- as.integer(as.Date(pars_obs$delta_start_date) - start_date)
  shift_end <- as.integer(as.Date(pars_obs$delta_start_date) -
                            start_date +
                            pars_obs$shift_duration)
  #if the epidemic starts before the end of the shift we just swap over the numbers
  if(shift_end <= 0){
    model_params$prob_hosp_multiplier <- pars_obs$prob_hosp_multiplier
  } else {
    #we must figure where along we are and fit that in, slowly increase the
    #modified parameter until we reach the end of the shift

    #update prob_hosp
    tt_prob_hosp_multiplier <- seq(shift_start, shift_end, by = 1)
    prob_hosp_multiplier <- seq(model_params$prob_hosp_multiplier,
                                pars_obs$prob_hosp_multiplier,
                                length.out = length(tt_prob_hosp_multiplier))
    if(!(0 %in% tt_prob_hosp_multiplier)){
      #since we've already covered before the start we must be after and just
      #change the first entry to 0
      tt_prob_hosp_multiplier[1] <- 0
    }
    model_params$tt_prob_hosp_multiplier <- tt_prob_hosp_multiplier
    model_params$prob_hosp_multiplier <- prob_hosp_multiplier

    tt_dur_R <- c(shift_start, shift_end)
    gamma_R <- c(2/pars[["delta_dur_R"]], model_params$gamma_R)
    if(shift_start > 0){
      tt_dur_R <- c(0, tt_dur_R)
      gamma_R <- c(model_params$gamma_R, gamma_R)
    }
    model_params$tt_dur_R <- tt_dur_R
    model_params$gamma_R <- gamma_R
  }

  run_deterministic_comparison_excess(data = data,
                                                   squire_model = squire_model, model_params = model_params,
                                                   model_start_date = start_date, obs_params = pars_obs,
                                                   return = return)

}
#'
#' @noRd
vaccine_efficacies_to_efficacy_ts <- function(
  ves, vaccine_efficacies, dose_ratio,
  date_vaccine_efficacy, delta_start_date, shift_duration, prob_hosp
){
  #get efficacies scaled by ves
  effs <- scale_ves(ves, vaccine_efficacies)
  #generate delta adjust times on the same frame as dose ratios
  delta_prop <- get_delta_prop(delta_start_date, shift_duration, date_vaccine_efficacy)
  #add values for starting (shouldn't matter as there won't be vaccine when these are in effect)
  #just makes it consitent with other dates and variables
  dose_ratio <- c(0, dose_ratio)
  delta_prop <- c(utils::head(delta_prop, 1), delta_prop)
  #calculate efficacy at each change of dose and delta
  ve_i <- (effs$ve_i_low * (1 - dose_ratio) + dose_ratio * effs$ve_i_high) * delta_prop +
    (effs$ve_i_low_d * (1 - dose_ratio) + dose_ratio * effs$ve_i_high_d) * (1 - delta_prop)
  ve_d <- (effs$ve_d_low * (1 - dose_ratio) + dose_ratio * effs$ve_d_high) * delta_prop +
    (effs$ve_d_low_d * (1 - dose_ratio) + dose_ratio * effs$ve_d_high_d) * (1 - delta_prop)
  #scale for break through
  ve_d <- (ve_d - ve_i)/(1-ve_i)
  #convert into nimue formats
  vaccine_efficacy_infection <- nimue:::format_ve_i_for_odin(purrr::map(ve_i, ~.x), c(0, seq_along(date_vaccine_efficacy)))
  vaccine_efficacy_disease <- nimue:::format_ve_d_for_odin(purrr::map(ve_d, ~.x), c(0, seq_along(date_vaccine_efficacy)),
                                                           prob_hosp = prob_hosp[1,,1])
  #return
  list(
    vaccine_efficacy_infection = vaccine_efficacy_infection,
    vaccine_efficacy_disease = vaccine_efficacy_disease
  )
}
#'
#' @noRd
scale_ves <- function(ves, vaccine_efficacies){
  purrr::map(vaccine_efficacies, function(eff){
    eff[1] - (pos((0.5-ves))/0.5)*(eff[1] - eff[2]) + (pos((ves-0.5))/0.5)*(eff[3] - eff[1])
  })
}
#'
#' @noRd
pos <- function(x){
  if(x<0) {
    return(0)
  } else {
    return(x)
  }}
#'
#' @noRd
get_delta_prop <- function(delta_start_date, shift_duration, date_vaccine_efficacy){
  if(delta_start_date + shift_duration < min(date_vaccine_efficacy)){
    rep(1, length(date_vaccine_efficacy))
  } else if (delta_start_date > max(date_vaccine_efficacy)){
    rep(0, length(date_vaccine_efficacy))
  } else {
    delta_dates <- seq(delta_start_date, delta_start_date + shift_duration, by = 1)
    delta_prop <- seq(0, 1, length.out = length(delta_dates))
    #only keep in date_vaccine_efficacy
    c(
      rep(0, sum(date_vaccine_efficacy < delta_start_date)),
      delta_prop[delta_dates %in% date_vaccine_efficacy],
      rep(1, sum(date_vaccine_efficacy > delta_start_date + shift_duration))
    )
  }
}
#'
#' @noRd
evaluate_Rt_pmcmc_simple <- function(R0 = R0, pars = pars){
  as.numeric(c(R0, unlist(pars[grepl("Rt_", names(pars))])))
}

