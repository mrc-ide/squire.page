#' Likelihood function for reported deaths with boosters
#' Not for individual use, only exported so it is easy to call in pmcmc functions.
#' @param pars Parameter inputs
#' @param data data to fit to
#' @param squire_model the model object
#' @param model_params the fixed parameters
#' @param pars_obs extra non fitted parameters
#' @param n_particles unused
#' @param forecast_days unused
#' @param return Specifies we are return likelihood for model outputs
#' @param Rt_args parameters for Rt
#' @param interventions parameters around vaccination
#' @param ... unused (for compatibility)
#' @export
calc_loglikelihood_booster <- function(pars, data, squire_model, model_params,
                                     pars_obs, n_particles,
                                     forecast_days = 0, return = "ll",
                                     Rt_args,
                                     interventions,
                                     ...) {
  #----------------..
  # (potentially redundant) assertion
  #----------------..
  squire:::assert_in(c("R0", "start_date"), names(pars),
                     message = "Must specify R0, start date to infer")

  #----------------..
  # unpack current params
  #----------------..
  R0 <- pars[["R0"]]
  start_date <- pars[["start_date"]]

  # reporting fraction par if in pars list
  if("rf" %in% names(pars)) {
    squire:::assert_numeric(pars[["rf"]])
    pars_obs$phi_death <- pars[["rf"]]
  }

  #----------------..
  # more assertions
  #----------------..
  squire:::assert_pos(R0)
  squire:::assert_date(start_date)

  #----------------..
  # setup model based on inputs and interventions
  #----------------..
  R0_change <- interventions$R0_change
  date_R0_change <- interventions$date_R0_change
  date_contact_matrix_set_change <- interventions$date_contact_matrix_set_change
  date_ICU_bed_capacity_change <- interventions$date_ICU_bed_capacity_change
  date_hosp_bed_capacity_change <- interventions$date_hosp_bed_capacity_change
  date_vaccine_change <- interventions$date_vaccine_change
  date_vaccine_efficacy_change <- interventions$date_vaccine_efficacy_change
  date_dur_V_change <- interventions$date_dur_V_change

  # change betas
  if (is.null(date_R0_change)) {
    tt_beta <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_R0_change,
                                                    change = R0_change,
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)
    model_params$tt_beta <- tt_list$tt
    R0_change <- tt_list$change
    date_R0_change <- tt_list$dates
  }

  # and contact matrixes
  if (is.null(date_contact_matrix_set_change)) {
    tt_contact_matrix <- 0
  } else {

    # here just provide positions for change and then use these to index mix_mat_set
    tt_list <- squire:::intervention_dates_for_odin(dates = date_contact_matrix_set_change,
                                                    change = seq_along(interventions$contact_matrix_set)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)


    model_params$tt_matrix <- tt_list$tt
    model_params$mix_mat_set <- model_params$mix_mat_set[tt_list$change,,]
  }

  # and icu beds
  if (is.null(date_ICU_bed_capacity_change)) {
    tt_ICU_beds <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_ICU_bed_capacity_change,
                                                    change = interventions$ICU_bed_capacity[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = interventions$ICU_bed_capacity[1])
    model_params$tt_ICU_beds <- tt_list$tt
    model_params$ICU_beds <- tt_list$change
  }

  # and hosp beds
  if (is.null(date_hosp_bed_capacity_change)) {
    tt_hosp_beds <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_hosp_bed_capacity_change,
                                                    change = interventions$hosp_bed_capacity[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = interventions$hosp_bed_capacity[1])
    model_params$tt_hosp_beds <- tt_list$tt
    model_params$hosp_beds <- tt_list$change
  }

  # and vaccine coverage
  if (is.null(date_vaccine_change)) {
    tt_vaccine <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_change,
                                                    change = seq_along(interventions$primary_doses)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)
    model_params$tt_primary_doses <- tt_list$tt
    model_params$primary_doses <- interventions$primary_doses[tt_list$change]
    model_params$tt_booster_doses <- tt_list$tt
    model_params$booster_doses <- interventions$booster_doses[tt_list$change]

    ##Delay dosing
    delayed <- apply_dose_delay_booster(seq(0, as.numeric(max(data$date) - start_date)),
                                        model_params$primary_doses, model_params$tt_primary_doses,
                                        model_params$booster_doses, model_params$tt_booster_doses,
                                        interventions$second_dose_delay, interventions$protection_delay_rate,
                                        interventions$protection_delay_shape)
    model_params$primary_doses <- delayed$primary_doses
    model_params$tt_primary_doses <- delayed$tt_primary_doses
    model_params$booster_doses <- delayed$booster_doses
    model_params$tt_booster_doses <- delayed$tt_booster_doses
    model_params$second_dose_delay <-delayed$second_dose_delay

  }

  # and vaccine efficacy infection
  if (is.null(date_vaccine_efficacy_change)) {
    tt_vaccine_efficacy_infection <- 0
  } else {

    # here we just pass the change as a position vector as we need to then
    # index the array of vaccine efficacies
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_change,
                                                    change = seq_along(interventions$vaccine_efficacy_infection)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)

    model_params$tt_vaccine_efficacy_infection <- tt_list$tt

    # here we have to not index the array by the postion vectors that are reutrned by intervention_dates_for_odin
    model_params$vaccine_efficacy_infection <- model_params$vaccine_efficacy_infection[tt_list$change,,]

    if(length(tt_list$tt) == 1){
      dim(model_params$vaccine_efficacy_infection) <- c(1, dim(model_params$vaccine_efficacy_infection))
    }
  }

  # and vaccine efficacy disease
  if (is.null(date_vaccine_efficacy_change)) {
    tt_vaccine_efficacy_disease <- 0
  } else {

    # here we just pass the change as a position vector as we need to then
    # index the array of vaccine efficacies
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_change,
                                                    change = seq_along(interventions$vaccine_efficacy_disease)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)

    model_params$tt_vaccine_efficacy_disease <- tt_list$tt

    # here we have to not index the array by the position vectors that are returned by intervention_dates_for_odin
    model_params$prob_hosp <- model_params$prob_hosp[tt_list$change,,]

    if(length(tt_list$tt) == 1){
      dim(model_params$prob_hosp) <- c(1, dim(model_params$prob_hosp))
    }
  }

  #and duration of vaccine protection
  if (is.null(date_dur_V_change)) {
    tt_dur_vaccine <- 0
  } else {

    # here we just pass the change as a position vector as we need to then
    # index the array of vaccine efficacies
    tt_list <- squire:::intervention_dates_for_odin(dates = date_dur_V_change,
                                                    change = seq_along(interventions$dur_V)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)

    model_params$tt_dur_vaccine <- tt_list$tt

    # here we have to not index the array by the position vectors that are returned by intervention_dates_for_odin
    model_params$gamma_vaccine <- model_params$gamma_vaccine[tt_list$change,]

    if(length(tt_list$tt) == 1){
      dim(model_params$gamma_vaccine) <- c(1, dim(model_params$gamma_vaccine))
    }
  }

  #nimue specific functions (currently stored in pars_obs, may one day be interventions with the rest)
  if(!is.null(pars_obs$variant_adjust)){
    #update dur_R if needed
    if(!is.null(pars_obs$variant_adjust$gamma_R)){
      tt_list <- squire:::intervention_dates_for_odin(dates = pars_obs$variant_adjust$date_dur_R_change,
                                                      change = seq_along(pars_obs$variant_adjust$gamma_R)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_dur_R <- tt_list$tt
      model_params$gamma_R <- pars_obs$variant_adjust$gamma_R[tt_list$change]
    }
    #update prob_hosp if needed
    if(!is.null(pars_obs$variant_adjust$prob_hosp_multiplier)){
      tt_list <- squire:::intervention_dates_for_odin(dates = pars_obs$variant_adjust$date_prob_hosp_multiplier_change,
                                                      change = seq_along(pars_obs$variant_adjust$prob_hosp_multiplier)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_prob_hosp_multiplier <- tt_list$tt
      model_params$prob_hosp_multiplier <- pars_obs$variant_adjust$prob_hosp_multiplier[tt_list$change]
    }
    #update prob_severe if needed
    if(!is.null(pars_obs$variant_adjust$prob_severe_multiplier)){
      tt_list <- squire:::intervention_dates_for_odin(dates = pars_obs$variant_adjust$date_prob_severe_multiplier_change,
                                                      change = seq_along(pars_obs$variant_adjust$prob_severe_multiplier)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_prob_severe_multiplier <- tt_list$tt
      model_params$prob_severe_multiplier <- pars_obs$variant_adjust$prob_severe_multiplier[tt_list$change]
    }
    #update dur_ICU if needed
    if(!is.null(pars_obs$variant_adjust$gamma_get_mv_survive)){
      tt_list <- squire:::intervention_dates_for_odin(dates = pars_obs$variant_adjust$date_dur_get_mv_survive_change,
                                                      change = seq_along(pars_obs$variant_adjust$gamma_get_mv_survive)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_dur_get_mv_survive  <- tt_list$tt
      model_params$gamma_get_mv_survive  <- pars_obs$variant_adjust$gamma_get_mv_survive[tt_list$change]
    }
    #update dur_ICU_death if needed
    if(!is.null(pars_obs$variant_adjust$gamma_get_mv_die)){
      tt_list <- squire:::intervention_dates_for_odin(dates = pars_obs$variant_adjust$date_dur_get_mv_die_change,
                                                      change = seq_along(pars_obs$variant_adjust$gamma_get_mv_die)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_dur_get_mv_die  <- tt_list$tt
      model_params$gamma_get_mv_die  <-pars_obs$variant_adjust$gamma_get_mv_die[tt_list$change]
    }
    #update dur_hosp if needed
    if(!is.null(pars_obs$variant_adjust$gamma_get_ox_survive)){
      tt_list <- squire:::intervention_dates_for_odin(dates = pars_obs$variant_adjust$date_dur_get_ox_survive_change,
                                                      change = seq_along(pars_obs$variant_adjust$gamma_get_ox_survive)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_dur_get_ox_survive  <- tt_list$tt
      model_params$gamma_get_ox_survive  <-pars_obs$variant_adjust$gamma_get_ox_survive[tt_list$change]
    }
    #update dur_hosp_death if needed
    if(!is.null(pars_obs$variant_adjust$gamma_get_ox_die)){
      tt_list <- squire:::intervention_dates_for_odin(dates = pars_obs$variant_adjust$date_dur_get_ox_die_change,
                                                      change = seq_along(pars_obs$variant_adjust$gamma_get_ox_die)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_dur_get_ox_die <- tt_list$tt
      model_params$gamma_get_ox_die  <-pars_obs$variant_adjust$gamma_get_ox_die[tt_list$change]
    }
  }

  #--------------------..
  # update new R0s based on R0_change and R0_date_change, and Meff_date_change
  #--------------------..
  # and now get new R0s for the R0
  R0 <- squire:::evaluate_Rt_pmcmc(R0_change = R0_change,
                                   R0 = R0,
                                   date_R0_change = date_R0_change,
                                   pars = pars,
                                   Rt_args = Rt_args)

  # which allow us to work out our beta
  #must adjust for changes to other tt's if neccessary
  if(length(c(model_params$tt_dur_ICase, model_params$tt_dur_IMild, model_params$tt_prob_hosp_multiplier)) > 3){
    tt_R0 <- unique(sort(c(model_params$tt_dur_ICase, model_params$tt_dur_IMild, model_params$tt_prob_hosp_multiplier, model_params$tt_beta)))
    R0 <- block_interpolate(tt_R0, R0, model_params$tt_beta)
  } else {
    tt_R0 <- model_params$tt_beta
  }
  beta_set <- beta_est(squire_model = squire_model,
                                model_params = model_params,
                                R0 = R0, tt_R0 = tt_R0)
  #----------------..
  # update the model params accordingly from new inputs
  #----------------..
  model_params$beta_set <- beta_set
  model_params$tt_beta <- tt_R0

  #some parameters we won't use
  pars_obs$treated_deaths_only <- FALSE

  #----------------..
  # run the deterministic comparison
  #----------------..
  pf_result <- run_deterministic_comparison_cases(
    data = data,
    squire_model = squire_model,
    model_params = model_params,
    model_start_date = start_date,
    obs_params = pars_obs,
    forecast_days = forecast_days,
    save_history = FALSE,
    return = return
  )
  # out
  pf_result
}
