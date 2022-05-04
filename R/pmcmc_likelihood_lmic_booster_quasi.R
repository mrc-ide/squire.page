# Rewritten version of the booster likelihood to limit the range of data considered
# for each rt trend value. This is drjacoby specific, as it relies on some features
# from that package
#'
#'@export
calc_loglikelihood_quasi_booster <- function(params, data, misc){
  #if return is full then block is 1
  if(misc$return == "full"){
    misc$block <- 1
  }

  #----------------..
  # (potentially redundant) assertion
  #----------------..
  squire:::assert_in(c("R0", "start_date"), names(params),
                     message = "Must specify R0, start date to infer")

  #----------------..
  # unpack current params
  #----------------..
  R0 <- params[["R0"]]
  start_date <- squire:::numeric_to_start_date(misc$first_date,
                                               params[["start_date"]])

  # reporting fraction par if in params list
  if("rf" %in% names(params)) {
    squire:::assert_numeric(params[["rf"]])
    misc$pars_obs$phi_death <- params[["rf"]]
  }

  #----------------..
  # more assertions
  #----------------..
  squire:::assert_pos(R0)
  squire:::assert_date(start_date)

  #----------------..
  # setup model based on inputs and interventions
  #----------------..
  R0_change <- misc$interventions$R0_change
  date_R0_change <- misc$interventions$date_R0_change
  date_contact_matrix_set_change <- misc$interventions$date_contact_matrix_set_change
  date_ICU_bed_capacity_change <- misc$interventions$date_ICU_bed_capacity_change
  date_hosp_bed_capacity_change <- misc$interventions$date_hosp_bed_capacity_change
  date_vaccine_change <- misc$interventions$date_vaccine_change
  date_vaccine_efficacy_change <- misc$interventions$date_vaccine_efficacy_change

  # change betas
  if (is.null(date_R0_change)) {
    tt_beta <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_R0_change,
                                                    change = R0_change,
                                                    start_date = start_date,
                                                    steps_per_day = round(1/misc$model_params$dt),
                                                    starting_change = 1)
    misc$model_params$tt_beta <- tt_list$tt
    R0_change <- tt_list$change
    date_R0_change <- tt_list$dates
  }

  # and contact matrixes
  if (is.null(date_contact_matrix_set_change)) {
    tt_contact_matrix <- 0
  } else {

    # here just provide positions for change and then use these to index mix_mat_set
    tt_list <- squire:::intervention_dates_for_odin(dates = date_contact_matrix_set_change,
                                                    change = seq_along(misc$interventions$contact_matrix_set)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/misc$model_params$dt),
                                                    starting_change = 1)


    misc$model_params$tt_matrix <- tt_list$tt
    misc$model_params$mix_mat_set <- misc$model_params$mix_mat_set[tt_list$change,,]
  }

  # and icu beds
  if (is.null(date_ICU_bed_capacity_change)) {
    tt_ICU_beds <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_ICU_bed_capacity_change,
                                                    change = misc$interventions$ICU_bed_capacity[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/misc$model_params$dt),
                                                    starting_change = misc$interventions$ICU_bed_capacity[1])
    misc$model_params$tt_ICU_beds <- tt_list$tt
    misc$model_params$ICU_beds <- tt_list$change
  }

  # and hosp beds
  if (is.null(date_hosp_bed_capacity_change)) {
    tt_hosp_beds <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_hosp_bed_capacity_change,
                                                    change = misc$interventions$hosp_bed_capacity[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/misc$model_params$dt),
                                                    starting_change = misc$interventions$hosp_bed_capacity[1])
    misc$model_params$tt_hosp_beds <- tt_list$tt
    misc$model_params$hosp_beds <- tt_list$change
  }

  # and vaccine coverage
  if (is.null(date_vaccine_change)) {
    tt_vaccine <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_change,
                                                    change = seq_along(misc$interventions$first_doses)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/misc$model_params$dt),
                                                    starting_change = 1)
    misc$model_params$tt_first_doses <- tt_list$tt
    misc$model_params$first_doses <- misc$interventions$first_doses[tt_list$change]
    misc$model_params$tt_second_doses <- tt_list$tt
    misc$model_params$second_doses <- misc$interventions$second_doses[tt_list$change]
    misc$model_params$tt_booster_doses <- tt_list$tt
    misc$model_params$booster_doses <- misc$interventions$booster_doses[tt_list$change]
  }

  # and vaccine efficacy infection
  if (is.null(date_vaccine_efficacy_change)) {
    tt_vaccine_efficacy_infection <- 0
  } else {

    # here we just pass the change as a position vector as we need to then
    # index the array of vaccine efficacies
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_change,
                                                    change = seq_along(misc$interventions$vaccine_efficacy_infection)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/misc$model_params$dt),
                                                    starting_change = 1)

    misc$model_params$tt_vaccine_efficacy_infection <- tt_list$tt

    # here we have to not index the array by the postion vectors that are reutrned by intervention_dates_for_odin
    misc$model_params$vaccine_efficacy_infection <- misc$model_params$vaccine_efficacy_infection[tt_list$change,,]

    if(length(tt_list$tt) == 1){
      dim(misc$model_params$vaccine_efficacy_infection) <- c(1, dim(misc$model_params$vaccine_efficacy_infection))
    }
  }

  # and vaccine efficacy disease
  if (is.null(date_vaccine_efficacy_change)) {
    tt_vaccine_efficacy_disease <- 0
  } else {

    # here we just pass the change as a position vector as we need to then
    # index the array of vaccine efficacies
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_change,
                                                    change = seq_along(misc$interventions$vaccine_efficacy_disease)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/misc$model_params$dt),
                                                    starting_change = 1)

    misc$model_params$tt_vaccine_efficacy_disease <- tt_list$tt

    # here we have to not index the array by the position vectors that are returned by intervention_dates_for_odin
    misc$model_params$prob_hosp <- misc$model_params$prob_hosp[tt_list$change,,]

    if(length(tt_list$tt) == 1){
      dim(misc$model_params$prob_hosp) <- c(1, dim(misc$model_params$prob_hosp))
    }
  }

  #nimue specific functions (currently stored in pars_obs, may one day be interventions with the rest)
  if(!is.null(misc$pars_obs$variant_adjust)){
    #update dur_R if needed
    if(!is.null(misc$pars_obs$variant_adjust$gamma_R)){
      tt_list <- squire:::intervention_dates_for_odin(dates = misc$pars_obs$variant_adjust$date_dur_R_change,
                                                      change = seq_along(misc$pars_obs$variant_adjust$gamma_R)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_dur_R <- tt_list$tt
      model_params$gamma_R <- misc$pars_obs$variant_adjust$gamma_R[tt_list$change]
    }
    #update prob_hosp if needed
    if(!is.null(misc$pars_obs$variant_adjust$prob_hosp_multiplier)){
      tt_list <- squire:::intervention_dates_for_odin(dates = misc$pars_obs$variant_adjust$date_prob_hosp_multiplier_change,
                                                      change = seq_along(misc$pars_obs$variant_adjust$prob_hosp_multiplier)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_prob_hosp_multiplier <- tt_list$tt
      model_params$prob_hosp_multiplier <- misc$pars_obs$variant_adjust$prob_hosp_multiplier[tt_list$change]
    }
    #update prob_severe if needed
    if(!is.null(misc$pars_obs$variant_adjust$prob_severe_multiplier)){
      tt_list <- squire:::intervention_dates_for_odin(dates = misc$pars_obs$variant_adjust$date_prob_severe_multiplier_change,
                                                      change = seq_along(misc$pars_obs$variant_adjust$prob_severe_multiplier)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_prob_severe_multiplier <- tt_list$tt
      model_params$prob_severe_multiplier <- misc$pars_obs$variant_adjust$prob_severe_multiplier[tt_list$change]
    }
    #update dur_ICU if needed
    if(!is.null(misc$pars_obs$variant_adjust$gamma_get_mv_survive)){
      tt_list <- squire:::intervention_dates_for_odin(dates = misc$pars_obs$variant_adjust$date_dur_get_mv_survive_change,
                                                      change = seq_along(misc$pars_obs$variant_adjust$gamma_get_mv_survive)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_dur_get_mv_survive  <- tt_list$tt
      model_params$gamma_get_mv_survive  <- misc$pars_obs$variant_adjust$gamma_get_mv_survive[tt_list$change]
    }
    #update dur_ICU_death if needed
    if(!is.null(misc$pars_obs$variant_adjust$gamma_get_mv_die)){
      tt_list <- squire:::intervention_dates_for_odin(dates = misc$pars_obs$variant_adjust$date_dur_get_mv_die_change,
                                                      change = seq_along(misc$pars_obs$variant_adjust$gamma_get_mv_die)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_dur_get_mv_die  <- tt_list$tt
      model_params$gamma_get_mv_die  <-misc$pars_obs$variant_adjust$gamma_get_mv_die[tt_list$change]
    }
    #update dur_hosp if needed
    if(!is.null(misc$pars_obs$variant_adjust$gamma_get_ox_survive)){
      tt_list <- squire:::intervention_dates_for_odin(dates = misc$pars_obs$variant_adjust$date_dur_get_ox_survive_change,
                                                      change = seq_along(misc$pars_obs$variant_adjust$gamma_get_ox_survive)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_dur_get_ox_survive  <- tt_list$tt
      model_params$gamma_get_ox_survive  <-misc$pars_obs$variant_adjust$gamma_get_ox_survive[tt_list$change]
    }
    #update dur_hosp_death if needed
    if(!is.null(misc$pars_obs$variant_adjust$gamma_get_ox_die)){
      tt_list <- squire:::intervention_dates_for_odin(dates = misc$pars_obs$variant_adjust$date_dur_get_ox_die_change,
                                                      change = seq_along(misc$pars_obs$variant_adjust$gamma_get_ox_die)[-1],
                                                      start_date = start_date,
                                                      steps_per_day = round(1/model_params$dt),
                                                      starting_change = 1)
      model_params$tt_dur_get_ox_die <- tt_list$tt
      model_params$gamma_get_ox_die  <-misc$pars_obs$variant_adjust$gamma_get_ox_die[tt_list$change]
    }
  }

  #--------------------..
  # update new R0s based on R0_change and R0_date_change, and Meff_date_change
  #--------------------..
  # and now get new R0s for the R0
  R0 <- squire:::evaluate_Rt_pmcmc(R0_change = R0_change,
                                   R0 = R0,
                                   date_R0_change = date_R0_change,
                                   pars = as.list(params),
                                   Rt_args = misc$Rt_args)

  mod_class <- class(misc$squire_model)
  class(misc$squire_model) <-  c("nimue_model", "squire_model")
  # which allow us to work out our beta
  beta_set <- squire:::beta_est(squire_model = misc$squire_model,
                                model_params = misc$model_params,
                                R0 = R0)
  class(misc$squire_model) <- mod_class
  #----------------..
  # update the model params accordingly from new inputs
  #----------------..
  misc$model_params$beta_set <- beta_set

  #some parameters we won't use
  misc$pars_obs$treated_deaths_only <- FALSE

  #blocking
  if(misc$block > 1){
    dates_to_fit_over <- misc$misc_blocks[[misc$block - 1]]$dates
    if(misc$pars_obs$cases_fitting){
      misc$pars_obs$cases_fitting <- misc$misc_blocks[[misc$block - 1]]$last
    }
  } else {
    dates_to_fit_over <- NULL
  }

  #----------------..
  # run the deterministic comparison
  #----------------..
  pf_result <- run_deterministic_comparison_cases(
    data = data$data,
    squire_model = misc$squire_model,
    model_params = misc$model_params,
    model_start_date = start_date,
    obs_params = misc$pars_obs,
    forecast_days = misc$forecast_days,
    save_history = FALSE,
    return = misc$return,
    dates_to_fit_over
  )

  if(misc$return == "ll"){
    pf_result$log_likelihood
  } else if(misc$return == "full"){
    pf_result
  }

}
