#taken from squire, has added ability to update prob hospital multiplier, and
# variant dependant variables and duration in R. Also incorporates adjustment
# to cases for last few days
#'
#'@export
calc_loglikelihood_variant <- function(pars, data, squire_model, model_params,
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
  date_vaccine_efficacy_infection_change <- interventions$date_vaccine_efficacy_infection_change
  date_vaccine_efficacy_disease_change <- interventions$date_vaccine_efficacy_disease_change

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
                                                    change = interventions$max_vaccine[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = interventions$max_vaccine[1])
    model_params$tt_vaccine <- tt_list$tt
    model_params$max_vaccine <- tt_list$change
  }

  # and vaccine efficacy infection
  if (is.null(date_vaccine_efficacy_infection_change)) {
    tt_vaccine_efficacy_infection <- 0
  } else {

    # here we just pass the change as a position vector as we need to then
    # index the array of vaccine efficacies
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_infection_change,
                                                    change = seq_along(interventions$vaccine_efficacy_infection)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)

    model_params$tt_vaccine_efficacy_infection <- tt_list$tt

    # here we have to not index the array by the postion vectors that are reutrned by intervention_dates_for_odin
    model_params$vaccine_efficacy_infection <- model_params$vaccine_efficacy_infection[tt_list$change,,]
  }

  # and vaccine efficacy disease
  if (is.null(date_vaccine_efficacy_disease_change)) {
    tt_vaccine_efficacy_disease <- 0
  } else {

    # here we just pass the change as a position vector as we need to then
    # index the array of vaccine efficacies
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_disease_change,
                                                    change = seq_along(interventions$vaccine_efficacy_disease)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)

    model_params$tt_vaccine_efficacy_disease <- tt_list$tt

    # here we have to not index the array by the position vectors that are returned by intervention_dates_for_odin
    model_params$prob_hosp <- model_params$prob_hosp[tt_list$change,,]
  }

  #nimue specific functions (currently stored in pars_obs, may one day be interventions with the rest)
  if(!is.null(pars_obs$variant_adjust)){
    #check if we are using NIMUE
    if(!("nimue_model" %in% class(nimue::nimue_deterministic_model()))){
      stop("Can only specify variant adjustments for Nimue")
    }

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
  beta_set <- squire:::beta_est(squire_model = squire_model,
                                model_params = model_params,
                                R0 = R0)

  #----------------..
  # update the model params accordingly from new inputs
  #----------------..
  model_params$beta_set <- beta_set

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

#'
#'@noRd
run_deterministic_comparison_cases <- function(data,
                                               squire_model,
                                               model_params,
                                               model_start_date = "2020-02-02",
                                               obs_params = list(phi_cases = 0.1,
                                                                 k_cases = 2,
                                                                 phi_death = 1,
                                                                 k_death = 2,
                                                                 exp_noise = 1e6),
                                               forecast_days = 0,
                                               save_history = FALSE,
                                               return = "ll") {

  # parameter checks
  if (!(return %in% c("full", "ll", "sample", "single"))) {
    stop("return argument must be full, ll, sample", "single")
  }
  if (as.Date(data$date[data$deaths > 0][1], "%Y-%m-%d") < as.Date(model_start_date, "%Y-%m-%d")) {
    stop("Model start date is later than data start date")
  }

  # convert data into particle-filter form
  data <- squire:::particle_filter_data(data = data,
                                        start_date = model_start_date,
                                        steps_per_day = round(1 / model_params$dt))

  model_params$tt_beta <- round(model_params$tt_beta*model_params$dt)
  model_params$tt_contact_matrix <- round(model_params$tt_contact_matrix*model_params$dt)
  model_params$tt_hosp_beds <- round(model_params$tt_hosp_beds*model_params$dt)
  model_params$tt_ICU_beds <- round(model_params$tt_ICU_beds*model_params$dt)

  #set up model
  model_func <- squire_model$odin_model(user = model_params,
                                        unused_user_action = "ignore")

  # steps for the deterministic
  steps <- c(0, data$day_end)
  fore_steps <- seq(data$day_end[nrow(data)], length.out = forecast_days + 1L)
  steps <- unique(c(steps,fore_steps))

  # model run
  if("atol" %in% names(obs_params) && "rtol" %in% names(obs_params)) {
    squire:::assert_numeric(obs_params$atol)
    atol <- obs_params$atol
    squire:::assert_numeric(obs_params$rtol)
    rtol <- obs_params$rtol
  } else {
    atol <- 1e-6
    rtol <- 1e-6
  }
  #if full we use a low tolerance
  if(return == "full") {
    atol <- 1e-8
    rtol <- 1e-8
  }

  out <- tryCatch(
    model_func$run(t = seq(0, utils::tail(steps,1), 1), atol = atol, rtol = rtol),
    error = function(x){"FAIL"}
  )
  if(!identical(out, "FAIL")){

  index <- squire:::odin_index(model_func)

  # get deaths for comparison
  Ds <- diff(rowSums(out[,index$D]))
  Ds <- Ds[data$day_end[-1]]
  Ds[Ds < 0] <- 0
  deaths <- data$deaths[-1]

  # calculate ll for deaths
  if (obs_params$treated_deaths_only) {

    Ds_heathcare <- diff(rowSums(out[,index$D_get]))
    Ds_heathcare <- Ds_heathcare[data$day_end[-1]]
    ll <- squire:::ll_nbinom(deaths, Ds_heathcare, obs_params$phi_death, obs_params$k_death, obs_params$exp_noise)

  } else {

    ll <- squire:::ll_nbinom(deaths, Ds, obs_params$phi_death, obs_params$k_death, obs_params$exp_noise)

  }

  if(obs_params$cases_fitting){
    # calculate ll for the cases for last few days
    #calculate the relevant dates
    final_date <- max(data$date)
    cases_fitting_start_date <- final_date - obs_params$cases_days
    cases_reporting_start_date <- cases_fitting_start_date -
      obs_params$cases_reporting
    #get the relevant data
    data_reporting <- data %>%
      dplyr::filter(
        date >= cases_reporting_start_date & date < cases_fitting_start_date
      ) %>%
      dplyr::mutate(
        cases = dplyr::if_else(
          .data$cases < 0,
          0,
          .data$cases
        )
      )
    data_fitting <- data %>%
      dplyr::filter(
        date >= cases_fitting_start_date
      ) %>%
      dplyr::mutate(
        cases = dplyr::if_else(
          .data$cases < 0,
          0,
          .data$cases
        )
      )
    #get infections
    model_infections <- rowSums(out[,index$E2]) * model_params$gamma_E
    #ensure non negative
    model_infections[model_infections < 0] <- 0
    #check that cases are formatted correctly and exists
    if(all(is.na(data_reporting$cases)) |
       all(is.na(data_fitting$cases))){
      stop("Data for cases not formatted correctly, please fix or disable fitting to cases")
    }
    #estimate reporting fraction
    est_reporting_fraction <- mean(
      data_reporting$cases/model_infections[data_reporting$day_end],
      na.rm = TRUE
    )
    #cap how high this can be (prevents Infinities and is unlikely)
    est_reporting_fraction <- min(c(est_reporting_fraction, 1000))
    #note that though this can be larger than 1, then the model output is small,
    #this should correct as the model fits, this should still follow the shape of
    #the infections
    #calculate log likelihood with estimated reporting fraction and given variance
    cases_data <- stats::na.omit(cbind(
      data_fitting$cases,
      model_infections[data_fitting$day_end]
    ))
    llc <- squire:::ll_nbinom(
      cases_data[,1], cases_data[,2], est_reporting_fraction,
      obs_params$k_cases, obs_params$exp_noise
    )
    #note that we still fit to deaths in this time period
  } else {
    llc <- 0
  }

  # format the out object
  date <- data$date[[1]] + seq_len(nrow(out)) - 1L
  rownames(out) <- as.character(date)
  attr(out, "date") <- date

  # format similar to particle_filter nomenclature
  # allow full return for simulations
  if (return == "ll") {
    ret <- list(log_likelihood = sum(ll) + sum(llc),
                sample_state = out[nrow(out), ])
  } else if(return == "full") {
    ret <- out
  }
  } else {
    #if the model failed
    if (return == "ll") {
      ret <- list(log_likelihood = -.Machine$double.xmax,
                  sample_state = NULL)
    } else if(return == "full") {
      ret <- NULL
    }
  }
  ret
}
