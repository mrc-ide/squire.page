#' Likelihood function for excess mortality data
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
excess_log_likelihood <- function(pars, data, squire_model, model_params, pars_obs, n_particles,
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
  date_vaccine_efficacy_infection_change <- interventions$date_vaccine_efficacy_infection_change
  date_vaccine_efficacy_disease_change <- interventions$date_vaccine_efficacy_disease_change
  if (is.null(date_Rt_change)) {
    tt_beta <- 0
  } else {
    #get the Rt values from R0 and the Rt_change values
    Rt <- evaluate_Rt_pmcmc_custom(R0 = R0, pars = pars, Rt_args = Rt_args)
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
  if (is.null(date_vaccine_efficacy_infection_change)) {
    tt_vaccine_efficacy_infection <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_infection_change,
                                                    change = seq_along(interventions$vaccine_efficacy_infection)[-1],
                                                    start_date = start_date, steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)
    model_params$tt_vaccine_efficacy_infection <- tt_list$tt
    model_params$vaccine_efficacy_infection <- model_params$vaccine_efficacy_infection[tt_list$change,
                                                                                       , ]
  }
  if (is.null(date_vaccine_efficacy_disease_change)) {
    tt_vaccine_efficacy_disease <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_disease_change,
                                                    change = seq_along(interventions$vaccine_efficacy_disease)[-1],
                                                    start_date = start_date, steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)
    model_params$tt_vaccine_efficacy_disease <- tt_list$tt
    model_params$prob_hosp <- model_params$prob_hosp[tt_list$change,
                                                     , ]
  }
  #calculate Beta from Rt
  beta_set <- squire:::beta_est(squire_model = squire_model, model_params = model_params,
                                R0 = Rt)
  model_params$beta_set <- beta_set

  #make the delta adjustments
  if("prob_hosp_multiplier" %in% names(pars_obs) |
     "dur_R" %in% names(pars_obs)) {
    #get the dates in the shift as t
    shift_start <- as.integer(as.Date(pars_obs$delta_start_date) - start_date)
    shift_end <- as.integer(as.Date(pars_obs$delta_start_date) -
                              start_date +
                              pars_obs$shift_duration)
    #if the epidemic starts before the end of the shift we just swap over the numbers
    if(shift_end <= 0 & "prob_hosp_multiplier" %in% names(pars_obs)){
      model_params$prob_hosp_multiplier <- pars_obs$prob_hosp_multiplier
    } else {
      #we must figure where along we are and fit that in, slowly increase the
      #modified parameter until we reach the end of the shift
      if("prob_hosp_multiplier" %in% names(pars_obs) & (
        pars_obs$prob_hosp_multiplier != model_params$prob_hosp_multiplier
        | is.null(model_params$prob_hosp_multiplier)
      )){
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
      }
      if("dur_R" %in% names(pars_obs) & (
        2/pars_obs$dur_R != model_params$gamma_R
      )){
        tt_dur_R <- c(shift_start, shift_end)
        gamma_R <- c(2/pars_obs$dur_R, model_params$gamma_R)
        if(shift_start > 0){
          tt_dur_R <- c(0, tt_dur_R)
          gamma_R <- c(model_params$gamma_R, gamma_R)
        }
        model_params$tt_dur_R <- tt_dur_R
        model_params$gamma_R <- gamma_R
      }
    }
  }
  run_deterministic_comparison_excess(data = data,
                                                     squire_model = squire_model, model_params = model_params,
                                                     model_start_date = start_date, obs_params = pars_obs,
                                                     return = return)

}


#'
#'@noRd
run_deterministic_comparison_excess <- function(data, squire_model, model_params, model_start_date = "2020-02-02",
                                                obs_params = list(
                                                  phi_cases = 0.1,
                                                  k_cases = 2,
                                                  phi_death = 1,
                                                  k_death = 2,
                                                  exp_noise = 1e+07,
                                                  likelihood = function(model_deaths, data_deaths){
                                                    squire:::ll_nbinom(data_deaths, model_deaths, obs_params$phi_death,
                                                                       obs_params$k_death,
                                                                       obs_params$exp_noise)
                                                  }
                                                ),
                                                return = "ll") {

  if (!(return %in% c("full", "ll", "sample", "single"))) {
    stop("return argument must be full, ll, sample", "single")
  }
  if (as.Date(data$week_start[1], "%Y-%m-%d") <
      as.Date(model_start_date, "%Y-%m-%d")) {
    stop("Model start date is later than data start date")
  }

  #set up tolerances
  #if tolerances are null use 1e-6
  if(is.null(obs_params$atol)) {
    atol <- 1e-6
  } else {
    atol <- obs_params$atol
  }
  if(is.null(obs_params$rtol)) {
    rtol <- 1e-6
  } else {
    rtol <- obs_params$rtol
  }
  #if full we use a low tolerance
  if(return == "full") {
    atol <- 1e-8
    rtol <- 1e-8
  }

  #set up to use our weekly data instead of per day
  model_params$tt_beta <- round(model_params$tt_beta * model_params$dt)
  model_params$tt_contact_matrix <- round(model_params$tt_contact_matrix *
                                          model_params$dt)
  model_params$tt_hosp_beds <- round(model_params$tt_hosp_beds *
                                     model_params$dt)
  model_params$tt_ICU_beds <- round(model_params$tt_ICU_beds *
                                    model_params$dt)

  #convert weeks into days relevant to our start_date
  data$date <- data$week_start
  data <- squire:::particle_filter_data(data = data, start_date = model_start_date,
                                        steps_per_day = round(1/model_params$dt))
  data$week_start <- data$day_start
  data$week_end <- data$day_end

  #set the last day to the same distance as the previous one
  data$week_end[nrow(data)] <- data$week_start[nrow(data)] +
    data$week_end[nrow(data)-1]  - data$week_start[nrow(data)-1]

  # run model with wrapper so that if it fails we try again with lower tolerances
  model_func <- squire_model$odin_model(user = model_params,
                                        unused_user_action = "ignore")
  out <- tryCatch(
    model_func$run(t = seq(0, utils::tail(data$week_end, 1), 1), atol = atol, rtol = rtol),
    error = function(x){"FAIL"}
  )
  if(!identical(out, "FAIL")){
    #if it worked
    index <- squire:::odin_index(model_func)

    #calculate the deaths for each week
    cumDs <- rowSums(out[, index$D])
    Ds <- cumDs[data$week_end[-1]] - cumDs[data$week_start[-1]]
    Ds[Ds < 0] <- 0
    deaths <- data$deaths[-1]
    #write into cpp?
    ll <- obs_params$likelihood(Ds, deaths, obs_params)
    # and wrap up as normal
    date <- data$date[[1]] + seq_len(nrow(out)) - 1L
    rownames(out) <- as.character(date)
    attr(out, "date") <- date
    # allow full return for simulations
    if (return == "ll") {
      ret <- list(log_likelihood = sum(ll),
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
  return(ret)
}

#'
#'@noRd
evaluate_Rt_pmcmc_custom <- function(R0, pars, Rt_args){
  #calculate the values
  Rt <- as.numeric(c(R0, R0*2*stats::plogis(cumsum(-unlist(pars[grepl("Rt_rw", names(pars))])))))
  return(Rt)
}
