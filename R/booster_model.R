#' Create a nimue model class for fitting with squire tools using the booster
#' dose framework.
#'
#' @title nimue model creation.
#' @param use_dde Logical for using dde to solve. Default = TRUE
#' We will use this structure to ensure that model fitting is flexible in the
#' future as more models are added
#'
#' @details Wraps the squire pmcmc fitting infrastructure.
#' @export
nimue_booster_model <- function(use_dde = TRUE) {

  model_class <- "booster_model"

  compare_model <- function(model, pars_obs, data) {
    squire:::compare_output(model, pars_obs, data, type="nimue_model")
  }

  # wrap param func in order to remove unused arguments (dt)
  # and then add in all the default that are passed to params usually
  # from run so have to add here
  parameters_func <- function(country, population, dt,
                              contact_matrix_set, tt_contact_matrix,
                              hosp_bed_capacity, tt_hosp_beds,
                              ICU_bed_capacity, tt_ICU_beds,

                              # vaccine defaults that are just empty in parms so declare here
                              dur_R = vaccine_pars_booster$dur_R,
                              tt_dur_R = vaccine_pars_booster$tt_dur_R,
                              dur_V = vaccine_pars_booster$dur_V,
                              tt_dur_V = vaccine_pars_booster$tt_dur_V,
                              vaccine_efficacy_infection = vaccine_pars_booster$vaccine_efficacy_infection,
                              tt_vaccine_efficacy_infection = vaccine_pars_booster$tt_vaccine_efficacy_infection,
                              vaccine_efficacy_disease = vaccine_pars_booster$vaccine_efficacy_disease,
                              tt_vaccine_efficacy_disease = vaccine_pars_booster$tt_vaccine_efficacy_disease,
                              primary_doses = vaccine_pars_booster$primary_doses,
                              tt_primary_doses = vaccine_pars_booster$tt_primary_doses,
                              second_dose_delay = vaccine_pars_booster$second_dose_delay,
                              booster_doses = vaccine_pars_booster$booster_doses,
                              tt_booster_doses = vaccine_pars_booster$tt_booster_doses,
                              vaccine_coverage_mat = vaccine_pars_booster$vaccine_coverage_mat,
                              vaccine_booster_follow_up_coverage = vaccine_pars_booster$vaccine_booster_follow_up_coverage,
                              protection_delay_rate = vaccine_pars_booster$protection_delay_rate,
                              protection_delay_shape = vaccine_pars_booster$protection_delay_shape,
                              protection_delay_time = NULL,

                              #disease parameters
                              prob_hosp = probs_booster$prob_hosp,
                              prob_hosp_multiplier = probs_booster$prob_hosp_multiplier,
                              tt_prob_hosp_multiplier = probs_booster$tt_prob_hosp_multiplier,
                              prob_severe = probs_booster$prob_severe,
                              prob_severe_multiplier = probs_booster$prob_severe_multiplier,
                              tt_prob_severe_multiplier = probs_booster$tt_prob_severe_multiplier,
                              prob_non_severe_death_treatment = probs_booster$prob_non_severe_death_treatment,
                              prob_non_severe_death_no_treatment = probs_booster$prob_non_severe_death_no_treatment,
                              prob_severe_death_treatment = probs_booster$prob_severe_death_treatment,
                              prob_severe_death_no_treatment = probs_booster$prob_severe_death_no_treatment,
                              p_dist = probs_booster$p_dist,
                              rel_infectiousness = probs_booster$rel_infectiousness,
                              rel_infectiousness_vaccinated = probs_booster$rel_infectiousness_vaccinated,

                              # durations
                              dur_E  = durs_booster$dur_E,
                              tt_dur_E = durs_booster$tt_dur_E,
                              dur_IMild = durs_booster$dur_IMild,
                              tt_dur_IMild = durs_booster$tt_dur_IMild,
                              dur_ICase = durs_booster$dur_ICase,
                              tt_dur_ICase = durs_booster$tt_dur_ICase,

                              # hospital durations
                              dur_get_ox_survive = durs_booster$dur_get_ox_survive,
                              tt_dur_get_ox_survive = durs_booster$tt_dur_get_ox_survive,
                              dur_get_ox_die = durs_booster$dur_get_ox_die,
                              tt_dur_get_ox_die = durs_booster$tt_dur_get_ox_die,
                              dur_not_get_ox_survive = durs_booster$dur_not_get_ox_survive,
                              dur_not_get_ox_die = durs_booster$dur_not_get_ox_die,

                              dur_get_mv_survive = durs_booster$dur_get_mv_survive,
                              tt_dur_get_mv_survive = durs_booster$tt_dur_get_mv_survive,
                              dur_get_mv_die = durs_booster$dur_get_mv_die,
                              tt_dur_get_mv_die = durs_booster$tt_dur_get_mv_die,
                              dur_not_get_mv_survive = durs_booster$dur_not_get_mv_survive,
                              dur_not_get_mv_die = durs_booster$dur_not_get_mv_die,

                              dur_rec = durs_booster$dur_rec,

                              # seeding cases default
                              seeding_cases = 5,

                              ...) {

    pars <- parameters_booster(
      country = country,
      population = population,
      contact_matrix_set = contact_matrix_set,
      tt_contact_matrix = tt_contact_matrix,
      hosp_bed_capacity = hosp_bed_capacity,
      tt_hosp_beds = tt_hosp_beds,
      ICU_bed_capacity = ICU_bed_capacity,
      tt_ICU_beds = tt_ICU_beds,
      dur_E = dur_E,
      tt_dur_E = tt_dur_E,
      dur_IMild = dur_IMild,
      tt_dur_IMild = tt_dur_IMild,
      dur_ICase = dur_ICase,
      tt_dur_ICase = tt_dur_ICase,
      dur_get_ox_survive = dur_get_ox_survive,
      tt_dur_get_ox_survive = tt_dur_get_ox_survive,
      dur_get_ox_die = dur_get_ox_die,
      tt_dur_get_ox_die = tt_dur_get_ox_die,
      dur_not_get_ox_survive = dur_not_get_ox_survive,
      dur_not_get_ox_die = dur_not_get_ox_die,
      dur_get_mv_survive = dur_get_mv_survive,
      tt_dur_get_mv_survive = tt_dur_get_mv_survive,
      dur_get_mv_die = dur_get_mv_die,
      tt_dur_get_mv_die = tt_dur_get_mv_die,
      dur_not_get_mv_survive = dur_not_get_mv_survive,
      dur_not_get_mv_die = dur_not_get_mv_die,
      dur_rec = dur_rec,
      dur_R = dur_R,
      tt_dur_R = tt_dur_R,
      dur_V = dur_V,
      tt_dur_V = tt_dur_V,
      vaccine_efficacy_infection = vaccine_efficacy_infection,
      tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection,
      vaccine_efficacy_disease = vaccine_efficacy_disease,
      tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
      primary_doses = primary_doses,
      tt_primary_doses = tt_primary_doses,
      second_dose_delay = second_dose_delay,
      booster_doses = booster_doses,
      tt_booster_doses = tt_booster_doses,
      vaccine_coverage_mat = vaccine_coverage_mat,
      vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
      protection_delay_rate = protection_delay_rate,
      protection_delay_shape = protection_delay_shape,
      protection_delay_time = protection_delay_time,
      seeding_cases = seeding_cases,
      prob_hosp = prob_hosp,
      prob_hosp_multiplier = prob_hosp_multiplier,
      tt_prob_hosp_multiplier = tt_prob_hosp_multiplier,
      prob_severe = prob_severe,
      prob_severe_multiplier = prob_severe_multiplier,
      tt_prob_severe_multiplier = tt_prob_severe_multiplier,
      prob_non_severe_death_treatment = prob_non_severe_death_treatment,
      prob_non_severe_death_no_treatment = prob_non_severe_death_no_treatment,
      prob_severe_death_treatment = prob_severe_death_treatment,
      prob_severe_death_no_treatment = prob_severe_death_no_treatment,
      p_dist = p_dist,
      rel_infectiousness = rel_infectiousness,
      rel_infectiousness_vaccinated = rel_infectiousness_vaccinated,
      ...)

    # append extra pars for fitting
    pars$dt <- dt
    pars$prob_hosp_baseline <- pars$prob_hosp[1, ,1]
    pars$use_dde <- use_dde

    class(pars) <- c("vaccine_parameters", "squire_parameters")
    return(pars)

  }

  # wrap run func correctly
  run_func <- function(country, population, dt,
                       contact_matrix_set, tt_contact_matrix,
                       hosp_bed_capacity, tt_hosp_beds,
                       ICU_bed_capacity, tt_ICU_beds,
                       replicates = 1,
                       day_return = TRUE,
                       time_period = 365,
                       ...) {

    out <- squire.page:::run_booster(country = country,
               contact_matrix_set = contact_matrix_set,
               tt_contact_matrix = tt_contact_matrix,
               hosp_bed_capacity = hosp_bed_capacity,
               tt_hosp_beds = tt_hosp_beds,
               ICU_bed_capacity = ICU_bed_capacity,
               tt_ICU_beds = tt_ICU_beds,
               population = population,
               replicates = 1,
               time_period = time_period,
               use_dde = use_dde,
               ...)

    return(out)

  }

  odin_model <- function(user, unused_user_action) {
    nimue_booster$new(user = user, use_dde = use_dde, unused_user_action = "ignore")
  }

  model <- list(odin_model = odin_model,
                generate_beta_func = beta_est_infectiousness,
                parameter_func = parameters_func,
                run_func = run_func,
                compare_model = compare_model,
                use_dde = use_dde)
  class(model) <- c(model_class, "nimue_model", "deterministic", "squire_model")
  model

}

#' Return the default probabilities for modelling defined in \code{squire}
#' For more info see \href{squire parameters vignette}{https://mrc-ide.github.io/squire/articles/parameters.html}
#' @noRd
default_probs_booster <- function() {
  c(squire::default_probs(),
    list(rel_infectiousness = rep(1, 17),
         rel_infectiousness_vaccinated = matrix(
           c(0.5, 0.5, 1, 0.5, 0.5, 1), ncol = 17, nrow = 6,
           dimnames = list(c("pV_1", "fV_1", "fV_2", "bV_1", "bV_2", "bV_3"))
        ),
         prob_hosp_multiplier = 1,
         tt_prob_hosp_multiplier = 0,
         prob_severe_multiplier = 1,
         tt_prob_severe_multiplier = 0))
}

probs_booster <- default_probs_booster()

#' Get the default durations from nimue and add the time varying element
#' @noRd
default_durs_booster <- function() {
  c(
    nimue:::default_durations(),
    list(
      tt_dur_E = 0,
      tt_dur_IMild = 0,
      tt_dur_ICase = 0
    )
  )
}


durs_booster <- default_durs_booster()

#' Return the default vaccine parameters for modelling
#' @noRd
default_vaccine_pars_booster <- function() {
  #scale VE for breakthrough
  d <- c(0.75, 0.9, 0.049041065, 0.95, 0.14865027, 0.02109197)
  i <- c(0.55, 0.75, 0, 0.8, 0.1285048, 0)
  d <- (d - i)/(1 - i)
  list(dur_R = Inf,
       tt_dur_R = 0,
       dur_V = c(1/0.007466205, 1/0.00807429, 1/0.03331390),
       tt_dur_V = 0,
       vaccine_efficacy_infection = matrix(
         i, ncol = 17, nrow = 6,
         dimnames = list(c("pV_1", "fV_1", "fV_2", "bV_1", "bV_2", "bV_3"))
        ),
       tt_vaccine_efficacy_infection = 0,
       vaccine_efficacy_disease = matrix(
         d, ncol = 17, nrow = 6,
         dimnames = list(c("pV_1", "fV_1", "fV_2", "bV_1", "bV_2", "bV_3"))
       ),
       tt_vaccine_efficacy_disease = 0,
       primary_doses = 1000,
       tt_primary_doses = 0,
       second_dose_delay = 60,
       booster_doses = 100,
       tt_booster_doses = 0,
       vaccine_coverage_mat = matrix(0.8, ncol = 17, nrow = 1),
       vaccine_booster_follow_up_coverage = NULL,
       protection_delay_rate = 1/7,
       protection_delay_shape = 2)
}

vaccine_pars_booster <- default_vaccine_pars_booster()

#' Vaccine parameters
#'
#' @details All durations are in days.
#' @noRd
parameters_booster <- function(

  # Demography
  country = NULL,
  population = NULL,
  tt_contact_matrix = 0,
  contact_matrix_set = NULL,

  # Transmission
  R0 = 3,
  tt_R0 = 0,
  beta_set = NULL,

  # Initial state, duration, reps
  time_period = 365,
  seeding_cases,
  seeding_age_order = NULL,
  init = NULL,

  # Parameters
  # Probabilities
  prob_hosp,
  prob_hosp_multiplier,
  tt_prob_hosp_multiplier,
  prob_severe,
  prob_severe_multiplier,
  tt_prob_severe_multiplier,
  prob_non_severe_death_treatment,
  prob_non_severe_death_no_treatment,
  prob_severe_death_treatment,
  prob_severe_death_no_treatment,
  p_dist,

  rel_infectiousness,
  rel_infectiousness_vaccinated,

  # Durations
  dur_E,
  tt_dur_E,
  dur_IMild,
  tt_dur_IMild,
  dur_ICase,
  tt_dur_ICase,

  dur_get_ox_survive,
  tt_dur_get_ox_survive,
  dur_get_ox_die,
  tt_dur_get_ox_die,
  dur_not_get_ox_survive,
  dur_not_get_ox_die,

  dur_get_mv_survive,
  tt_dur_get_mv_survive,
  dur_get_mv_die,
  tt_dur_get_mv_die,
  dur_not_get_mv_survive,
  dur_not_get_mv_die,

  dur_rec,
  dur_R,
  tt_dur_R,

  # Vaccine
  dur_V,
  tt_dur_V,
  vaccine_efficacy_infection,
  tt_vaccine_efficacy_infection,
  vaccine_efficacy_disease,
  tt_vaccine_efficacy_disease,
  primary_doses,
  tt_primary_doses,
  booster_doses,
  tt_booster_doses,
  second_dose_delay,
  vaccine_coverage_mat,
  vaccine_booster_follow_up_coverage,
  protection_delay_rate,
  protection_delay_shape,
  protection_delay_time,

  # Health system capacity
  hosp_bed_capacity,
  ICU_bed_capacity,
  tt_hosp_beds,
  tt_ICU_beds


) {

  # Handle country population args
  cpm <- squire:::parse_country_population_mixing_matrix(country = country,
                                                         population = population,
                                                         contact_matrix_set = contact_matrix_set)
  country <- cpm$country
  population <- cpm$population
  contact_matrix_set <- cpm$contact_matrix_set

  # Standardise contact matrix set
  if(is.matrix(contact_matrix_set)){
    contact_matrix_set <- list(contact_matrix_set)
  }

  # populate contact matrix set if not provided
  if (length(contact_matrix_set) == 1) {
    baseline <- contact_matrix_set[[1]]
    contact_matrix_set <- vector("list", length(tt_contact_matrix))
    for(i in seq_along(tt_contact_matrix)) {
      contact_matrix_set[[i]] <- baseline
    }
  }


  # populate hospital and ICU bed capacity if not provided
  if (is.null(hosp_bed_capacity)) {
    if (!is.null(country)) {
      beds <- squire::get_healthcare_capacity(country)
      hosp_beds <- beds$hosp_beds
      hosp_bed_capacity <- rep(round(hosp_beds * sum(population)/1000), length(tt_hosp_beds))
    } else {
      hosp_bed_capacity <- round(5 * sum(population)/1000)
    }
  }
  if (is.null(ICU_bed_capacity)) {
    if (!is.null(country)) {
      beds <- squire::get_healthcare_capacity(country)
      ICU_beds <- beds$ICU_beds
      ICU_bed_capacity <- rep(round(ICU_beds * sum(population)/1000), length(tt_ICU_beds))
    } else {
      ICU_bed_capacity <- round(3 * hosp_bed_capacity/100)
    }
  }

  # Initial state and matrix formatting
  # ----------------------------------------------------------------------------

  # Initialise initial conditions
  mod_init <- nimue:::init(population, seeding_cases, seeding_age_order, init) %>%
  #add extra columns
  purrr::map(~cbind(.x, rep(0, nrow(.x))))
  ##remove extra columns
  #mod_init <- purrr::map(mod_init, ~.x[,6])

  # Convert contact matrices to input matrices
  matrices_set <- squire:::matrix_set_explicit(contact_matrix_set, population)

  # If a vector is put in for matrix targeting
  if(is.vector(vaccine_coverage_mat)){
    vaccine_coverage_mat <- matrix(vaccine_coverage_mat, ncol = 17)
  }

  # Input checks
  # ----------------------------------------------------------------------------
  mc <- squire:::matrix_check(population[-1], contact_matrix_set)
  stopifnot(length(R0) == length(tt_R0))
  stopifnot(length(contact_matrix_set) == length(tt_contact_matrix))
  stopifnot(length(hosp_bed_capacity) == length(tt_hosp_beds))
  stopifnot(length(ICU_bed_capacity) == length(tt_ICU_beds))
  stopifnot(length(primary_doses) == length(tt_primary_doses))
  stopifnot(length(booster_doses) == length(tt_booster_doses))
  stopifnot(length(prob_hosp_multiplier) == length(tt_prob_hosp_multiplier))
  stopifnot(length(prob_severe_multiplier) == length(tt_prob_severe_multiplier))
  stopifnot(length(dur_R) == length(tt_dur_R))
  stopifnot(length(dur_get_ox_survive) == length(tt_dur_get_ox_survive))
  stopifnot(length(dur_get_ox_die) == length(tt_dur_get_ox_die))
  stopifnot(length(dur_get_mv_survive) == length(tt_dur_get_mv_survive))
  stopifnot(length(dur_E) == length(tt_dur_E))
  stopifnot(length(dur_IMild) == length(tt_dur_IMild))
  stopifnot(length(dur_ICase) == length(tt_dur_ICase))
  stopifnot(ncol(vaccine_coverage_mat) == 17)

  nimue:::assert_pos(dur_E)
  nimue:::assert_pos(dur_IMild)
  nimue:::assert_pos(dur_ICase)
  nimue:::assert_pos(dur_get_ox_survive)
  nimue:::assert_pos(dur_get_ox_die)
  nimue:::assert_pos(dur_not_get_ox_survive)
  nimue:::assert_pos(dur_not_get_ox_die)
  nimue:::assert_pos(dur_get_mv_survive)
  nimue:::assert_pos(dur_get_mv_die)
  nimue:::assert_pos(dur_not_get_mv_survive)
  nimue:::assert_pos(dur_not_get_mv_die)
  nimue:::assert_pos(dur_R)
  nimue:::assert_pos(time_period)
  nimue:::assert_pos(hosp_bed_capacity)
  nimue:::assert_pos(ICU_bed_capacity)
  nimue:::assert_pos(primary_doses)
  nimue:::assert_pos(booster_doses)
  nimue:::assert_pos(second_dose_delay)
  nimue:::assert_pos(prob_hosp_multiplier)
  nimue:::assert_pos(prob_severe_multiplier)

  #cannot have 1 in coverage matrix
  if(any(vaccine_coverage_mat == 1)){
    stop("vaccine_coverage_mat cannot have any element == 1")
  }
  #check inclusive, each row is greater than or equal to the previous
  if(any(vaccine_coverage_mat[seq_len(nrow(vaccine_coverage_mat) - 1),] > vaccine_coverage_mat[seq_len(nrow(vaccine_coverage_mat))[-1],])){
    stop("each row of vaccine_coverage_mat must include the previous row")
  }

  nimue:::assert_length(prob_hosp, length(population))
  nimue:::assert_length(prob_severe, length(population))
  nimue:::assert_length(prob_non_severe_death_treatment, length(population))
  nimue:::assert_length(prob_non_severe_death_no_treatment, length(population))
  nimue:::assert_length(prob_severe_death_treatment, length(population))
  nimue:::assert_length(prob_severe_death_no_treatment, length(population))
  nimue:::assert_length(rel_infectiousness, length(population))
  nimue:::assert_length(p_dist, length(population))

  nimue:::assert_numeric(prob_hosp, length(population))
  nimue:::assert_numeric(prob_severe, length(population))
  nimue:::assert_numeric(prob_non_severe_death_treatment, length(population))
  nimue:::assert_numeric(prob_non_severe_death_no_treatment, length(population))
  nimue:::assert_numeric(prob_severe_death_treatment, length(population))
  nimue:::assert_numeric(prob_severe_death_no_treatment, length(population))
  nimue:::assert_numeric(rel_infectiousness, length(population))
  nimue:::assert_numeric(p_dist, length(population))

  nimue:::assert_leq(prob_hosp, 1)
  nimue:::assert_leq(prob_severe, 1)
  nimue:::assert_leq(prob_non_severe_death_treatment, 1)
  nimue:::assert_leq(prob_non_severe_death_no_treatment, 1)
  nimue:::assert_leq(prob_severe_death_treatment, 1)
  nimue:::assert_leq(prob_severe_death_no_treatment, 1)
  nimue:::assert_leq(rel_infectiousness, 1)
  nimue:::assert_leq(p_dist, 1)

  nimue:::assert_greq(prob_hosp, 0)
  nimue:::assert_greq(prob_severe, 0)
  nimue:::assert_greq(prob_non_severe_death_treatment, 0)
  nimue:::assert_greq(prob_non_severe_death_no_treatment, 0)
  nimue:::assert_greq(prob_severe_death_treatment, 0)
  nimue:::assert_greq(prob_severe_death_no_treatment, 0)
  nimue:::assert_greq(rel_infectiousness, 0)
  nimue:::assert_greq(p_dist, 0)

  if(is.null(vaccine_booster_follow_up_coverage)){
    vaccine_booster_follow_up_coverage <- rep(1, 17)
  } else {
    if(any(!vaccine_booster_follow_up_coverage %in% c(0, 1))){
      stop("vaccine_booster_follow_up_coverage must be NULL for a vector of 0s and 1s with length = N_age")
    }
  }


  # Convert and Generate Parameters As Required
  # ----------------------------------------------------------------------------

  # durations
  gamma_E = 2 * 1/dur_E
  gamma_IMild = 1/dur_IMild
  gamma_ICase = 2 * 1/dur_ICase
  gamma_get_ox_survive = 2 * 1/dur_get_ox_survive
  gamma_get_ox_die = 2 * 1/dur_get_ox_die
  gamma_not_get_ox_survive = 2 * 1/dur_not_get_ox_survive
  gamma_not_get_ox_die = 2 * 1/dur_not_get_ox_die
  gamma_get_mv_survive = 2 * 1/dur_get_mv_survive
  gamma_get_mv_die = 2 * 1/dur_get_mv_die
  gamma_not_get_mv_survive = 2 * 1/dur_not_get_mv_survive
  gamma_not_get_mv_die = 2 * 1/dur_not_get_mv_die
  gamma_rec = 2 * 1/dur_rec
  gamma_R <- 2 * 1/dur_R

  if (is.null(beta_set)) {
    baseline_matrix <- squire:::process_contact_matrix_scaled_age(contact_matrix_set[[1]], population)
    #check for time changing parameters
    if(length(c(tt_dur_ICase, tt_dur_IMild, tt_prob_hosp_multiplier)) > 3){
      tt_R0_old <- tt_R0
      tt_R0 <- unique(sort(c(tt_dur_ICase, tt_dur_IMild, tt_prob_hosp_multiplier, tt_R0_old)))
      R0 <- block_interpolate(tt_R0, R0, tt_R0_old)
    }
    beta_set <- beta_est_booster(
      R0 = R0, tt_R0 = tt_R0, prob_hosp_multiplier = prob_hosp_multiplier,
      tt_prob_hosp_multiplier = tt_prob_hosp_multiplier,
      prob_hosp_baseline = prob_hosp, dur_ICase = dur_ICase,
      tt_dur_ICase = tt_dur_ICase, dur_IMild = dur_IMild,
      tt_dur_IMild = tt_dur_IMild, rel_infectiousness = rel_infectiousness,
      mixing_matrix = baseline_matrix
    )
  }

  # normalise to sum to 1
  p_dist <- matrix(rep(p_dist, 5), nrow = 17, ncol = 5)
  p_dist <- p_dist/mean(p_dist)

  # Format vaccine-specific parameters
  if(typeof(dur_V) == "list"){
    gamma_vaccine <- purrr::map(dur_V, ~c(0, 0, 1/.x[1], 0, 1/.x[-1], 0)) %>%
      unlist() %>%
      matrix(ncol = 7, nrow = length(tt_dur_V), byrow = TRUE)
    tt_dur_vaccine <- tt_dur_V
  } else {
    gamma_vaccine <- matrix(c(0, 0, 1/dur_V[1], 0, 1/dur_V[-1], 0), nrow = 1)
    tt_dur_vaccine <- 0
  }

  rel_infectiousness_vaccinated <- format_rel_inf_vacc_for_odin_booster(rel_infectiousness_vaccinated)

  # Vaccine efficacies are now time changing (if specified),
  # so we need to convert these to be interpolated by odin
  # These functions also check that efficacies are correct length
  # both in terms of age groups and in terms of required timepoints

  # First the vaccine efficacy infection
  vaccine_efficacy_infection_odin_array <- format_ve_i_for_odin_booster(
    vaccine_efficacy_infection = vaccine_efficacy_infection,
    tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection
  )

  # Second the vaccine efficacy disease affecting prob_hosp
  prob_hosp_odin_array <- format_ve_d_for_odin_booster(
    vaccine_efficacy_disease = vaccine_efficacy_disease,
    tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
    prob_hosp = prob_hosp
  )

  ##Delay dosing
  if(is.null(protection_delay_time)){
    protection_delay_time <- time_period
  }
  delayed <- apply_dose_delay_booster(seq_len(protection_delay_time) - 1, primary_doses, tt_primary_doses,
                                      booster_doses, tt_booster_doses,
                                      second_dose_delay, protection_delay_rate,
                                      protection_delay_shape)
  primary_doses <- delayed$primary_doses
  tt_primary_doses <- delayed$tt_primary_doses
  booster_doses <- delayed$booster_doses
  tt_booster_doses <- delayed$tt_booster_doses
  second_dose_delay <-delayed$second_dose_delay

  # Collate Parameters Into List
  pars <- c(mod_init,
            list(N_age = length(population),
                 gamma_E = gamma_E,
                 tt_dur_E = tt_dur_E,
                 gamma_IMild = gamma_IMild,
                 tt_dur_IMild = tt_dur_IMild,
                 gamma_ICase = gamma_ICase,
                 tt_dur_ICase = tt_dur_ICase,
                 gamma_get_ox_survive = gamma_get_ox_survive,
                 tt_dur_get_ox_survive = tt_dur_get_ox_survive,
                 gamma_get_ox_die = gamma_get_ox_die,
                 tt_dur_get_ox_die = tt_dur_get_ox_die,
                 gamma_not_get_ox_survive = gamma_not_get_ox_survive,
                 gamma_not_get_ox_die = gamma_not_get_ox_die,
                 gamma_get_mv_survive = gamma_get_mv_survive,
                 tt_dur_get_mv_survive = tt_dur_get_mv_survive,
                 gamma_get_mv_die = gamma_get_mv_die,
                 tt_dur_get_mv_die = tt_dur_get_mv_die,
                 gamma_not_get_mv_survive = gamma_not_get_mv_survive,
                 gamma_not_get_mv_die = gamma_not_get_mv_die,
                 gamma_rec = gamma_rec,
                 gamma_R = gamma_R,
                 tt_dur_R = tt_dur_R,
                 prob_hosp = prob_hosp_odin_array,
                 prob_hosp_multiplier = prob_hosp_multiplier,
                 tt_prob_hosp_multiplier = tt_prob_hosp_multiplier,
                 prob_severe = prob_severe,
                 prob_severe_multiplier = prob_severe_multiplier,
                 tt_prob_severe_multiplier = tt_prob_severe_multiplier,
                 prob_non_severe_death_treatment = prob_non_severe_death_treatment,
                 prob_non_severe_death_no_treatment = prob_non_severe_death_no_treatment,
                 prob_severe_death_treatment = prob_severe_death_treatment,
                 prob_severe_death_no_treatment = prob_severe_death_no_treatment,
                 rel_infectiousness = rel_infectiousness,
                 rel_infectiousness_vaccinated = rel_infectiousness_vaccinated,
                 p_dist = p_dist,
                 hosp_beds = hosp_bed_capacity,
                 ICU_beds = ICU_bed_capacity,
                 tt_hosp_beds = tt_hosp_beds,
                 tt_ICU_beds = tt_ICU_beds,
                 tt_matrix = tt_contact_matrix,
                 mix_mat_set = matrices_set,
                 tt_beta = tt_R0,
                 beta_set = beta_set,
                 population = population,
                 contact_matrix_set = contact_matrix_set,
                 primary_doses = primary_doses,
                 tt_primary_doses = tt_primary_doses,
                 booster_doses = booster_doses,
                 tt_booster_doses = tt_booster_doses,
                 second_dose_delay = second_dose_delay,
                 vaccine_efficacy_infection = vaccine_efficacy_infection_odin_array,
                 tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection,
                 tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
                 vaccine_coverage_mat = vaccine_coverage_mat,
                 vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
                 N_vaccine = 7,
                 N_prioritisation_steps = nrow(vaccine_coverage_mat),
                 gamma_vaccine = gamma_vaccine,
                 tt_dur_vaccine = tt_dur_vaccine))

  class(pars) <- c("booster_vaccine_parameters", "vaccine_parameters", "nimue_parameters")

  return(pars)
}

#' @noRd
apply_dose_delay_booster <- function(t, primary_doses, tt_primary_doses,
                                     booster_doses, tt_booster_doses,
                                     second_dose_delay, protection_delay_rate,
                                     protection_delay_shape){
  if(!is.null(protection_delay_rate) & !is.null(protection_delay_shape)){
    #interpolate
    if(any(primary_doses > 0)){
      primary_doses_int <- block_interpolate(t, primary_doses, tt_primary_doses)
      primary_doses <- diff(c(0,
                              purrr::map_dbl(seq_along(primary_doses_int), function(t){
                                sum(stats::pgamma(seq(t - 1, 0), shape = protection_delay_shape, rate = protection_delay_rate) *
                                      primary_doses_int[seq_len(t)])
                              })
      ))
      tt_primary_doses <- t
    }
    if(any(booster_doses > 0)){
      booster_doses_int <- block_interpolate(t, booster_doses, tt_booster_doses)
      booster_doses <- diff(c(0,
                              purrr::map_dbl(seq_along(booster_doses_int), function(t){
                                sum(stats::pgamma(seq(t - 1, 0), shape = protection_delay_shape, rate = protection_delay_rate) *
                                      booster_doses_int[seq_len(t)])
                              })
      ))
      tt_booster_doses <- t
    }
    #add mean to second dose delay
    second_dose_delay <- second_dose_delay + protection_delay_shape/protection_delay_rate
  }
  return(
    list(
      primary_doses = primary_doses, tt_primary_doses = tt_primary_doses,
      booster_doses = booster_doses, tt_booster_doses = tt_booster_doses,
      second_dose_delay = second_dose_delay
    )
  )
}

#' @noRd
format_rel_inf_vacc_for_odin_booster <- function(rel_inf_vacc) {

  #if only have one assume this holds for all ages/compartments
  if(length(rel_inf_vacc) == 1){
    rel_inf_vacc <- matrix(rel_inf_vacc, ncol = 17, nrow = 6)
  } else if (is.numeric(rel_inf_vacc) & length(rel_inf_vacc) == 6){
    #expand across age groups
    rel_inf_vacc <- matrix(rel_inf_vacc, ncol = 17, nrow = 6)
  } else if (is.numeric(rel_inf_vacc) & length(rel_inf_vacc) == 17){
    #expand across vaccine comparments
    rel_inf_vacc <- matrix(rel_inf_vacc, ncol = 17, nrow = 6, byrow = TRUE)
  } else if (!(nrow(rel_inf_vacc) == 6 & ncol(rel_inf_vacc) == 17)){
    stop("rel_infectiousness_vaccinated must be a single value, a vector of length 6 or 17, or a matrix with 6 rows and 17 columns")
  }


  #add 1 for unvaccinated, then rotate to match model requirements
  return(
    t(rbind(
      rep(1, 17),
      rel_inf_vacc
    ))
  )

}

#' @noRd
format_ve_i_for_odin_booster <- function(vaccine_efficacy_infection,
                                 tt_vaccine_efficacy_infection) {

  # If just provided as a vector then we put into a list ready for formatting
  if(!is.list(vaccine_efficacy_infection)){
    vaccine_efficacy_infection <- list(vaccine_efficacy_infection)
  }

  # check that the correct length agreement between tt_vaccine_efficacy_infection
  nimue:::assert_length(vaccine_efficacy_infection, length(vaccine_efficacy_infection))

  # now check that each vaccine efficacy is correct number of columns (1 or 17)
  vaccine_efficacy_infection <- lapply(vaccine_efficacy_infection, function(ve_i) {

    #if numeric vector
    if(any(class(ve_i) == "numeric")) {
      if(length(ve_i) == 6) {
        #make into matrix
        ve_i <- matrix(ve_i, ncol = 17, nrow = 6)
      } else {
        stop("If element of vaccine_efficacy_infection is a vector, it must have 6 values corresponding to a first dose, second dose, a waned compartment, booster dose, and the two waning levels")
      }
    }

    if(ncol(ve_i) != 17 | nrow(ve_i) != 6){
      stop("Parameter vaccine_efficacy_infection must be vector of length 6 or a matrix with ncol = 17, nrow = 6")
    }

    return(ve_i)

  })

  # and now format so each list is the vaccine_efficacy_infection at each time
  # point for the 5 vaccine classes
  ve_i_list <- lapply(seq_along(tt_vaccine_efficacy_infection), function(ve_i_index) {
    ve_i <-
      rbind(
        vaccine_efficacy_infection[[ve_i_index]]
      )
    #add 0 for unvaccinated
    ve_i <- rbind(
      rep(0, ncol(ve_i)),
      ve_i
    )
    ve_i = 1 - ve_i
  })

  # and use this list to create an array that is in right format for odin
  vaccine_efficacy_infection_odin_array <- aperm(
    array(unlist(ve_i_list), dim = c(dim(ve_i_list[[1]]), length(ve_i_list))),
    c(3, 2, 1)
  )

  return(vaccine_efficacy_infection_odin_array)

}


#' @noRd
format_ve_d_for_odin_booster <- function(vaccine_efficacy_disease,
                                 tt_vaccine_efficacy_disease,
                                 prob_hosp) {


  # If just provided as a vector then we put into a list ready for formatting
  if(!is.list(vaccine_efficacy_disease)){
    vaccine_efficacy_disease <- list(vaccine_efficacy_disease)
  }

  # check that the correct length agreement between tt_vaccine_efficacy_disease
  nimue:::assert_length(vaccine_efficacy_disease, length(vaccine_efficacy_disease))

  # now check that each vaccine efficacy is correct length (1 or 17)
  vaccine_efficacy_disease <- lapply(vaccine_efficacy_disease, function(ve_d) {

    #if numeric vector
    if(any(class(ve_d) == "numeric")) {
      if(length(ve_d) == 6) {
        #make into matrix
        ve_d <- matrix(ve_d, ncol = 17, nrow = 6)
      } else {
        stop("If element of vaccine_efficacy_disease is a vector, it must have 6 values corresponding to a first dose, a waned compartment, second dose, and the two waning levels")
      }
    }

    if(ncol(ve_d) != 17 | nrow(ve_d) != 6){
      stop("Parameter vaccine_efficacy_disease must be vector of length 6 or a matrix with ncol = 17, nrow = 6")
    }

    return(ve_d)

  })
  # and now format so each list is the prob_hosp at each time
  # point for the 5 vaccine classes
  prob_hosp_list <- lapply(seq_along(tt_vaccine_efficacy_disease), function(ve_d_index) {

    ve_d <-
      vaccine_efficacy_disease[[ve_d_index]]


    prob_hosp_vaccine <-
      sweep(
        (1 - ve_d),
        MARGIN = 2,
        STATS = prob_hosp,
        FUN = "*"
      )

    # add the baseline for unvaccinated
    prob_hosp <- rbind(
      prob_hosp,
      prob_hosp_vaccine
    )

    return(prob_hosp)

  })

  # and use this list to create an array that is in right format for odin
  prob_hosp_odin_array <- aperm(
    array(unlist(prob_hosp_list), dim = c(dim(prob_hosp_list[[1]]), length(prob_hosp_list))),
    c(3, 2, 1)
  )

  return(prob_hosp_odin_array)

}

#' Run the vaccine model
#'
#' @param population Population vector (for each age group). Default = NULL,
#'   which will cause population to be sourced from \code{country}
#' @param country Character for country beign simulated. WIll be used to
#'   generate \code{population} and \code{contact_matrix_set} if
#'   unprovided. Either \code{country} or \code{population} and
#'   \code{contact_matrix_set} must be provided.
#' @param contact_matrix_set Contact matrices used in simulation. Default =
#'   NULL, which will generate this based on the \code{country}.
#' @param tt_contact_matrix Time change points for matrix change. Default = 0
#' @param R0 Basic Reproduction Number. Default = 3
#' @param tt_R0 Change time points for R0. Default = 0
#' @param beta_set Alternative parameterisation via beta rather than R0.
#'   Default = NULL, which causes beta to be estimated from R0
#' @param time_period Length of simulation. Default = 365
#' @param replicates  Number of replicates. Default = 10
#' @param seeding_cases Initial number of cases seeding the epidemic
#' @param seed Random seed used for simulations. Deafult = runif(1, 0, 10000)
#' @param prob_hosp probability of hospitalisation by age.
#'   Default = c(0.000744192, 0.000634166,0.001171109, 0.002394593, 0.005346437,
#'   0.010289885, 0.016234604, 0.023349169, 0.028944623, 0.038607042,
#'   0.057734879, 0.072422135, 0.101602458, 0.116979814, 0.146099064,
#'   0.176634654 ,0.180000000)
#' @param prob_hosp_multiplier Time varying multiplier to probability of developing
#' severe symptoms. Default = 1, which is no change to provided prob_hosp.
#' @param tt_prob_hosp_multiplier Timing of changes to multiplier of probability of
#' developing severe symptoms. Default = 0
#' @param prob_severe_multiplier Time varying multiplier to probability of
#'   hospitalisation. Default = 1, which is no change to provided prob_hosp.
#' @param tt_prob_severe_multiplier Timing of changes to multiplier of probability
#'   of hospitalisation. Default = 0
#' @param prob_severe Probability of developing severe symptoms by age.
#'   Default = c(0.05022296,	0.05022296,	0.05022296,	0.05022296,	0.05022296,
#'   0.05022296,	0.05022296,	0.053214942, 0.05974426,	0.074602879,
#'   0.103612417, 0.149427991, 0.223777304,	0.306985918,
#'   0.385779555, 0.461217861, 0.709444444)
#' @param prob_non_severe_death_treatment Probability of death from non severe
#'   treated infection.
#'   Default = c(0.0125702,	0.0125702,	0.0125702,	0.0125702,
#'   0.0125702,	0.0125702,	0.0125702,	0.013361147,
#'   0.015104687,	0.019164124,	0.027477519,	0.041762108,
#'   0.068531658,	0.105302319,	0.149305732,	0.20349534,	0.5804312)
#' @param prob_severe_death_treatment Probability of death from severe infection
#'   that is treated. Default = rep(0.5, 17)
#' @param prob_non_severe_death_no_treatment Probability of death in non severe
#'   hospital inections that aren't treated
#' @param prob_severe_death_no_treatment Probability of death from severe infection
#'   that is not treated. Default = rep(0.95, 17)
#' @param p_dist Preferentiality of age group receiving treatment relative to
#'   other age groups when demand exceeds healthcare capacity.
#' @param rel_infectiousness Relative infectiousness per age category relative
#'   to maximum infectiousness category. Default = rep(1, 17)
#' @param rel_infectiousness_vaccinated  Relative infectiousness per age
#'   category  of vaccinated individuals relative to unvaccinated individuals.
#'   Default = rep(1, 17), which is no impact of vaccination on onwards
#'   transmissions
#' @param dur_E Mean duration of incubation period (days). Default = 4.6
#' @param tt_dur_E Times at which dur_E changes, default = 0.
#' @param dur_IMild Mean duration of mild infection (days). Default = 2.1
#' @param tt_dur_IMild Times at which dur_IMild changes, default = 0.
#' @param dur_ICase Mean duration from symptom onset to hospital admission (days).
#'   Default = 4.5
#' @param tt_dur_ICase Times at which dur_ICase changes, default = 0.
#' @param dur_get_ox_survive Mean duration of oxygen given survive. Default = 5. Can be
#'   time varying, with timing of changes given by tt_dur_get_ox_survive.
#' @param tt_dur_get_ox_survive Timing of changes in duration of  oxygen given survive.
#' @param dur_get_ox_die Mean duration of oxygen given death. Default = 5. Can be
#'   time varying, with timing of changes given by tt_dur_get_ox_die.
#' @param tt_dur_get_ox_die Timing of changes in duration of  oxygen given death.
#' @param dur_not_get_ox_survive Mean duration without oxygen given survive.
#'   Default = 5
#' @param dur_not_get_ox_die Mean duration without  oxygen given death.
#'  Default = 5
#' @param dur_get_mv_survive Mean duration of ventilation given survive.
#'   Default = 7.3. Can be time varying, with timing of changes given by tt_dur_get_mv_survive.
#' @param tt_dur_get_mv_survive Timing of changes in duration of ventilation given survive.
#' @param dur_get_mv_die Mean duration of ventilation given death. Default = 6. Can be
#'   time varying, with timing of changes given by tt_dur_get_mv_die.
#' @param tt_dur_get_mv_die Timing of changes in duration of ventilation given death.
#' @param dur_not_get_mv_survive Mean duration without ventilation given
#'   survive. Default = 7.3
#' @param dur_not_get_mv_die Mean duration without ventilation given
#'   death. Default = 1
#' @param dur_rec Duration of recovery after coming off ventilation. Default = 2
#' @param hosp_bed_capacity General bed capacity. Can be single number of vector if capacity time-varies.
#' @param ICU_bed_capacity ICU bed capacity. Can be single number of vector if capacity time-varies.
#' @param tt_hosp_beds Times at which hospital bed capacity changes (Default = 0 = doesn't change)
#' @param tt_ICU_beds Times at which ICU bed capacity changes (Default = 0 = doesn't change)
#' @param seeding_cases Initial number of cases seeding the epidemic
#' @param seeding_age_order Vector specifying the order in which seeds are allocated to ages.
#'   If NULL, seeds are distributed randomly within working ages. If specified, must be a vector
#'   of length 17 specifying the order seeds are allocated, e.g. 1:17 will allocate first seed
#'   to the youngest age group, then the second youngest and so on. Default = NULL
#' @param init Initial conditions for simulation provided. Allows overriding
#'   if initial conditions start with an already infected population etc.
#'   Default = NULL.
#' @param dur_R Mean duration of naturally acquired immunity (days). Can be
#'   time varying, with timing of changes given by tt_dur_R.
#' @param tt_dur_R Timing of changes in duration of natural immunity.
#' @param dur_V Mean duration of vaccine-derived immunity (days) for partial protection and full protection. Should be a
#'   numeric vector of length 3, corresponding to the duration of time in each waned compartmenet after recieving a first dose and then for the two second dose compartments.
#'   Alternatively can be a list of values if this changes over time.
#' @param tt_dur_V List of change times for dur_V.
#' @param vaccine_efficacy_infection Efficacy of vaccine against infection.
#'   This parameter must either be a length 6 numeric (a single efficacy for
#'   each vaccine state (first dose, second dose, waned second dose, booster dose, and two waned second dose compartments))
#'   or numeric vector with 17 columns and 6 rows
#'   (an efficacy for each age group and vaccine state).
#'   An efficacy of 1 will reduce FOI by 100 percent, an efficacy of 0.2 will
#'   reduce FOI by 20 percent etc.
#'   To specify changes in vaccine efficacy over time, vaccine efficacies must
#'   be provided as a list, with each list element being the efficacy at each
#'   time point specified by \code{tt_vaccine_efficacy_infection}. These
#'   efficacies must also be length 6 numeric or 6x17 numeric matrix.
#' @param tt_vaccine_efficacy_infection Timing of changes in vaccine efficacy
#'   against infection. Default = 0, which assumes fixed efficacy over time.
#'   Must be the same length as the length of \code{vaccine_efficacy_infection}
#'   when provided as a list. Time changing efficacies can occur in response to
#'   changing vaccines being  given and dosing strategy changes.
#' @param vaccine_efficacy_disease Efficacy of partial vaccination against severe
#'   (requiring hospitilisation) disease (by age). This parameter must either be
#'   length 6 numeric (a single efficacy for each vaccine state (first dose, second dose, waned second dose, booster dose, and two waned second dose compartments)) or numeric vector with 17 columns and 6 rows
#'   (an efficacy for each age group and vaccine state). An efficacy of 1 will
#'   reduce the probability of hospitalisation by 100 percent, an efficacy of
#'   0.2 will reduce the probability of hospitalisation by 20 percent etc.
#'   To specify changes in vaccine efficacy over time, vaccine efficacies must
#'   be provided as a list, with each list element being the efficacy at each
#'   time point specified by \code{tt_vaccine_efficacy_disease}. These
#'   efficacies must also be length 6 numeric or 5x17 numeric matrix.
#' @param tt_vaccine_efficacy_disease Timing of changes in vaccine efficacy
#'   against disease. Default = 0, which assumes fixed efficacy over time.
#'   Must be the same length as the length of \code{vaccine_efficacy_disease}
#'   when provided as a list. Time changing efficacies can occur in response to
#'   changing vaccines being  given and dosing strategy changes.
#' @param primary_doses The maximum number of individuals who can be vaccinated with their first dose per day.
#' @param tt_primary_doses Time change points for vaccine capacity (\code{first_doses}).
#' @param second_dose_delay Delay between first dose and second dose in the initial series, this model assume all who get first doses get a second dose, default = 60 (days).
#' @param booster_doses The maximum number of individuals who can be vaccinated with their booster dose per day.
#' @param tt_booster_doses Time change points for vaccine capacity (\code{booster_doses}).
#' @param vaccine_coverage_mat Vaccine coverage targets by age (columns) and priority (row)
#' @param vaccine_booster_follow_up_coverage Age group eligibility for follow-up boosters (i.e. 2nd, 3rd, ... booster doses),
#' default = NULL means all are eligible. Format: 0 indicates not eligible, 1 indicates eligible.
#' @param protection_delay_rate Rate for the delay in development of vaccine protection, applied via gamma/erlang distribution,
#' default = 1/14. If NULL no delay is applied.
#' @param protection_delay_shape Shape for the delay in development of vaccine protection, applied via gamma/erlang distribution,
#' default = 1/14. If NULL no delay is applied.
#' @param use_dde Use the dde solver (default is \code{TRUE})
#' @param ... Additional arguments for solver
#'
#' @return Simulation output
#' @noRd
run_booster <- function(

  # demography
  country = NULL,
  population = NULL,
  tt_contact_matrix = 0,
  contact_matrix_set = NULL,

  # transmission
  R0 = 3,
  tt_R0 = 0,
  beta_set = NULL,

  # initial state, duration, reps
  time_period = 365,
  replicates = 10,
  seed = stats::runif(1, 0, 100000000),

  # parameters
  # probabilities
  prob_hosp = probs_booster$prob_hosp,
  prob_hosp_multiplier = probs_booster$prob_hosp_multiplier,
  tt_prob_hosp_multiplier = probs_booster$tt_prob_hosp_multiplier,
  prob_severe = probs_booster$prob_severe,
  prob_severe_multiplier = probs_booster$prob_severe_multiplier,
  tt_prob_severe_multiplier = probs_booster$tt_prob_severe_multiplier,
  prob_non_severe_death_treatment = probs_booster$prob_non_severe_death_treatment,
  prob_non_severe_death_no_treatment = probs_booster$prob_non_severe_death_no_treatment,
  prob_severe_death_treatment = probs_booster$prob_severe_death_treatment,
  prob_severe_death_no_treatment = probs_booster$prob_severe_death_no_treatment,
  p_dist = probs_booster$p_dist,

  # onward infectiousness
  rel_infectiousness = probs_booster$rel_infectiousness,
  rel_infectiousness_vaccinated = probs_booster$rel_infectiousness_vaccinated,

  # durations
  dur_E  = durs_booster$dur_E,
  tt_dur_E  = durs_booster$tt_dur_E,
  dur_IMild = durs_booster$dur_IMild,
  tt_dur_IMild = durs_booster$tt_dur_IMild,
  dur_ICase = durs_booster$dur_ICase,
  tt_dur_ICase = durs_booster$tt_dur_ICase,

  # hospital durations
  dur_get_ox_survive = durs_booster$dur_get_ox_survive,
  tt_dur_get_ox_survive = durs_booster$tt_dur_get_ox_survive,
  dur_get_ox_die = durs_booster$dur_get_ox_die,
  tt_dur_get_ox_die = durs_booster$tt_dur_get_ox_die,
  dur_not_get_ox_survive = durs_booster$dur_not_get_ox_survive,
  dur_not_get_ox_die = durs_booster$dur_not_get_ox_die,

  dur_get_mv_survive = durs_booster$dur_get_mv_survive,
  tt_dur_get_mv_survive = durs_booster$tt_dur_get_mv_survive,
  dur_get_mv_die = durs_booster$dur_get_mv_die,
  tt_dur_get_mv_die = durs_booster$tt_dur_get_mv_die,
  dur_not_get_mv_survive = durs_booster$dur_not_get_mv_survive,
  dur_not_get_mv_die = durs_booster$dur_not_get_mv_die,

  dur_rec = durs_booster$dur_rec,

  # vaccine
  dur_R = vaccine_pars_booster$dur_R,
  tt_dur_R = vaccine_pars_booster$tt_dur_R,
  dur_V = vaccine_pars_booster$dur_V,
  tt_dur_V = vaccine_pars_booster$tt_dur_V,
  vaccine_efficacy_infection = vaccine_pars_booster$vaccine_efficacy_infection,
  tt_vaccine_efficacy_infection = vaccine_pars_booster$tt_vaccine_efficacy_infection,
  vaccine_efficacy_disease = vaccine_pars_booster$vaccine_efficacy_disease,
  tt_vaccine_efficacy_disease = vaccine_pars_booster$tt_vaccine_efficacy_disease,
  primary_doses = vaccine_pars_booster$primary_doses,
  tt_primary_doses = vaccine_pars_booster$tt_primary_doses,
  booster_doses = vaccine_pars_booster$booster_doses,
  tt_booster_doses = vaccine_pars_booster$tt_booster_doses,
  second_dose_delay = vaccine_pars_booster$second_dose_delay,
  vaccine_coverage_mat = vaccine_pars_booster$vaccine_coverage_mat,
  vaccine_booster_follow_up_coverage = vaccine_pars_booster$vaccine_booster_follow_up_coverage,
  protection_delay_rate = vaccine_pars_booster$protection_delay_rate,
  protection_delay_shape = vaccine_pars_booster$protection_delay_shape,

  # health system capacity
  hosp_bed_capacity = NULL,
  ICU_bed_capacity = NULL,
  tt_hosp_beds = 0,
  tt_ICU_beds = 0,

  seeding_cases = 20,
  seeding_age_order = NULL,
  init = NULL,
  use_dde = TRUE,
  ...
) {

  # Grab function arguments
  args <- as.list(environment())
  set.seed(seed)

  # create parameter list
  pars <- parameters_booster(country = country,
                     population = population,
                     tt_contact_matrix = tt_contact_matrix,
                     contact_matrix_set = contact_matrix_set,
                     R0 = R0,
                     tt_R0 = tt_R0 ,
                     beta_set = beta_set,
                     time_period = time_period,
                     seeding_cases = seeding_cases,
                     seeding_age_order = seeding_age_order,
                     prob_hosp = prob_hosp,
                     tt_prob_hosp_multiplier = tt_prob_hosp_multiplier,
                     prob_hosp_multiplier = prob_hosp_multiplier,
                     prob_severe = prob_severe,
                     prob_severe_multiplier = prob_severe_multiplier,
                     tt_prob_severe_multiplier = tt_prob_severe_multiplier,
                     prob_non_severe_death_treatment = prob_non_severe_death_treatment,
                     prob_non_severe_death_no_treatment = prob_non_severe_death_no_treatment,
                     prob_severe_death_treatment = prob_severe_death_treatment,
                     prob_severe_death_no_treatment = prob_severe_death_no_treatment,
                     p_dist = p_dist,
                     rel_infectiousness = rel_infectiousness,
                     rel_infectiousness_vaccinated = rel_infectiousness_vaccinated,
                     dur_E = dur_E,
                     tt_dur_E = tt_dur_E,
                     dur_IMild = dur_IMild,
                     tt_dur_IMild = tt_dur_IMild,
                     dur_ICase = dur_ICase,
                     tt_dur_ICase = tt_dur_ICase,
                     dur_get_ox_survive = dur_get_ox_survive,
                     tt_dur_get_ox_survive = tt_dur_get_ox_survive,
                     dur_get_ox_die = dur_get_ox_die,
                     tt_dur_get_ox_die = tt_dur_get_ox_die,
                     dur_not_get_ox_survive = dur_not_get_ox_survive,
                     dur_not_get_ox_die = dur_not_get_ox_die,
                     dur_get_mv_survive = dur_get_mv_survive,
                     tt_dur_get_mv_survive = tt_dur_get_mv_survive,
                     dur_get_mv_die = dur_get_mv_die,
                     tt_dur_get_mv_die = tt_dur_get_mv_die,
                     dur_not_get_mv_survive = dur_not_get_mv_survive,
                     dur_not_get_mv_die = dur_not_get_mv_die,
                     dur_rec = dur_rec,
                     dur_R = dur_R,
                     tt_dur_R = tt_dur_R,
                     hosp_bed_capacity = hosp_bed_capacity,
                     ICU_bed_capacity = ICU_bed_capacity,
                     tt_hosp_beds = tt_hosp_beds,
                     tt_ICU_beds = tt_ICU_beds,
                     dur_V = dur_V,
                     tt_dur_V = tt_dur_V,
                     vaccine_efficacy_infection = vaccine_efficacy_infection,
                     tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection,
                     vaccine_efficacy_disease = vaccine_efficacy_disease,
                     tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
                     primary_doses = primary_doses,
                     tt_primary_doses = tt_primary_doses,
                     second_dose_delay = second_dose_delay,
                     booster_doses = booster_doses,
                     tt_booster_doses = tt_booster_doses,
                     protection_delay_rate = protection_delay_rate,
                     protection_delay_shape = protection_delay_shape,
                     protection_delay_time = time_period,
                     vaccine_coverage_mat = vaccine_coverage_mat,
                     vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
                     init = init)

  # Set model type
  replicates <- 1
  mod_gen = nimue_booster

  # Running the Model
  mod <- mod_gen$new(user = pars, unused_user_action = "ignore",
                     use_dde = use_dde)

  # Daily output by default
  t <- round(seq(from = 1, to = time_period))

  results <- mod$run(t, ...)

  # coerce to array
  results <- array(results, dim = c(dim(results), 1), dimnames = dimnames(results))

  # Summarise inputs
  parameters <- args
  parameters$population <- pars$population
  parameters$hosp_bed_capacity <- pars$hosp_beds
  parameters$ICU_bed_capacity <- pars$ICU_beds
  parameters$beta_set <- pars$beta_set
  parameters$seeding_cases <- pars$E1_0
  parameters$contact_matrix_set <- pars$contact_matrix_set

  out <- list(output = results, parameters = parameters, model = mod, odin_parameters = pars)
  out <- structure(out, class = c("lmic_booster_nimue_simulation", "nimue_simulation"))
  return(out)

}
