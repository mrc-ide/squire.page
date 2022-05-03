#'
#'@export
pmcmc_drjacoby <- function(data,
                           replicates = 100,
                          n_mcmc,
                          n_burnin,
                          date_vaccine_change,
                          first_doses,
                          second_doses,
                          booster_doses,
                          n_chains = 1,
                          log_likelihood = NULL,
                          log_prior = NULL,
                          drjacoby_list = list(),
                          squire_model = nimue_booster_model(),
                          pars_obs = list(phi_cases = 1,
                                          k_cases = 2,
                                          phi_death = 1,
                                          k_death = 2,
                                          exp_noise = 1e6,
                                          cases_fitting = TRUE,
                                          cases_days = 21,
                                          cases_reporting = 21,
                                          variant_adjust = NULL),
                          pars_init = list('start_date'     = as.Date("2020-02-07"),
                                           'R0'             = 2.5,
                                           'Meff'           = 2,
                                           'Meff_pl'        = 3),
                          pars_min = list('start_date'      = as.Date("2020-02-01"),
                                          'R0'              = 0,
                                          'Meff'            = 1,
                                          'Meff_pl'         = 2),
                          pars_max = list('start_date'      = as.Date("2020-02-20"),
                                          'R0'              = 5,
                                          'Meff'            = 3,
                                          'Meff_pl'         = 4),
                          pars_discrete = list('start_date' = TRUE,
                                               'R0'         = FALSE,
                                               'Meff'       = FALSE,
                                               'Meff_pl'    = FALSE),
                          country = NULL,
                          population = NULL,
                          contact_matrix_set = NULL,
                          baseline_contact_matrix = NULL,
                          date_contact_matrix_set_change = NULL,
                          R0_change = NULL,
                          date_R0_change = NULL,
                          hosp_bed_capacity = NULL,
                          baseline_hosp_bed_capacity = NULL,
                          date_hosp_bed_capacity_change = NULL,
                          ICU_bed_capacity = NULL,
                          baseline_ICU_bed_capacity = NULL,
                          date_ICU_bed_capacity_change = NULL,
                          baseline_first_doses = 0,
                          baseline_second_doses = 0,
                          baseline_booster_doses = 0,
                          baseline_vaccine_efficacy_infection = vaccine_pars_booster$vaccine_efficacy_infection,
                          vaccine_efficacy_infection = NULL,
                          baseline_vaccine_efficacy_disease = vaccine_pars_booster$vaccine_efficacy_disease,
                          vaccine_efficacy_disease = NULL,
                          date_vaccine_efficacy_change = NULL,
                          Rt_args = NULL,
                          ...
) {

  #------------------------------------------------------------
  # Section 1 of pMCMC Wrapper: Checks & Setup
  #------------------------------------------------------------

  #--------------------
  # assertions & checks
  #--------------------

  squire:::check_drjacoby_list(drjacoby_list)

  # we work with pars_init being a list of inital conditions for starting
  if(any(c("start_date", "R0") %in% names(pars_init))) {
    pars_init <- list(pars_init)
  }

  # make it same length as chains, which allows us to pass in multiple starting points
  if(length(pars_init) != n_chains) {
    pars_init <- rep(pars_init, n_chains)
    pars_init <- pars_init[seq_len(n_chains)]
  }

  # data assertions
  squire:::assert_dataframe(data)
  squire:::assert_in("date", names(data))
  squire:::assert_in("deaths", names(data))
  squire:::assert_date(data$date)
  squire:::assert_increasing(as.numeric(as.Date(data$date)),
                             message = "Dates must be in increasing order")

  # check input pars df
  squire:::assert_list(pars_init)
  squire:::assert_list(pars_init[[1]])
  squire:::assert_list(pars_min)
  squire:::assert_list(pars_max)
  squire:::assert_list(pars_discrete)
  squire:::assert_eq(names(pars_init[[1]]), names(pars_min))
  squire:::assert_eq(names(pars_min), names(pars_max))
  squire:::assert_eq(names(pars_max), names(pars_discrete))
  squire:::assert_in(c("R0", "start_date"),names(pars_init[[1]]),
                     message = "Params to infer must include R0, start_date")
  squire:::assert_date(pars_init[[1]]$start_date)
  squire:::assert_date(pars_min$start_date)
  squire:::assert_date(pars_max$start_date)
  if (pars_max$start_date >= as.Date(data$date[1])-1) {
    stop("Maximum start date must be at least 2 days before the first date in data")
  }

  # check date variables are as Date class
  for(i in seq_along(pars_init)) {
    pars_init[[i]]$start_date <- as.Date(pars_init[[i]]$start_date)
  }
  pars_min$start_date <- as.Date(pars_min$start_date)
  pars_max$start_date <- as.Date(pars_max$start_date)

  # check bounds
  for(var in names(pars_init[[1]])) {

    squire:::assert_bounded(as.numeric(pars_init[[1]][[var]]),
                            left = as.numeric(pars_min[[var]]),
                            right = as.numeric(pars_max[[var]]),
                            name = paste(var, "init"))

    squire:::assert_single_numeric(as.numeric(pars_min[[var]]), name = paste(var, "min"))
    squire:::assert_single_numeric(as.numeric(pars_max[[var]]), name = paste(var, "max"))
    squire:::assert_single_numeric(as.numeric(pars_init[[1]][[var]]), name = paste(var, "init"))

  }

  # additonal checks that R0 is positive as undefined otherwise
  squire:::assert_pos(pars_min$R0)
  squire:::assert_pos(pars_max$R0)
  squire:::assert_pos(pars_init[[1]]$R0)
  squire:::assert_bounded(pars_init[[1]]$R0, left = pars_min$R0, right = pars_max$R0)

  # check likelihood items
  if ( !(is.null(log_likelihood) | inherits(log_likelihood, "function")) ) {
    stop("Log Likelihood (log_likelihood) must be null or a user specified function")
  }
  if ( !(is.null(log_prior) | inherits(log_prior, "function")) ) {
    stop("Log Likelihood (log_likelihood) must be null or a user specified function")
  }
  squire:::assert_logical(unlist(pars_discrete))
  squire:::assert_list(pars_obs)
  squire:::assert_in(c("phi_cases", "k_cases", "phi_death", "k_death", "exp_noise"), names(pars_obs))
  squire:::assert_numeric(unlist(pars_obs[c("phi_cases", "k_cases", "phi_death", "k_death", "exp_noise")]))

  # mcmc items
  squire:::assert_pos_int(n_mcmc)
  squire:::assert_pos_int(n_chains)
  squire:::assert_int(n_burnin)

  # squire and odin
  squire:::assert_custom_class(squire_model, "squire_model")
  squire:::assert_pos_int(replicates)

  # date change items
  squire:::assert_same_length(R0_change, date_R0_change)
  # checks that dates are not in the future compared to our data
  if (!is.null(date_R0_change)) {
    squire:::assert_date(date_R0_change)
    if(as.Date(utils::tail(date_R0_change,1)) > as.Date(utils::tail(data$date, 1))) {
      stop("Last date in date_R0_change is greater than the last date in data")
    }
  }

  # ------------------------------------
  # checks on odin interacting variables
  # ------------------------------------

  if(!is.null(contact_matrix_set)) {
    squire:::assert_list(contact_matrix_set)
  }
  squire:::assert_same_length(contact_matrix_set, date_contact_matrix_set_change)
  squire:::assert_same_length(ICU_bed_capacity, date_ICU_bed_capacity_change)
  squire:::assert_same_length(hosp_bed_capacity, date_hosp_bed_capacity_change)
  squire:::assert_same_length(first_doses, date_vaccine_change)
  squire:::assert_same_length(second_doses, date_vaccine_change)
  squire:::assert_same_length(booster_doses, date_vaccine_change)
  squire:::assert_same_length(vaccine_efficacy_infection, date_vaccine_efficacy_change)
  squire:::assert_same_length(vaccine_efficacy_disease, date_vaccine_efficacy_change)

  # handle contact matrix changes
  if(!is.null(date_contact_matrix_set_change)) {

    squire:::assert_date(date_contact_matrix_set_change)
    squire:::assert_list(contact_matrix_set)

    if(is.null(baseline_contact_matrix)) {
      stop("baseline_contact_matrix can't be NULL if date_contact_matrix_set_change is provided")
    }
    if(as.Date(utils::tail(date_contact_matrix_set_change,1)) > as.Date(utils::tail(data$date, 1))) {
      stop("Last date in date_contact_matrix_set_change is greater than the last date in data")
    }

    # Get in correct format
    if(is.matrix(baseline_contact_matrix)) {
      baseline_contact_matrix <- list(baseline_contact_matrix)
    }

    tt_contact_matrix <- c(0, seq_len(length(date_contact_matrix_set_change)))
    contact_matrix_set <- append(baseline_contact_matrix, contact_matrix_set)

  } else {
    tt_contact_matrix <- 0
    contact_matrix_set <- baseline_contact_matrix
  }

  # handle ICU changes
  if(!is.null(date_ICU_bed_capacity_change)) {

    squire:::assert_date(date_ICU_bed_capacity_change)
    squire:::assert_vector(ICU_bed_capacity)
    squire:::assert_numeric(ICU_bed_capacity)

    if(is.null(baseline_ICU_bed_capacity)) {
      stop("baseline_ICU_bed_capacity can't be NULL if date_ICU_bed_capacity_change is provided")
    }
    squire:::assert_numeric(baseline_ICU_bed_capacity)
    if(as.Date(utils::tail(date_ICU_bed_capacity_change,1)) > as.Date(utils::tail(data$date, 1))) {
      stop("Last date in date_ICU_bed_capacity_change is greater than the last date in data")
    }

    tt_ICU_beds <- c(0, seq_len(length(date_ICU_bed_capacity_change)))
    ICU_bed_capacity <- c(baseline_ICU_bed_capacity, ICU_bed_capacity)

  } else {
    tt_ICU_beds <- 0
    ICU_bed_capacity <- baseline_ICU_bed_capacity
  }

  # handle vaccine changes
  if(!is.null(date_vaccine_change)) {

    squire:::assert_date(date_vaccine_change)
    squire:::assert_vector(first_doses)
    squire:::assert_numeric(first_doses)
    squire:::assert_numeric(baseline_first_doses)

    if(is.null(baseline_first_doses)) {
      stop("baseline_first_doses can't be NULL if date_vaccine_change is provided")
    }
    if(as.Date(utils::tail(date_vaccine_change,1)) > as.Date(utils::tail(data$date, 1))) {
      stop("Last date in date_vaccine_change is greater than the last date in data")
    }

    tt_first_doses <- c(0, seq_len(length(date_vaccine_change)))
    first_doses <- c(baseline_first_doses, first_doses)
    tt_second_doses <- tt_first_doses
    second_doses <- c(baseline_second_doses, second_doses)
    tt_booster_doses <- tt_first_doses
    booster_doses <- c(baseline_booster_doses, booster_doses)

  } else {
    tt_vaccine <- 0
    if(!is.null(baseline_first_doses)) {
      first_doses <- baseline_first_doses
    } else {
      first_doses <- 0
    }
    if(!is.null(baseline_second_doses)) {
      second_doses <- baseline_second_doses
    } else {
      second_doses <- 0
    }
    if(!is.null(baseline_booster_doses)) {
      booster_doses <- baseline_booster_doses
    } else {
      booster_doses <- 0
    }
  }

  # handle vaccine efficacy disease changes
  if(!is.null(date_vaccine_efficacy_change)) {

    squire:::assert_date(date_vaccine_efficacy_change)
    if(!is.list(vaccine_efficacy_infection)) {
      vaccine_efficacy_infection <- list(vaccine_efficacy_infection)
    }
    squire:::assert_vector(vaccine_efficacy_infection[[1]])
    squire:::assert_numeric(vaccine_efficacy_infection[[1]])
    squire:::assert_numeric(baseline_vaccine_efficacy_infection)

    if(is.null(baseline_vaccine_efficacy_infection)) {
      stop("baseline_vaccine_efficacy_infection can't be NULL if date_vaccine_efficacy_change is provided")
    }
    if(as.Date(utils::tail(date_vaccine_efficacy_change,1)) > as.Date(utils::tail(data$date, 1))) {
      stop("Last date in date_vaccine_efficacy_change is greater than the last date in data")
    }

    tt_vaccine_efficacy_infection <- c(0, seq_len(length(date_vaccine_efficacy_change)))
    vaccine_efficacy_infection <- c(list(baseline_vaccine_efficacy_infection), vaccine_efficacy_infection)

  } else {
    tt_vaccine_efficacy_infection <- 0
    if(!is.null(baseline_vaccine_efficacy_infection)) {
      vaccine_efficacy_infection <- baseline_vaccine_efficacy_infection
    } else {
      vaccine_efficacy_infection <- rep(0.8, 17)
    }
  }

  # handle vaccine efficacy disease changes
  if(!is.null(date_vaccine_efficacy_change)) {

    squire:::assert_date(date_vaccine_efficacy_change)
    if(!is.list(vaccine_efficacy_disease)) {
      vaccine_efficacy_disease <- list(vaccine_efficacy_disease)
    }
    squire:::assert_vector(vaccine_efficacy_disease[[1]])
    squire:::assert_numeric(vaccine_efficacy_disease[[1]])
    squire:::assert_numeric(baseline_vaccine_efficacy_disease)

    if(is.null(baseline_vaccine_efficacy_disease)) {
      stop("baseline_vaccine_efficacy_disease can't be NULL if date_vaccine_efficacy_change is provided")
    }
    if(as.Date(utils::tail(date_vaccine_efficacy_change,1)) > as.Date(utils::tail(data$date, 1))) {
      stop("Last date in date_vaccine_efficacy_change is greater than the last date in data")
    }

    tt_vaccine_efficacy_disease <- c(0, seq_len(length(date_vaccine_efficacy_change)))
    vaccine_efficacy_disease <- c(list(baseline_vaccine_efficacy_disease), vaccine_efficacy_disease)

  } else {
    tt_vaccine_efficacy_disease <- 0
    if(!is.null(baseline_vaccine_efficacy_disease)) {
      vaccine_efficacy_disease <- baseline_vaccine_efficacy_disease
    } else {
      vaccine_efficacy_disease <- rep(0.95, 17)
    }
  }


  # handle hosp bed changed
  if(!is.null(date_hosp_bed_capacity_change)) {

    squire:::assert_date(date_hosp_bed_capacity_change)
    squire:::assert_vector(hosp_bed_capacity)
    squire:::assert_numeric(hosp_bed_capacity)

    if(is.null(baseline_hosp_bed_capacity)) {
      stop("baseline_hosp_bed_capacity can't be NULL if date_hosp_bed_capacity_change is provided")
    }
    squire:::assert_numeric(baseline_hosp_bed_capacity)
    if(as.Date(utils::tail(date_hosp_bed_capacity_change,1)) > as.Date(utils::tail(data$date, 1))) {
      stop("Last date in date_hosp_bed_capacity_change is greater than the last date in data")
    }

    tt_hosp_beds <- c(0, seq_len(length(date_hosp_bed_capacity_change)))
    hosp_bed_capacity <- c(baseline_hosp_bed_capacity, hosp_bed_capacity)

  } else {
    tt_hosp_beds <- 0
    hosp_bed_capacity <- baseline_hosp_bed_capacity
  }

  #----------------
  # Generate Odin items
  #----------------

  # make the date definitely a date
  data$date <- as.Date(as.character(data$date))

  # build model parameters
  model_params <- squire_model$parameter_func(
    country = country,
    population = population,
    dt = 1,
    contact_matrix_set = contact_matrix_set,
    tt_contact_matrix = tt_contact_matrix,
    hosp_bed_capacity = hosp_bed_capacity,
    tt_hosp_beds = tt_hosp_beds,
    ICU_bed_capacity = ICU_bed_capacity,
    tt_ICU_beds = tt_ICU_beds,
    first_doses = first_doses,
    tt_first_doses = tt_first_doses,
    second_doses = second_doses,
    tt_second_doses = tt_second_doses,
    booster_doses = booster_doses,
    tt_booster_doses = tt_booster_doses,
    vaccine_efficacy_infection = vaccine_efficacy_infection,
    tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection,
    vaccine_efficacy_disease = vaccine_efficacy_disease,
    tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
    ...)

  # collect interventions for odin model likelihood
  interventions <- list(
    R0_change = R0_change,
    date_R0_change = date_R0_change,
    date_contact_matrix_set_change = date_contact_matrix_set_change,
    contact_matrix_set = contact_matrix_set,
    date_ICU_bed_capacity_change = date_ICU_bed_capacity_change,
    ICU_bed_capacity = ICU_bed_capacity,
    date_hosp_bed_capacity_change = date_hosp_bed_capacity_change,
    hosp_bed_capacity = hosp_bed_capacity,
    date_vaccine_change = date_vaccine_change,
    first_doses = first_doses,
    second_doses = second_doses,
    booster_doses = booster_doses,
    date_vaccine_efficacy_change = date_vaccine_efficacy_change,
    vaccine_efficacy_disease = vaccine_efficacy_disease,
    date_vaccine_efficacy_change = date_vaccine_efficacy_change,
    vaccine_efficacy_infection = vaccine_efficacy_infection
  )

  #----------------..
  # Collect Odin and MCMC Inputs
  #----------------..
  inputs <- list(
    data = data,
    n_mcmc = n_mcmc,
    model_params = model_params,
    interventions = interventions,
    pars_obs = pars_obs,
    Rt_args = Rt_args,
    squire_model = squire_model,
    pars = list(pars_obs = pars_obs,
                pars_init = pars_init,
                pars_min = pars_min,
                pars_max = pars_max,
                pars_discrete = pars_discrete))


  #----------------
  # create prior and likelihood functions given the inputs
  #----------------

  if(is.null(log_prior)) {
    # set improper, uninformative prior
    log_prior <- function(params, misc) log(1e-10)
  }
  calc_lprior <- log_prior
  squire:::check_drjacoby_logprior(calc_lprior)

  if(is.null(log_likelihood)) {
    log_likelihood <- squire:::calc_loglikelihood_drjacoby()
  }
  calc_ll <- log_likelihood
  squire:::check_drjacoby_loglike(calc_ll)

  #----------------
  # set run_mcmc_func here to out drjacoby one
  #----------------
  run_mcmc_func <- squire:::run_drjacoby_mcmc

  # needs to be a vector to pass to reflecting boundary function
  pars_min <- unlist(pars_min)
  pars_max <- unlist(pars_max)
  pars_discrete <- unlist(pars_discrete)

  #--------------------------------------------------------
  # Section 2 of pMCMC Wrapper: Run pMCMC
  #--------------------------------------------------------

  # Are we debuggine
  if (Sys.getenv("SQUIRE_PARALLEL_DEBUG") == "TRUE") {

    # if debug remove the cluster
    drjacoby_list$cl <- NULL

  }

  # run drjacoby
  message("Running drjacoby...")
  mcmc_out <- squire:::run_drjacoby_mcmc(loglike = calc_ll,
                                         logprior = calc_lprior,
                                         inputs = inputs,
                                         burnin = n_burnin,
                                         chains = n_chains,
                                         drjacoby_list = drjacoby_list)

  # process output to play with rest of squire
  chains <- squire:::convert_drjacoby_mcmc(mcmc_out)
  if(n_chains > 1) {
    pmcmc <- list(inputs = inputs,
                  chains = chains,
                  drjacoby_out = mcmc_out)

    pmcmc$inputs$n_mcmc <- n_mcmc
    class(pmcmc) <- 'squire_pmcmc_list'
  } else {
    pmcmc <- chains$chain1
    pmcmc$inputs <- inputs
    pmcmc$drjacoby_out <- mcmc_out
    class(pmcmc) <- 'squire_pmcmc'
  }

  #--------------------------------------------------------
  # Section 3 of pMCMC Wrapper: Sample PMCMC Results
  #--------------------------------------------------------
  pmcmc_samples <- squire::sample_drjacoby(pmcmc_results = pmcmc,
                                           burnin = n_burnin,
                                           n_chains = n_chains,
                                           n_trajectories = replicates,
                                           log_likelihood = log_likelihood,
                                           forecast_days = 0)

  #--------------------------------------------------------
  # Section 4 of pMCMC Wrapper: Tidy Output
  #--------------------------------------------------------

  #----------------
  # Pull Sampled results and "recreate" squire models
  #----------------
  # create a fake run object and fill in the required elements
  r <- squire_model$run_func(country = country,
                             contact_matrix_set = contact_matrix_set,
                             tt_contact_matrix = tt_contact_matrix,
                             hosp_bed_capacity = hosp_bed_capacity,
                             tt_hosp_beds = tt_hosp_beds,
                             ICU_bed_capacity = ICU_bed_capacity,
                             tt_ICU_beds = tt_ICU_beds,
                             first_doses = first_doses,
                             tt_first_doses = tt_first_doses,
                             second_doses = second_doses,
                             tt_second_doses = tt_second_doses,
                             booster_doses = booster_doses,
                             tt_booster_doses = tt_booster_doses,
                             vaccine_efficacy_infection = vaccine_efficacy_infection,
                             tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection,
                             vaccine_efficacy_disease = vaccine_efficacy_disease,
                             tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
                             population = population,
                             replicates = 1,
                             day_return = TRUE,
                             time_period = nrow(pmcmc_samples$trajectories),
                             ...)

  # and add the parameters that changed between each simulation, i.e. posterior draws
  r$replicate_parameters <- pmcmc_samples$sampled_PMCMC_Results

  # as well as adding the pmcmc chains so it's easy to draw from the chains again in the future
  r$pmcmc_results <- pmcmc

  # then let's create the output that we are going to use
  names(pmcmc_samples)[names(pmcmc_samples) == "trajectories"] <- "output"
  dimnames(pmcmc_samples$output) <- list(dimnames(pmcmc_samples$output)[[1]], dimnames(r$output)[[2]], NULL)
  r$output <- pmcmc_samples$output

  # and adjust the time as before
  full_row <- match(0, apply(r$output[,"time",],2,function(x) { sum(is.na(x)) }))
  saved_full <- r$output[,"time",full_row]
  for(i in seq_len(replicates)) {
    na_pos <- which(is.na(r$output[,"time",i]))
    full_to_place <- saved_full - which(rownames(r$output) == as.Date(max(data$date))) + 1L
    if(length(na_pos) > 0) {
      full_to_place[na_pos] <- NA
    }
    r$output[,"time",i] <- full_to_place
  }

  # second let's recreate the output
  r$model <- pmcmc_samples$inputs$squire_model$odin_model(
    user = pmcmc_samples$inputs$model_params, unused_user_action = "ignore"
  )
  # we will add the interventions here so that we know what times are needed for projection
  r$interventions <- interventions

  # and fix the replicates
  r$parameters$replicates <- replicates
  r$parameters$time_period <- as.numeric(diff(as.Date(range(rownames(r$output)))))
  r$parameters$dt <- model_params$dt

  #--------------------..
  # out
  #--------------------..
  return(r)
}
