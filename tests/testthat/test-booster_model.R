library(dplyr)
library(tidyr)
test_that("Run function", {
  #just see if it works
  model_out <-
    squire.page:::run_booster(country = "United Kingdom")
  expect_s3_class(model_out, "lmic_booster_nimue_simulation")
  expect_s3_class(model_out, "nimue_simulation")
  model_out <-
    squire.page:::run_booster(country = "Gambia")
  expect_true(!is.null(model_out$output))
  expect_equal(names(model_out), c("output", "parameters", "model", "odin_parameters"))
})

test_that("Format", {
  model_out <-
    squire.page:::run_booster(country = "United Kingdom")
  #functionality
  deaths <- nimue_format(model_out, "deaths")
  expect_equal(as.character(unique(deaths$compartment)), "deaths")
  expect_true(is.numeric(deaths$y))
  #booster specific
  booster_data <- nimue_format(model_out, c("vaccinated_first_dose", "vaccinated_second_dose", "vaccinated_booster_dose", "vaccinated_second_waned", "vaccinated_booster_waned"))
  expect_equal(sort(as.character(unique(booster_data$compartment))), c("vaccinated_booster_dose", "vaccinated_booster_waned", "vaccinated_first_dose",
                                                                       "vaccinated_second_dose", "vaccinated_second_waned"))
  expect_true(is.numeric(booster_data$y))
  #warning if not available
  expect_warning(
    nimue_format(model_out, c("vaccinated", "deaths"))
  )
})

test_that("Model Functionality", {
  get_total <- function(out, measure){
    out %>% nimue_format(measure) %>% pull(y) %>% sum(na.rm = TRUE)
  }
  base <-
    squire.page:::run_booster(country = "United Kingdom",
                              primary_doses = 4000,
                              booster_doses = 2000,
                              dur_R = 365, time_period = 3*365)
  #compare to no boosters
  no_booster <-
    squire.page:::run_booster(country = "United Kingdom",
                              primary_doses = 4000,
                              booster_doses = 0,
                              dur_R = 365, time_period = 3*365)
  expect_gt(get_total(no_booster, "deaths"), get_total(base, "deaths"))
  expect_gt(get_total(no_booster, "vaccinated_second_waned"), get_total(base, "vaccinated_second_waned"))
  #no first doses
  no_first <-
    squire.page:::run_booster(country = "United Kingdom",
                              primary_doses = 0,
                              booster_doses = 2000,
                              dur_R = 365, time_period = 3*365)
  expect_gt(get_total(no_first, "deaths"), get_total(no_booster, "deaths"))
  expect_equal(get_total(no_first, "vaccinated_booster_dose"), 0)
  expect_equal(get_total(no_first, "vaccinated_booster_waned"), 0)
  expect_equal(get_total(no_first, "vaccinated_second_dose"), 0)
  expect_equal(get_total(no_first, "vaccinated_second_waned"), 0)
  expect_equal(get_total(no_first, "vaccinated_first_dose"), 0)
  #no population growth
  expect_equal(nimue_format(base, c("S", "E",
                       "IMild", "ICase", "IICU", "IHospital",
                       "IRec", "R", "D")) %>%
    group_by(t) %>%
    summarise(N = sum(y, na.rm = TRUE)) %>%
    pull(N), nimue_format(base, "N") %>% pull(y))
  expect_true(nimue_format(base, c("S", "E",
                                   "IMild", "ICase", "IICU", "IHospital",
                                   "IRec", "R", "D")) %>%
                group_by(t) %>%
                summarise(N = sum(y, na.rm = TRUE)) %>%
                pull(N) %>%
                diff() %>%
                abs() %>%
                `<`(1) %>%
                all())
  #with these numbers expect doses to be increasing all the time
  expect_true(all(nimue_format(base, "vaccinated_first_dose") %>% pull(y) %>% diff() >= 0))
  expect_true(all(nimue_format(base, "vaccinated_second_dose") %>% pull(y) %>% diff() >= 0))
  #allow error of a maximum of five doses
  expect_true(all(nimue_format(base, "vaccinated_booster_dose") %>% pull(y) %>% diff() %>% `+`(5) >= 0))
  #doses make sense
  expect_true(all(nimue_format(base, "first_doses_given") %>% pull(y) <= 4005))
  expect_true(all(nimue_format(base, "second_doses_given") %>% pull(y) <= 4005))
  expect_true(all(nimue_format(base, "booster_doses_given") %>% pull(y) <= 2005))
  # #check vaccination rollout in zero infection epidemic
  # zero_inf <-
  #   squire.page:::run_booster(country = "United Kingdom",
  #                             primary_doses = 4000,
  #                             booster_doses = 2000,
  #                             dur_R = 365, time_period = 3*365,
  #                             beta = 0)
})

test_that("LMIC Booster Likelihood Function", {
  set.seed(100)
  #just see if it works
  data <- data.frame(
    date = seq(as.Date("2020-03-01"), as.Date("2020-03-01") + 99, by = 1),
    deaths = rpois(100, 50),
    cases = rpois(100, 100)
  )
  pmcmc_output <- suppressMessages(pmcmc_drjacoby(
    data = data,
    replicates = 3,
    n_mcmc = 5,
    n_burnin = 5,
    log_likelihood = squire:::convert_log_likelihood_func_for_drjacoby(calc_loglikelihood_booster),
    country = "United Kingdom",
    primary_doses = 4000,
    booster_doses = 2000,
    date_vaccine_change = "2020-05-01",
    R0_change = 1,
    date_R0_change = "2020-03-01",
    second_dose_delay = 60,
    protection_delay_rate = NULL,
    protection_delay_shape = NULL,
    drjacoby_list = list(silent = TRUE)
  ))

  expect_s3_class(pmcmc_output, "lmic_booster_nimue_simulation")
  expect_true(nrow(pmcmc_output$replicate_parameters) == 3)

  #check format works on this
  output <- nimue_format(pmcmc_output, "deaths")
  expect_true(length(unique(output$replicate)) == 3)
})
