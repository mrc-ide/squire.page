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
  booster_data <- nimue_format(model_out, c("vaccinated_first_dose", "vaccinated_second_dose", "vaccinated_first_waned", "vaccinated_second_waned"))
  expect_equal(sort(as.character(unique(booster_data$compartment))), c("vaccinated_first_dose", "vaccinated_first_waned", "vaccinated_second_dose",
                                                                       "vaccinated_second_waned"))
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
                              first_doses = 4000,
                              second_doses = 3000,
                              booster_doses = 2000,
                              dur_R = 365, time_period = 3*365)
  #compare to no boosters
  no_booster <-
    squire.page:::run_booster(country = "United Kingdom",
                              first_doses = 4000,
                              second_doses = 3000,
                              booster_doses = 0,
                              dur_R = 365, time_period = 3*365)
  expect_gt(get_total(no_booster, "deaths"), get_total(base, "deaths"))
  expect_gt(get_total(no_booster, "vaccinated_second_waned"), get_total(base, "vaccinated_second_waned"))
  #no second doses
  no_second <-
    squire.page:::run_booster(country = "United Kingdom",
                              first_doses = 4000,
                              second_doses = 0,
                              booster_doses = 2000,
                              dur_R = 365, time_period = 3*365)
  expect_gt(get_total(no_second, "deaths"), get_total(no_booster, "deaths"))
  expect_equal(get_total(no_second, "vaccinated_second_waned"), 0)
  expect_gt(get_total(no_second, "vaccinated_first_waned"), get_total(no_booster, "vaccinated_first_waned"))
  #no first doses
  no_first <-
    squire.page:::run_booster(country = "United Kingdom",
                              first_doses = 0,
                              second_doses = 3000,
                              booster_doses = 2000,
                              dur_R = 365, time_period = 3*365)
  expect_gt(get_total(no_first, "deaths"), get_total(no_booster, "deaths"))
  expect_equal(get_total(no_first, "vaccinated_first_waned"), 0)
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
  expect_true(all(nimue_format(base, "vaccinated_first_dose") %>% pull(y) %>% diff() > 0))
  expect_true(all(nimue_format(base, "vaccinated_second_dose") %>% pull(y) %>% diff() > 0))
  expect_true(all(nimue_format(base, "vaccinated_second_waned") %>% pull(y) %>% diff() > 0))
  #doses make sense
  expect_true(all(nimue_format(base, "first_doses_given") %>% pull(y) <= 4001))
  expect_true(all(nimue_format(base, "second_doses_given") %>% pull(y) <= 3001))
  expect_true(all(nimue_format(base, "booster_doses_given") %>% pull(y) <= 2001))
  #check vaccination rollout in zero infection epidemic
  zero_inf <-
    squire.page:::run_booster(country = "United Kingdom",
                              first_doses = 4000,
                              second_doses = 3000,
                              booster_doses = 2000,
                              dur_R = 365, time_period = 3*365,
                              beta = 0)
  expect_true(all(nimue_format(zero_inf, "vaccinated_first_dose") %>% pull(y) %>% diff() %>% `-`(1000) %>% abs() %>% tail(3*365 - 10) < 1))
  expect_true(all(nimue_format(zero_inf, "vaccinated_second_dose") %>% pull(y) %>% diff() %>% `-`(3000) %>% abs() %>% tail(3*365 - 20) < 1))
  expect_true(all(nimue_format(zero_inf, "booster_doses_given") %>% pull(y) %>% `-`(2000) %>% abs() %>% tail(3*365 - 100) < 1))
})

test_that("LMIC Booster Likelihood Function", {
  set.seed(100)
  #just see if it works
  pmcmc_output <- suppressMessages(pmcmc_booster_drjacoby(
    data = data.frame(
      date = seq(as.Date("2020-03-01"), as.Date("2020-03-01") + 99, by = 1),
      deaths = rpois(100, 50),
      cases = rpois(100, 100)
    ),
    replicates = 10,
    n_mcmc = 100,
    n_burnin = 50,
    quasi_likelihood = FALSE,
    country = "United Kingdom",
    first_doses = 4000,
    second_doses = 3000,
    booster_doses = 2000,
    date_vaccine_change = "2020-05-01",
    R0_change = 1,
    date_R0_change = "2020-03-01"
  ))

  expect_s3_class(pmcmc_output, "lmic_booster_nimue_simulation")
  expect_true(nrow(pmcmc_output$replicate_parameters) == 10)

  #check format works on this
  output <- nimue_format(pmcmc_output, "deaths")
  expect_true(length(unique(output$replicate)) == 10)

})

test_that("LMIC Booster Quasi-Likelihood Function", {
  set.seed(100)
  #just see if it works
  pmcmc_output <- suppressMessages(pmcmc_booster_drjacoby(
    data = data.frame(
      date = seq(as.Date("2020-03-01"), as.Date("2020-03-01") + 99, by = 1),
      deaths = rpois(100, 50),
      cases = rpois(100, 100)
    ),
    replicates = 10,
    n_mcmc = 100,
    n_burnin = 50,
    quasi_likelihood = TRUE,
    country = "United Kingdom",
    first_doses = 4000,
    second_doses = 3000,
    booster_doses = 2000,
    date_vaccine_change = "2020-05-01",
    R0_change = 1,
    date_R0_change = "2020-03-01"
  ))

  expect_s3_class(pmcmc_output, "lmic_booster_nimue_simulation")
  expect_true(nrow(pmcmc_output$replicate_parameters) == 10)

  #check format works on this
  output <- nimue_format(pmcmc_output, "deaths")
  expect_true(length(unique(output$replicate)) == 10)

})
