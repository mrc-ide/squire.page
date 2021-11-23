test_that("Correct common outputs for standard fit", {
  #test using Afghanistan fit
  df <- prepare_input_json_df(afg_fit, ox_interventions = NULL)
  expect_true("Rt" %in% names(df))
  expect_true("tt_beta" %in% names(df))
  expect_true("date" %in% names(df))
  expect_true("beta_set" %in% names(df))
  expect_true("beta_set_min" %in% names(df))
  expect_true("beta_set_max" %in% names(df))
  expect_true("Rt_min" %in% names(df))
  expect_true("Rt_max" %in% names(df))
  expect_true("deaths" %in% names(df))
  expect_type(df$date, "double")
})
test_that("Correct vaccine outputs for standard fit", {
  #test using Afghanistan fit
  df <- prepare_input_json_df(afg_fit, ox_interventions = NULL)
  expect_true("vaccine_efficacy_disease" %in% names(df))
  expect_true("vaccine_strategy" %in% names(df))
  expect_true("vaccine_coverage" %in% names(df))
  expect_true("vaccine_efficacy_infection" %in% names(df))
  expect_true("max_vaccine" %in% names(df))
  expect_equal(unique(df$vaccines_available), 0.2)
  expect_equal(unique(df$vaccine_coverage), 0.8)
})
test_that("Extension correct", {
  #test using Afghanistan fit
  df <- prepare_input_json_df(afg_fit, ox_interventions = NULL)
  expect_true(
    nrow(afg_fit$pmcmc_results$inputs$data) + 240 <=
    nrow(df)
  )
})
test_that("Errors", {
  #test using Afghanistan fit
  temp <- afg_fit
  temp$pmcmc_results <- NULL
  expect_error(prepare_input_json_df(temp, ox_interventions = NULL))
  temp <- afg_fit
  temp$interventions$vaccine_strategy <- NULL
  expect_error(prepare_input_json_df(temp, ox_interventions = NULL))
})
