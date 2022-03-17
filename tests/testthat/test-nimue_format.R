#generate data to use
data <- generate_draws(afg_fit, generate_parameters(afg_fit, draws = 3), draws = 3)
test_that("Summaries", {
  summaries <- c("deaths",
                 "infections",
                 "hospital_occupancy",
                 "ICU_occupancy",
                 "hospital_demand",
                 "ICU_demand",
                 "hospital_incidence",
                 "ICU_incidence",
                 "hospitalisations",
                 "vaccines",
                 "unvaccinated",
                 "vaccinated",
                 "priorvaccinated",
                 "long_covid")
  expect_true(
    all(
      summaries %in%
        nimue_format(
          data,
          var_select = summaries
        )$compartment
    )
  )
})
test_that("options", {
  expect_true(
    all(
      c("E1", "E2") %in%
        nimue_format(
          data,
          var_select = "E",
          reduce_age = TRUE,
          combine_compartments = FALSE
        )$compartment
    )
  )
  expect_true(
    "E" %in%
      nimue_format(
        data,
        var_select = "E",
        reduce_age = TRUE,
        combine_compartments = TRUE
      )$compartment
  )
  expect_true(
    "age_group" %in%
      names(nimue_format(
          data,
          var_select = "E",
          reduce_age = FALSE,
          combine_compartments = TRUE
        )
      )
  )
})
remove(data)
