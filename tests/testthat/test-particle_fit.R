test_that("basic functionality", {
  data <- tibble(
    deaths = rpois(100, 100),
    date_start = seq(as.Date("2020-01-01"), by = 2, length.out = 100),
    date_end = seq(as.Date("2020-01-01"), by = 2, length.out = 101)[-1]
  )
  start_date <- as.Date("2020-01-01") - 30
  distribution <- map(
    seq_len(3), ~list(
      dur_R = rpois(1, 365),
      ICU_bed_capacity  = rpois(1, 7805)
    )
  )
  parameters <- list(
    country = "Turkey"
  )
  squire_model <- squire:::deterministic_model()

  out <- particle_fit(data, distribution, squire_model, parameters, start_date)

  expect_s3_class(out, "particle_fit")

  expect_s3_class(dp_plot(out), "ggplot")
  expect_s3_class(cdp_plot(out), "ggplot")
  expect_s3_class(ar_plot(out), "ggplot")

})
