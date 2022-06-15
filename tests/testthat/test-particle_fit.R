test_that("basic functionality", {
  data <- tibble(
    deaths = rpois(50, 1000),
    date_start = seq(as.Date("2020-01-01"), by = 2, length.out = 50),
    date_end = seq(as.Date("2020-01-01"), by = 2, length.out = 51)[-1]
  )
  start_date <- min(data$date_start) - 20
  distribution <- map(
    seq_len(2), ~list(
      dur_R = rpois(1, 365),
      ICU_bed_capacity  = rpois(1, 7805)
    )
  )
  parameters <- list(
    country = "Turkey"
  )
  squire_model <- squire:::deterministic_model()

  out <- rt_optimise(data, distribution, squire_model, parameters, start_date,
                      k = 7, n_particles = 10)

  expect_s3_class(out, "rt_optimised")

  expect_s3_class(dp_plot(out), "ggplot")
  expect_s3_class(cdp_plot(out), "ggplot")
  expect_s3_class(ar_plot(out), "ggplot")
})

test_that("nimue functionality", {
  data <- tibble(
    deaths = rpois(50, 1000),
    date_start = seq(as.Date("2020-01-01"), by = 2, length.out = 50),
    date_end = seq(as.Date("2020-01-01"), by = 2, length.out = 51)[-1]
  )
  start_date <- min(data$date_start) - 20
  distribution <- map(
    seq_len(2), ~list(
      dur_R = rpois(1, 365),
      ICU_bed_capacity  = rpois(1, 7805)
    )
  )
  parameters <- list(
    country = "Turkey"
  )
  squire_model <- nimue::nimue_deterministic_model()

  out <- rt_optimise(data, distribution, squire_model, parameters, start_date,
                      k = 7, n_particles = 10)

  expect_s3_class(out, "rt_optimised")
  expect_s3_class(out, "nimue_simulation")

  expect_s3_class(dp_plot(out), "ggplot")
  expect_s3_class(cdp_plot(out), "ggplot")
  expect_s3_class(ar_plot(out), "ggplot")

})

test_that("nimue booster functionality", {
  data <- tibble(
    deaths = rpois(50, 1000),
    date_start = seq(as.Date("2020-01-01"), by = 2, length.out = 50),
    date_end = seq(as.Date("2020-01-01"), by = 2, length.out = 51)[-1]
  )
  start_date <- min(data$date_start) - 20
  distribution <- map(
    seq_len(2), ~list(
      dur_R = rpois(1, 365),
      ICU_bed_capacity  = rpois(1, 7805)
    )
  )
  parameters <- list(
    country = "Turkey"
  )
  squire_model <- nimue_booster_model()

  out <- rt_optimise(data, distribution, squire_model, parameters, start_date,
                      k = 7, n_particles = 10)

  expect_s3_class(out, "rt_optimised")
  expect_s3_class(out, "nimue_simulation")

  expect_s3_class(dp_plot(out), "ggplot")
  expect_s3_class(cdp_plot(out), "ggplot")
  expect_s3_class(ar_plot(out), "ggplot")

})
