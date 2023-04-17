n_points <- 10

test_that("basic functionality", {
  data <- tibble(
    deaths = rpois(n_points, 1000),
    date_start = seq(as.Date("2020-01-01"), by = 2, length.out = n_points),
    date_end = seq(as.Date("2020-01-01"), by = 2, length.out = n_points + 1)[-1]
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
    deaths = rpois(n_points, 1000),
    date_start = seq(as.Date("2020-01-01"), by = 2, length.out = n_points),
    date_end = seq(as.Date("2020-01-01"), by = 2, length.out = n_points + 1)[-1]
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
    deaths = rpois(n_points, 1000),
    date_start = seq(as.Date("2020-01-01"), by = 2, length.out = n_points),
    date_end = seq(as.Date("2020-01-01"), by = 2, length.out = n_points + 1)[-1]
  )
  start_date <- min(data$date_start) - 20
  distribution <- map(
    seq_len(2), ~list(
      dur_R = rpois(1, 365),
      ICU_bed_capacity  = rpois(1, 7805)
    )
  )
  parameters <- list(
    country = "Turkey",
    protection_delay_time = n_points
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

test_that("nimue booster functionality difference", {
  data <- tibble(
    deaths = rpois(n_points, 1000),
    date_start = seq(as.Date("2020-01-01"), by = 2, length.out = n_points),
    date_end = seq(as.Date("2020-01-01"), by = 2, length.out = n_points + 1)[-1]
  )
  start_date <- min(data$date_start) - 20
  distribution <- map(
    seq_len(2), ~list(
      dur_R = rpois(1, 365),
      ICU_bed_capacity  = rpois(1, 7805)
    )
  )
  parameters <- list(
    country = "Turkey",
    protection_delay_time = n_points
  )
  squire_model <- nimue_booster_model(use_difference = TRUE)

  out <- rt_optimise(data, distribution, squire_model, parameters, start_date,
                     k = 7, n_particles = 10, dt = 0.5)

  expect_s3_class(out, "rt_optimised")
  expect_s3_class(out, "nimue_simulation")

  expect_s3_class(dp_plot(out), "ggplot")
  expect_s3_class(cdp_plot(out), "ggplot")
  expect_s3_class(ar_plot(out), "ggplot")

})
