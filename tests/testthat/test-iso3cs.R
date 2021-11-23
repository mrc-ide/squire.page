test_that("Simple list functions", {
  expect_type(get_covax_iso3c(), "character")
  expect_type(get_gavi_iso3c(), "character")
  expect_type(get_global_fund_iso3c(), "character")
})
test_that("Classifiers", {
  expect_warning(
    get_income_group(c("GBR", "KKLALL"))
  )
  expect_warning(
    get_WHO_region(c("GBR", "KKLALL"))
  )
  expect_type(
    get_income_group("GBR"), "integer"
  )
  expect_type(
    get_WHO_region("GBR"), "character"
  )
})
