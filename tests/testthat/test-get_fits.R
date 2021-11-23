test_that("Errors", {
  #can't check other functions without adding orderly repo to package
  expect_error(get_fits(file.path("piiaiiisi"), date = "2021-01-01"))
})
