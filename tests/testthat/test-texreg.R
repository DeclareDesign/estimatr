# This file is ignored by .Rbuildignore to keep from suggesting broom

context("S3 - tidy and broom compatability")

test_that("estimatr::tidy works loaded before or after after broom", {

  model2 <- lm_robust(extra~group, sleep, clusters = "ID")

  expect_error(
    extract(model2),
    "could not find function"
  )

  library(texreg)

  expect_is(
    extract(model2),
    "texreg"
  )

})
