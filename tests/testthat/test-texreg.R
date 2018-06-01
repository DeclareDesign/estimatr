# This file is ignored by .Rbuildignore to keep from suggesting texreg

context("S3 - texreg builds")

test_that("texreg extension works", {

  model2 <- lm_robust(extra~group, sleep, clusters = ID)

  capture.output(treg <- extract(model2, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE))
  expect_is(
    treg,
    "texreg"
  )

  # Defaults to having CIs
  expect_true(grepl("-0.53;", texreg::texreg(model2)))

  # Remove to get SEs
  expect_true(grepl("\\(0.57\\)", texreg::texreg(model2, include.ci = FALSE)))

})
