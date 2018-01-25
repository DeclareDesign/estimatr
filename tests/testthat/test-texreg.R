# This file is ignored by .Rbuildignore to keep from suggesting texreg

context("S3 - texreg builds")

test_that("texreg extension works", {

  model2 <- lm_robust(extra~group, sleep, clusters = "ID")

  library(texreg)

  capture.output(treg <- extract(model2, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE))
  expect_is(
    treg,
    "texreg"
  )

})
