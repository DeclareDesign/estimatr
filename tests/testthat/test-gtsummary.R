# This file is ignored by .Rbuildignore to keep from suggesting gt and gtsummary

context("S3 - gtsummary works")

test_that("gtsummary works with glance", {

  library(gt)
  library(gtsummary)

  model1 <- lm_robust(mpg ~ am, mtcars)
  model2 <- lm_robust(mpg ~ am, mtcars, clusters = cyl)
  model3 <- lm_lin(mpg ~ am, ~ cyl, mtcars)

  gto <- gtsummary(list(model1, model2, model3))

  expect_equal(ncol(gto), 3L)

  expect_equal(nrow(gto), 15L)

})
