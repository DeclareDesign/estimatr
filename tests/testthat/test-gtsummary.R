# This file is ignored by .Rbuildignore to keep from suggesting gt and gtsummary

context("S3 - gtsummary works")

test_that("gtsummary works with glance", {

  library(gt)
  library(gtsummary)

  model1 <- lm_robust(mpg ~ am, mtcars)
  model2 <- lm_robust(mpg ~ am, mtcars, clusters = cyl)
  model3 <- lm_lin(mpg ~ am, ~ cyl, mtcars)

  gto <- gtsummary(list(model1, model2, model3))

  expect_equal(colnames(gto), c("Model 1", "Model 2", "Model 3"))

  expect_equal(nrow(gto), 15L)

  # iv_robust
  model1 <- iv_robust(mpg ~ am | gear, mtcars)
  model2 <- iv_robust(mpg ~ am | gear, mtcars, clusters = cyl, diagnostics = TRUE)

  gtsummary(list(model1, model2), gof_omit = c("N|[sS]tatistic|p.value|p{1}"))

  # difference_in_means
  model1 <- difference_in_means(mpg ~ am, mtcars)
  model2 <- difference_in_means(mpg ~ am, mtcars, blocks = vs)

  gtsummary(list(model1, model2))

  # horvitz_thompson
  model1 <- horvitz_thompson(mpg ~ am, mtcars)
  model2 <- horvitz_thompson(mpg ~ am, mtcars, blocks = vs)

  gtsummary(list(model1, model2))
})
