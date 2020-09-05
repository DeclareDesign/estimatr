# This file is ignored by .Rbuildignore to keep from suggesting gt and modelsummary

context("S3 - modelsummary works")

test_that("modelsummary works with glance", {

  library(modelsummary)

  model1 <- lm_robust(mpg ~ am, mtcars)
  model2 <- lm_robust(mpg ~ am, mtcars, clusters = cyl)
  model3 <- lm_lin(mpg ~ am, ~ cyl, mtcars)

  mso <- modelsummary::modelsummary(list(model1, model2, model3), output = "data.frame")

  expect_equal(colnames(mso), c("group", "term", "statistic",
                                "Model 1", "Model 2", "Model 3"))

  expect_equal(nrow(mso), 12L)

  expect_equal(ncol(mso), 6L)

  # iv_robust
  model1 <- iv_robust(mpg ~ am | gear, mtcars)
  model2 <- iv_robust(mpg ~ am | gear, mtcars, clusters = cyl, diagnostics = TRUE)

  mso <- modelsummary::modelsummary(list(model1, model2),
                               gof_omit = c("N|[sS]tatistic|p.value|p{1}"),
                               output = "data.frame")

  expect_equal(nrow(mso), 6)

  expect_equal(ncol(mso), 5)

  # difference_in_means
  model1 <- difference_in_means(mpg ~ am, mtcars)
  model2 <- difference_in_means(mpg ~ am, mtcars, blocks = vs)
  mso <- modelsummary:::extract_models(list(model1, model2))


  # horvitz_thompson
  model1 <- horvitz_thompson(mpg ~ am, mtcars)
  model2 <- horvitz_thompson(mpg ~ am, mtcars, blocks = vs)

  mso <- modelsummary::modelsummary(list(model1, model2), output = "data.frame")

  expect_equal(nrow(mso), 6)

  expect_equal(ncol(mso), 5)

})
