context("S3 - emmeans")



test_that("emmeans can work with lm_robust objects", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  lmr <- lm_robust(mpg ~ factor(cyl) * hp + wt, data = mtcars)

  rg <- emmeans::ref_grid(lmr)
  expect_equal(class(rg)[1], "emmGrid")

  grid <- rg@grid
  expect_equal(nrow(grid), 3)
  expect_equal(sum(grid$.wgt.), 32)
  expect_equal(predict(rg)[1], 17.424, tolerance = .01)
})

test_that("lm_robust multivariate model works with emmeans", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  lmr <- lm_robust(yield ~ Block + Variety, data = emmeans::MOats)
  emm <- emmeans(lmr, "rep.meas")
  expect_equal(summary(emm)$emmean[4], 123.4, tolerance = 0.1)
})

test_that("lm_robust model with rank deficiency works with emmeans", {
  skip_if_not_installed("emmeans")
  library(emmeans)
  lmr <- lm_robust(log(breaks) ~ wool * tension, data = warpbreaks, subset = -(19:30))
  pred <- predict(ref_grid(lmr))
  expect_true(is.na(pred[5]))
  expect_equal(length(pred), 6)
  expect_equal(sum(is.na(pred)), 1)
})

# Not testing emmeans package capabilities themselves. If we can construct the
# reference grid correctly, we are basically OK.
# Pretty much anything else that could fail would happen in the emmeans package,
# not in the support methods in this package.

