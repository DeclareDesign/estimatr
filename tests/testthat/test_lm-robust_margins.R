context("lm robust margins")

skip_if_not_installed("margins")


test_that("lm robust can work with margins",{

  x <- lm(mpg ~ cyl * hp + wt, data = mtcars)
  x2 <- lm(mpg ~ cyl * hp + wt, data = mtcars)
  expect_identical(marginal_effects(x), marginal_effects(x2))

})

test_that("lm robust + weights can work with margins",{

  x <- lm(mpg ~ cyl * hp, data = mtcars, weights = wt)
  x2 <- lm(mpg ~ cyl * hp, data = mtcars, weights = wt)
  expect_identical(marginal_effects(x), marginal_effects(x2))

})

test_that("lm robust + cluster can work with margins",{

  skip("Doesn't work yet")
  x <- lm(mpg ~ cyl * hp + wt, data = mtcars)
  x2 <- lm_robust(mpg ~ cyl * hp + wt, data = mtcars, clusters = am)
  expect_identical(marginal_effects(x), marginal_effects(x2))

})

