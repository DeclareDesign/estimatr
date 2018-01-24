context("Helper - significance testing")

test_that("Errors properly", {
  dat <- data.frame(
    mps = rep(1:4, each = 2),
    y = rnorm(8),
    z = c(0, 1)
  )

  expect_warning(
    lm_lin(y ~ z, ~ factor(mps), data = dat),
    "Some degrees of freedom have been estimated as negative or zero"
  )

  expect_error(
    lm_robust(y ~ z, data = dat, alpha = 10),
    "`alpha` must be numeric between 0 and 1"
  )
})
