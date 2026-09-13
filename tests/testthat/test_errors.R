library(estimatr)

# Every refusal, with the words that explain it.
#
# A refusal is only as useful as its message, and a message nothing asserts
# drifts: two of the branches below said something false before this file
# (review C8). Each test names the call, the reason it cannot be answered, and
# the fragment of the message that tells the user so. expect_error() with a
# pattern rather than a snapshot, because snapshots are skipped under R CMD
# check on CRAN and these branches should run everywhere.
#
# Refusals that belong to a specific design are asserted beside that design:
# degenerate fits in test_degenerate.R, blocked designs in
# test_blocked_variance.R, fixed-effects se_type combinations in
# test_fixed_effects.R.

set.seed(1)
n <- 40
err_data <- data.frame(
  y = rnorm(n),
  x = rnorm(n),
  z = rep(0:1, 20),
  z4 = rep(0:3, 10),
  cl = rep(1:10, each = 4),
  pair = rep(1:20, each = 2),
  w = runif(n, 0.5, 2)
)
# Alternates within each cluster of four.
err_data$z_within <- rep(0:1, 20)

test_that("alpha outside (0, 1) is refused", {
  expect_error(lm_robust(y ~ x, data = err_data, alpha = 1.5), "`alpha` must be numeric between 0 and 1")
  expect_error(lm_robust(y ~ x, data = err_data, alpha = 0), "`alpha` must be numeric between 0 and 1")
  expect_error(difference_in_means(y ~ z, data = err_data, alpha = 1),
               "`alpha` must be numeric between 0 and 1")
})

test_that("an se_type off the menu is refused, and the menu depends on clustering", {
  expect_error(lm_robust(y ~ x, data = err_data, se_type = "HC5"),
               "with no `clusters`.\nYou passed: HC5", fixed = TRUE)
  expect_error(lm_robust(y ~ x, data = err_data, se_type = "CR2"),
               "reserved for a case with clusters")
  expect_error(iv_robust(y ~ x | z, data = err_data, se_type = "CR0"),
               "reserved for a case with clusters")
  expect_error(lm_robust(y ~ x, data = err_data, clusters = cl, se_type = "HC2"),
               "when `clusters` are specified.\nYou passed: HC2", fixed = TRUE)
})

test_that("difference_in_means refuses designs it cannot estimate, and says why", {
  # Values the treatment never takes. This used to filter every row away and
  # then ask for both conditions "within each block".
  expect_error(difference_in_means(y ~ z, data = err_data, condition1 = 5, condition2 = 7),
               "`condition1` and `condition2` must be values found in the treatment")
  expect_error(difference_in_means(y ~ z4, data = err_data),
               "Treatment has > 2 values; must specify both `condition1` and `condition2`")
  expect_error(difference_in_means(y ~ z_within, data = err_data, clusters = cl),
               "All units within a cluster must have the same treatment condition")
  expect_error(difference_in_means(y ~ z, data = err_data, blocks = pair, weights = w),
               "Cannot use `weights` with matched pairs design")
})

test_that("lm_lin refuses a covariates formula it cannot use", {
  expect_error(lm_lin(y ~ z, covariates = y ~ x, data = err_data),
               "Must not specify a response variable in `covariates` formula")
  expect_error(lm_lin(y ~ z, covariates = ~ 1, data = err_data),
               "`covariates` must have a variable on the right-hand side, not 0 or 1")
})

test_that("horvitz_thompson refuses assignment probabilities it cannot match to units", {
  expect_error(horvitz_thompson(y ~ z, data = err_data),
               "Must supply `condition_prs`")
  expect_error(horvitz_thompson(y ~ z, data = err_data, condition_prs = "p"),
               "Unrecognised `condition_prs` format")
  # Named, but for conditions the treatment does not have. This used to be
  # reported as an unrecognised format.
  expect_error(horvitz_thompson(y ~ z, data = err_data, condition_prs = c(a = 0.5, b = 0.5)),
               "its names \\(a, b\\) do not include both conditions \\(0, 1\\)")
  expect_error(horvitz_thompson(y ~ z, data = err_data,
                                condition_prs = cbind(a = rep(0.5, n), b = rep(0.5, n))),
               "condition1/condition2 not found in condition_prs column names")
  expect_error(horvitz_thompson(y ~ z, data = err_data, condition_prs = c("0" = 0.5, "1" = 0.5),
                                condition1 = 0, condition2 = 2),
               "No units observed in one of the two conditions")
  skip_if_not_installed("randomizr")
  expect_error(horvitz_thompson(y ~ z, data = err_data,
                                condition_prs = randomizr::declare_ra(N = n + 2)),
               "declares 42 units but the data has 40 rows")
})

test_that("a formula that does not say what the estimator needs is refused", {
  expect_error(iv_robust(mpg ~ hp + cyl, data = mtcars),
               "Must specify a `formula` with both regressors and instruments")
  expect_error(lm_lin(y ~ z + x, covariates = ~ x, data = err_data),
               "must only have the treatment variable on the right-hand side")
  expect_error(lm_lin(y ~ z, err_data$x, data = err_data),
               "must be specified as a formula")
  expect_error(difference_in_means(y ~ z + x, data = err_data),
               "must have only one variable on the right-hand side")
})

test_that("blocks that cannot hold the design are refused", {
  set.seed(42)
  N <- 100
  d <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, 0.5),
                  bl = rep(1:10, each = 10), crossing_cl = rep(1:10, 10))
  expect_error(difference_in_means(Y ~ Z, blocks = bl, clusters = crossing_cl, data = d),
               "All `clusters` must be contained within `blocks`")
  d$one_unit_block <- c(1, rep(2:10, length.out = N - 1))
  expect_error(difference_in_means(Y ~ Z, blocks = one_unit_block, data = d),
               "All `blocks` must have multiple units")
})

test_that("weights that are infinite, or all zero, are refused", {
  # Both returned every coefficient as NA under a warning that called the
  # regressors collinear, in 1.0.6 as well. lm() refuses an infinite weight.
  infinite <- err_data
  infinite$w[3] <- Inf
  expect_error(lm_robust(y ~ x, data = infinite, weights = w), "`weights` must be finite")
  expect_error(iv_robust(y ~ x | z, data = infinite, weights = w), "`weights` must be finite")
  expect_error(difference_in_means(y ~ z, data = infinite, weights = w),
               "`weights` must be finite")

  zeroed <- err_data
  zeroed$w <- 0
  expect_error(lm_robust(y ~ x, data = zeroed, weights = w), "Every weight is zero")
  expect_error(lm_lin(y ~ z, covariates = ~ x, data = zeroed, weights = w),
               "Every weight is zero")

  # A missing weight is still a dropped row, not a refusal.
  missing <- err_data
  missing$w[3] <- NA
  expect_warning(lm_robust(y ~ x, data = missing, weights = w), "missingness in the weights")
})

test_that("horvitz_thompson refuses probabilities outside [0, 1]", {
  # A pair of probabilities like these came back as an ordinary estimate.
  expect_error(horvitz_thompson(y ~ z, data = err_data, condition_prs = c("0" = -0.2, "1" = 1.2)),
               "outside \\[0, 1\\]")
  expect_error(horvitz_thompson(y ~ z, data = err_data, condition_prs = rep(1.5, n)),
               "outside \\[0, 1\\]")
})

test_that("#297: lh_robust refuses a multivariate outcome and says why", {
  skip_if_not_installed("carData")
  expect_error(
    lh_robust(cbind(mpg, am) ~ cyl + gear, data = mtcars, linear_hypothesis = "cyl = 2"),
    "multiple outcomes"
  )
})
