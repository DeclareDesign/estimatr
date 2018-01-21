context("lm robust se")

test_that("lm robust se",{

  N <- 100
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5), X = rnorm(N), W = runif(N))

  lm_robust(Y ~ Z, data = dat)

  lm_robust(Y ~ Z, se_type = "none", data = dat)

  lm_robust(Y ~ Z + X, data = dat)
  lm_robust(Y ~ Z + X, coefficient_name = "X", data = dat)
  lm_robust(Y ~ Z + X, coefficient_name = c("Z", "X"), data = dat)
  lm_robust(Y ~ Z + X, coefficient_name = c("(Intercept)", "Z", "X"), data = dat)
  lm_robust(Y ~ Z*X, coefficient_name = "Z:X", data = dat)

  lm_robust(Y ~ Z + X, data = dat, subset = W > 0.5)
  # Works with subset
  expect_identical(
    lm_robust(Y ~ Z + X, data = dat, subset = W > 0.5),
    lm_robust(Y ~ Z + X, data = dat[dat$W > 0.5, ])
  )

  # we gotta figure out no quoting....
  expect_error(lm_robust(Y ~ Z + X, coefficient_name = c(Z, X), data = dat))
  expect_error(lm_robust(Y ~ Z + X, coefficient_name = c((Intercept), Z, X), data = dat))
  expect_error(lm_robust(Y ~ Z*X, coefficient_name = Z:X, data = dat))


  lm_robust(Y ~ Z, weights = W, data = dat)

  #matches.
  #commarobust::commarobust(lm(Y ~ Z, weights = W, data = dat))

  # To easily do with and without weights
  test_lm_robust_variance <- function(w) {
    # Test other estimators
    lm_hc0 <- lm_robust(Y ~ Z + X, data = dat, weights = w, se_type = "HC0")
    lm_hc1 <- lm_robust(Y ~ Z + X, data = dat, weights = w, se_type = "HC1")
    lm_hc2 <- lm_robust(Y ~ Z + X, data = dat, weights = w, se_type = "HC2")
    lm_hc3 <- lm_robust(Y ~ Z + X, data = dat, weights = w, se_type = "HC3")
    lm_stata <- lm_robust(Y ~ Z + X, data = dat, weights = w, se_type = "stata")

    # Stata is the same as HC1
    expect_identical(
      lm_hc1,
      lm_stata
    )

    expect_false(all(lm_hc0$se == lm_hc1$se))
    expect_false(all(lm_hc0$se == lm_hc2$se))
    expect_false(all(lm_hc0$se == lm_hc3$se))
    expect_false(all(lm_hc1$se == lm_hc2$se))
    expect_false(all(lm_hc1$se == lm_hc3$se))
    expect_false(all(lm_hc2$se == lm_hc3$se))

    expect_equivalent(
      lm_hc0$df,
      lm_hc1$df,
      lm_hc2$df,
      lm_hc3$df,
      lm_stata$df
    )

    expect_equivalent(
      lm_hc0$se^2,
      lm_hc1$se^2 * ((N - length(lm_hc1$est)) / N)
    )
  }

  # No weights first
  test_lm_robust_variance(NULL)
  test_lm_robust_variance(dat$W)


})

test_that("lm robust works with missingness",{

  dat <- data.frame(Y = rnorm(100),
                   Z = rbinom(100, 1, .5),
                   X = rnorm(100),
                   W = runif(100))

  dat$X[23] <- NA

  expect_identical(
    lm_robust(Y ~ Z + X, data = dat),
    lm_robust(Y ~ Z + X, data = dat[-23, ])
  )
  lm_robust(Y ~ Z + X, coefficient_name = "X", data = dat)
  lm_robust(Y ~ Z + X, coefficient_name = c("Z", "X"), data = dat)
  lm_robust(Y ~ Z + X, coefficient_name = c("(Intercept)", "Z", "X"), data = dat)
  lm_robust(Y ~ Z*X, coefficient_name = "Z:X", data = dat)

  ## Outcome missingness
  dat$Y[35] <- NA

  estimatr_missout_out <- lm_robust(Y ~ Z + X, data = dat)

  lm_missout_out <- lm(Y ~ Z + X, data = dat)
  lm_missout_hc2 <- cbind(
    lm_missout_out$coefficients,
    sqrt(diag(sandwich::vcovHC(lm_missout_out, type = "HC2")))
  )

  expect_equivalent(
    as.matrix(tidy(estimatr_missout_out)[, c("est", "se")]),
    lm_missout_hc2
  )



})

test_that("lm robust works with weights",{

  N <- 100
  dat <- data.frame(Y = rnorm(N),
                   Z = rbinom(N, 1, .5),
                   X = rnorm(N),
                   W = runif(N))

  ## Make sure weighting works
  expect_error(
    estimatr_out <- lm_robust(Y ~ Z * X, weights = W, data = dat),
    NA
  )

  # Compare to lm output
  lm_out <- lm(Y ~ Z * X, weights = W, data = dat)
  lmo_hc2 <- cbind(lm_out$coefficients,
                   sqrt(diag(sandwich::vcovHC(lm_out, type = 'HC2'))))

  expect_equivalent(
    as.matrix(tidy(estimatr_out)[, c('est', 'se')]),
    lmo_hc2
  )

  ## Make sure weighting works with missingness
  dat$W[39] <- NA

  expect_warning(
    estimatr_miss_out <- lm_robust(Y ~ Z * X, weights = W, data = dat),
    'missing'
  )

  expect_identical(
    estimatr_miss_out,
    lm_robust(Y ~ Z * X, weights = W, data = dat[-39, ])
  )

  # Compare to lm output
  lm_miss_out <- lm(Y ~ Z * X, weights = W, data = dat)
  lmo_miss_hc2 <- cbind(lm_miss_out$coefficients,
                        sqrt(diag(sandwich::vcovHC(lm_miss_out, type = 'HC2'))))

  expect_equivalent(
    as.matrix(tidy(estimatr_miss_out)[, c('est', 'se')]),
    lmo_miss_hc2
  )

  expect_error(
    lm_robust(Y ~ Z, data = dat, weights = c(0, runif(N - 1))),
    "weights must all be positive and non-zero."
  )

})

test_that("lm robust works with large data", {
  N <- 75000
  dat <- data.frame(Y = rbinom(N, 1, .5),
                   X1 = rnorm(N),
                   X2 = rnorm(N),
                   X3 = rnorm(N))
  expect_error(
    lm_robust(Y ~ X1 + X2 + X3, data = dat, se_type = 'none'),
    NA
  )

})

test_that("lm robust works with rank-deficient X", {
  N <- 100
  dat <- data.frame(Y = rbinom(N, 1, .5),
                   X1 = rnorm(N),
                   X2 = rnorm(N),
                   X3 = rnorm(N))

  dat$Z1 <- dat$X1
  sum_lm <- summary(lm(Y ~ X1 + X2 + Z1 + X3, data = dat))
  ## manually build vector of coefficients, can't extract from summary.lm
  out_sumlm <- matrix(NA, nrow = length(sum_lm$aliased), ncol = 2)
  j <- 1
  for(i in seq_along(sum_lm$aliased)) {
    if (!sum_lm$aliased[i]) {
      out_sumlm[i, ] <- sum_lm$coefficients[j, 1:2]
      j <- j + 1
    }
  }

  ## order sometimes is different! Not stable order!
  # expect_equivalent(
  #   as.matrix(tidy(lm_robust(Y ~ X1 + X2 + Z1 + X3, data = dat, se_type = 'classical'))[, c('est', 'se')]),
  #   out_sumlm
  # )

  dat$Z1 <- dat$X1 + 5

  ## Not the same as LM! Different QR decompositions when dependency isn't just equivalency
  expect_equivalent(
    as.matrix(tidy(lm_robust(Y ~ X1 + X2 + Z1 + X3, data = dat, se_type = 'classical'))[, c('est', 'se')]),
    as.matrix(summary(RcppEigen::fastLm(Y ~ X1 + X2 + Z1 + X3, data = dat))$coefficients[, 1:2])
  )

})
