context("lm robust se")

test_that("lm robust se",{

df <- data.frame(Y = rnorm(100), Z = rbinom(100, 1, .5), X = rnorm(100), W = runif(100))

lm_robust(Y ~ Z, data = df)

lm_robust(Y ~ Z, se_type = "none", data = df)

lm_robust(Y ~ Z + X, data = df)
lm_robust(Y ~ Z + X, coefficient_name = "X", data = df)
lm_robust(Y ~ Z + X, coefficient_name = c("Z", "X"), data = df)
lm_robust(Y ~ Z + X, coefficient_name = c("(Intercept)", "Z", "X"), data = df)
lm_robust(Y ~ Z*X, coefficient_name = "Z:X", data = df)

# we gotta figure out no quoting....
expect_error(lm_robust(Y ~ Z + X, coefficient_name = X, data = df))
expect_error(lm_robust(Y ~ Z + X, coefficient_name = c(Z, X), data = df))
expect_error(lm_robust(Y ~ Z + X, coefficient_name = c((Intercept), Z, X), data = df))
expect_error(lm_robust(Y ~ Z*X, coefficient_name = Z:X, data = df))


lm_robust(Y ~ Z, weights = W, data = df)

#matches.
#commarobust::commarobust(lm(Y ~ Z, weights = W, data = df))


})

context("lm robust se")

test_that("lm robust works with missingness",{

  df <- data.frame(Y = rnorm(100),
                   Z = rbinom(100, 1, .5),
                   X = rnorm(100),
                   W = runif(100))

  df$X[23] <- NA

  expect_identical(
    lm_robust(Y ~ Z + X, data = df),
    lm_robust(Y ~ Z + X, data = df[-23, ])
  )
  lm_robust(Y ~ Z + X, coefficient_name = "X", data = df)
  lm_robust(Y ~ Z + X, coefficient_name = c("Z", "X"), data = df)
  lm_robust(Y ~ Z + X, coefficient_name = c("(Intercept)", "Z", "X"), data = df)
  lm_robust(Y ~ Z*X, coefficient_name = "Z:X", data = df)

})


context("lm robust se")

test_that("lm robust works with weights",{

  df <- data.frame(Y = rnorm(100),
                   Z = rbinom(100, 1, .5),
                   X = rnorm(100),
                   W = runif(100))

  ## Make sure weighting works
  expect_error(
    estimatr_out <- lm_robust(Y ~ Z * X, weights = W, data = df),
    NA
  )

  # Compare to lm output
  lm_out <- lm(Y ~ Z * X, weights = W, data = df)
  lmo_hc2 <- cbind(lm_out$coefficients,
                   sqrt(diag(sandwich::vcovHC(lm_out, type = 'HC2'))))

  expect_equivalent(
    as.matrix(estimatr_out[, c('est', 'se')]),
    lmo_hc2
  )

  ## Make sure weighting works with missingness
  df$W[39] <- NA

  expect_warning(
    estimatr_miss_out <- lm_robust(Y ~ Z * X, weights = W, data = df),
    'missing'
  )

  expect_identical(
    estimatr_miss_out,
    lm_robust(Y ~ Z * X, weights = W, data = df[-39, ])
  )

  # Compare to lm output
  lm_miss_out <- lm(Y ~ Z * X, weights = W, data = df)
  lmo_miss_hc2 <- cbind(lm_miss_out$coefficients,
                        sqrt(diag(sandwich::vcovHC(lm_miss_out, type = 'HC2'))))

  expect_equivalent(
    as.matrix(estimatr_miss_out[, c('est', 'se')]),
    lmo_miss_hc2
  )

})
