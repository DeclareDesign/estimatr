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

  df <- data.frame(Y = rnorm(100), Z = rbinom(100, 1, .5), X = rnorm(100), W = runif(100))

  df$X[23] <- NA

  lm_robust(Y ~ Z + X, data = df)
  lm_robust(Y ~ Z + X, coefficient_name = "X", data = df)
  lm_robust(Y ~ Z + X, coefficient_name = c("Z", "X"), data = df)
  lm_robust(Y ~ Z + X, coefficient_name = c("(Intercept)", "Z", "X"), data = df)
  lm_robust(Y ~ Z*X, coefficient_name = "Z:X", data = df)

})
