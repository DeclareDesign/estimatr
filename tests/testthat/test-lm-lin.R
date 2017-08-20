context("lm lin")

test_that("Test LM Lin",{

  df <- data.frame(Y = rnorm(100),
                   Z = rbinom(100, 1, .5),
                   X1 = rnorm(100),
                   X2 = rbinom(100, 1, .5),
                   cluster = sample(1:10, size = 100, replace = T))

  expect_error(
    lm_lin_out <- lm_lin(Y ~ Z, data = df, covariates = ~ X1 + X2),
    NA
  )

  expect_error(
    lm_lin(Y ~ Z + X1,
           covariates = ~ X2
           data = df),
    'right-hand side'
  )

  df2$Z <- rnorm(100)
  expect_error(
    lm_lin(Y ~ Z,
           covariates = ~ X1,
           data = df2),
    'binary'
  )

  expect_error(
    lm_lin(Y ~ Z,
           covariates = Y ~ X1,
           data = df),
    'right-hand sided equation'
  )

  df$X1_bar <- df$X1 - mean(df$X1)
  df$X2_bar <- df$X2 - mean(df$X2)

  lm_rob_out <- lm_robust(Y ~ Z + Z*X1_bar + Z*X2_bar, data = df)

  expect_identical(
    lm_lin_out,
    lm_rob_out
  )

  expect_equivalent(
    lm_lin_out[, c('est')],
    lm(Y ~ Z + Z*X1_bar + Z*X2_bar, data = df)$coefficients
  )

  expect_identical(
    lm_lin(Y ~ Z,
           covariates = ~ X1 + X2,
           data = df,
           cluster_variable_name = cluster),
    lm_robust(Y ~ Z + Z*X1_bar + Z*X2_bar,
              data = df,
              cluster_variable_name = cluster)
  )

})
