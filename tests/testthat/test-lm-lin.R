context("lm lin")

test_that("Test LM Lin",{

  df <- data.frame(Y = rnorm(100),
                   Z = rbinom(100, 1, .5),
                   X1 = rnorm(100),
                   X2 = rbinom(100, 1, .5),
                   cluster = sample(1:10, size = 100, replace = T))
  lm_lin_out <- lm_lin(Y ~ Z, data = df, covariates = ~ X1 + X2)
  expect_error(
    lm_lin_out <- lm_lin(Y ~ Z, data = df, covariates = ~ X1 + X2),
    NA
  )

  expect_error(
    lm_lin(Y ~ Z + X1,
           covariates = ~ X2,
           data = df),
    'right-hand side'
  )

  df2 <- df
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

  # works with one covar
  expect_error(
    lm_lin(Y ~ Z,
           covariates = ~ X1,
           data = df),
    NA
  )

  df$X1_bar <- df$X1 - mean(df$X1)
  df$X2_bar <- df$X2 - mean(df$X2)

  lm_rob_out <- lm_robust(Y ~ Z + Z*X1_bar + Z*X2_bar, data = df)

  expect_identical(
    lm_lin_out,
    lm_rob_out
  )

  expect_equivalent(
    lm_lin_out$est,
    lm(Y ~ Z + Z*X1_bar + Z*X2_bar, data = df)$coefficients
  )


  ## Works with missing data in treatment
  df$Z[23] <- NA
  df$X1_bar <- df$X1 - mean(df$X1[-23])
  df$X2_bar <- df$X2 - mean(df$X2[-23])

  expect_identical(
    lm_lin(Y ~ Z,
           covariates = ~ X1 + X2,
           data = df),
    lm_robust(Y ~ Z + Z*X1_bar + Z*X2_bar,
              data = df)
  )

  ## Works with no intercept
  expect_identical(
    lm_lin(Y ~ Z + 0,
           covariates = ~ X1 + X2,
           data = df),
    lm_robust(Y ~ Z + Z*X1_bar + Z*X2_bar + 0,
              data = df)
  )

  ## Test cluster passes through
  expect_identical(
    lm_lin(Y ~ Z,
           covariates = ~ X1 + X2,
           data = df,
           cluster_variable_name = cluster),
    lm_robust(Y ~ Z + Z*X1_bar + Z*X2_bar,
              data = df,
              cluster_variable_name = cluster)
  )

  ## Test that it works with subset
  keep <- setdiff(which(df$Y > 0), 23)
  df$X1_bar <- df$X1 - mean(df$X1[keep])
  df$X2_bar <- df$X2 - mean(df$X2[keep])

  expect_identical(
    lm_lin(Y ~ Z,
           covariates = ~ X1 + X2,
           data = df,
           cluster_variable_name = cluster,
           subset = Y > 0),
    lm_robust(Y ~ Z + Z*X1_bar + Z*X2_bar,
              data = df,
              cluster_variable_name = cluster,
              subset = Y > 0)
  )

  # Works with factors
  df <- data.frame(Y = rnorm(100),
                   treat = as.factor(rbinom(100, 1, .5)),
                   X1 = rnorm(100),
                   X2 = as.factor(rbinom(100, 1, .5)),
                   cluster = sample(1:10, size = 100, replace = T))

  df$X1_bar <- df$X1 - mean(df$X1)
  df$X21_bar <- as.numeric(df$X2==1) - mean(df$X2 == 1)

  expect_equivalent(
    lm_lin(Y ~ treat,
           covariates = ~ X1 + X2,
           data = df,
           cluster_variable_name = cluster),
    lm_robust(Y ~ treat + treat*X1_bar + treat*X21_bar,
              data = df,
              cluster_variable_name = cluster)
  )

  ## Works with a factor with spaces in the name (often the case for clusters)
  df$X2 <- as.factor(sample(c("This is a level", "Get on my level"),
                            size = 100,
                            replace = T))
  ## for lm_robust
  df$X2_bar <-
    as.numeric(df$X2=="This is a level") - mean(df$X2 == "This is a level")

  ## Names will differ
  expect_equivalent(
    tidy(
      lm_lin(Y ~ treat,
             covariates = ~ X1 + X2,
             data = df,
             cluster_variable_name = cluster)
    )[, -1],
    tidy(
      lm_robust(Y ~ treat + treat*X1_bar + treat*X2_bar,
                data = df,
                cluster_variable_name = cluster)
    )[, -1]
  )

  ## Works with missingness on cluster
  df$cluster[40] <- NA
  df$X1_bar <- df$X1 - mean(df$X1[-40])
  df$X2_bar <-
    as.numeric(df$X2=="This is a level") - mean(df$X2[-40] == "This is a level")

  expect_warning(
    lin_out <- lm_lin(Y ~ treat,
                      covariates = ~ X1 + X2,
                      data = df,
                      cluster_variable_name = cluster),
    'missingness in the cluster'
  )

  expect_warning(
    rob_out <- lm_robust(Y ~ treat + treat*X1_bar + treat*X2_bar,
                         data = df,
                         cluster_variable_name = cluster),
    'missingness in the cluster'
  )

  # drop coefficient name because names will differ
  expect_equivalent(
    lin_out[-which(names(lin_out) == 'coefficient_name')],
    rob_out[-which(names(rob_out) == 'coefficient_name')]
  )
})
