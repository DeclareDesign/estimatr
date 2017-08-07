context("lm cluster se")

test_that("lm cluster se", {

  df <- data.frame(Y = rnorm(100),
                   Z = rbinom(100, 1, .5),
                   X = rnorm(100),
                   J = sample(1:10, 100, replace = T),
                   W = runif(100))


  ## Test functionality
  lm_robust_se(Y ~ Z, cluster_variable_name = J, data = df)

  lm_robust_se(Y ~ Z + X, cluster_variable_name = J, data = df)

  lm_robust_se(
    Y ~ Z + X,
    cluster_variable_name = J,
    coefficient_name = c("Z", "X"),
    data = df
  )

  ## Test equality
  lm_interact <-
    lm_robust_se(
      Y ~ Z*X,
      cluster_variable_name = J,
      coefficient_name = "Z:X",
      data = df
    )

  lm_interact_simple <- lm(Y ~ Z*X, data = df)

  bm_interact <-
    BMlmSE(
      lm_interact_simple,
      clustervar = as.factor(df$J),
      IK = FALSE
    )

  bm_interact_interval <-
    lm_interact_simple$coefficients["Z:X"] +
    qt(0.975, df = bm_interact$dof["Z:X"]) * bm_interact$se["Z:X"] * c(-1, 1)

  expect_equivalent(
    lm_interact[c("se", "ci_lower", "ci_upper")],
    c(bm_interact$se["Z:X"], bm_interact_interval)
  )

  lm_full <-
    lm_robust_se(
      Y ~ Z + X,
      cluster_variable_name = J,
      coefficient_name = c("(Intercept)", "Z", "X"),
      data = df
    )

  lm_full_simple <- lm(Y ~ Z + X, data = df)

  bm_full <-
    BMlmSE(
      lm_full_simple,
      clustervar = as.factor(df$J),
      IK = FALSE
    )

  bm_full_moe <- qt(0.975, df = bm_full$dof) * bm_full$se
  bm_full_lower <- lm_full_simple$coefficients - bm_full_moe
  bm_full_upper <- lm_full_simple$coefficients + bm_full_moe

  expect_equivalent(
    as.matrix( lm_full[, c("se", "ci_lower", "ci_upper")] ),
    cbind(bm_full$se, bm_full_lower, bm_full_upper)
  )

  ## Test error handling
  expect_error(
    lm_robust_se(
      Y ~ Z,
      cluster_variable_name = J,
      weights = W,
      data = df
    ),
    'weights'
  )

  expect_error(
    lm_robust_se(
      Y ~ Z,
      cluster_variable_name = J,
      se_type = 'HC2',
      data = df
    ),
    'BM'
  )

  expect_error(
    lm_robust_se(
      Y ~ Z,
      se_type = 'BM',
      data = df
    ),
    'BM'
  )
})

