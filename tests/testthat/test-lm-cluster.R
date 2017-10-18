context("lm cluster se")

test_that("lm cluster se", {

  dat <- data.frame(Y = rnorm(100),
                    Z = rbinom(100, 1, .5),
                    X = rnorm(100),
                    J = sample(1:10, 100, replace = T),
                    W = runif(100))


  ## Test functionality
  lm_robust(Y ~ Z, cluster_variable_name = J, data = dat)

  lm_robust(Y ~ Z + X, cluster_variable_name = J, data = dat)

  lm_robust(
    Y ~ Z + X,
    cluster_variable_name = J,
    coefficient_name = c("Z", "X"),
    data = dat
  )


  lm_robust(
    Y ~ Z + X,
    cluster_variable_name = J,
    coefficient_name = c("Z", "X"),
    se_type = 'stata',
    data = dat
  )

  expect_equivalent(
    unlist(
      lm_robust(
        Y ~ X + Z,
        cluster_variable_name = J,
        ci = FALSE,
        coefficient_name = c("X", "Z"),
        data = dat
      )[c("p", "ci_lower", "ci_upper")]
    ),
    rep(NA, 3)
  )

  ## Test equality
  lm_interact <-
    lm_robust(
      Y ~ Z*X,
      cluster_variable_name = J,
      coefficient_name = "Z:X",
      data = dat
    )

  lm_interact_stata <-
    lm_robust(
      Y ~ Z*X,
      cluster_variable_name = J,
      se_type = 'stata',
      coefficient_name = "Z:X",
      data = dat
    )

  lm_interact_simple <- lm(Y ~ Z*X, data = dat)

  bm_interact <-
    BMlmSE(
      lm_interact_simple,
      clustervar = as.factor(dat$J),
      IK = FALSE
    )

  bm_interact_interval <-
    lm_interact_simple$coefficients["Z:X"] +
    qt(0.975, df = bm_interact$dof["Z:X"]) * bm_interact$se["Z:X"] * c(-1, 1)

  bm_interact_stata_interval <-
    lm_interact_simple$coefficients["Z:X"] +
    qt(0.975, df = length(unique(dat$J))- 1) * bm_interact$se.Stata["Z:X"] * c(-1, 1)

  expect_equivalent(
    tidy(lm_interact)[c("se", "ci_lower", "ci_upper")],
    c(bm_interact$se["Z:X"], bm_interact_interval)
  )

  expect_equivalent(
    tidy(lm_interact_stata)[c("se", "ci_lower", "ci_upper")],
    c(bm_interact$se.Stata["Z:X"], bm_interact_stata_interval)
  )


  lm_full <-
    lm_robust(
      Y ~ Z + X,
      cluster_variable_name = J,
      data = dat
    )

  lm_full_simple <- lm(Y ~ Z + X, data = dat)

  bm_full <-
    BMlmSE(
      lm_full_simple,
      clustervar = as.factor(dat$J),
      IK = FALSE
    )

  bm_full_moe <- qt(0.975, df = bm_full$dof) * bm_full$se
  bm_full_lower <- lm_full_simple$coefficients - bm_full_moe
  bm_full_upper <- lm_full_simple$coefficients + bm_full_moe

  expect_equivalent(
    as.matrix( tidy(lm_full)[, c("se", "ci_lower", "ci_upper")] ),
    cbind(bm_full$se, bm_full_lower, bm_full_upper)
  )

  ## Test error handling
  expect_error(
    lm_robust(
      Y ~ Z,
      cluster_variable_name = J,
      weights = W,
      data = dat
    ),
    'weights'
  )

  expect_error(
    lm_robust(
      Y ~ Z,
      cluster_variable_name = J,
      se_type = 'HC2',
      data = dat
    ),
    'BM'
  )

  expect_error(
    lm_robust(
      Y ~ Z,
      se_type = 'BM',
      data = dat
    ),
    'BM'
  )
})

test_that("lm cluster se with missingness", {
  dat <- data.frame(Y = rnorm(100),
                    Z = rbinom(100, 1, .5),
                    X = rnorm(100),
                    J = sample(1:10, 100, replace = T),
                    W = runif(100))

  dat$X[23] <- NA
  dat$J[63] <- NA


  expect_warning(
    estimatr_cluster_out <- lm_robust(
      Y ~ Z + X,
      cluster_variable_name = J,
      data = dat
    ),
    'missingness in the cluster'
  )

  expect_identical(
    estimatr_cluster_out,
    lm_robust(
      Y ~ Z + X,
      cluster_variable_name = J,
      data = dat[-c(23, 63),]
    )
  )


})

test_that("lm works with quoted or unquoted vars and withor without factor clusters", {
  dat <- data.frame(Y = rnorm(100),
                    Z = rbinom(100, 1, .5),
                    X = rnorm(100),
                    J = sample(1:10, 100, replace = T),
                    W = runif(100))

  expect_identical(
    lm_robust(Y~Z, data = dat, weights = W),
    lm_robust(Y~Z, data = dat, weights = 'W')
  )

  # works with char
  dat$J <- as.character(dat$J)

  expect_identical(
    lm_robust(Y~Z, data = dat, cluster_variable_name = J),
    lm_robust(Y~Z, data = dat, cluster_variable_name = 'J')
  )


  # works with num
  dat$J <- as.numeric(dat$J)

  expect_identical(
    lm_robust(Y~Z, data = dat, cluster_variable_name = J),
    lm_robust(Y~Z, data = dat, cluster_variable_name = 'J')
  )

  dat$J_fac <- as.factor(dat$J)
  expect_identical(
    lm_robust(Y~Z, data = dat, cluster_variable_name = J_fac),
    lm_robust(Y~Z, data = dat, cluster_variable_name = J)
  )
})
