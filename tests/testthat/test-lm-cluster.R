context("Estimator - lm_robust, clustered")

test_that("lm cluster se", {
  N <- 100
  dat <- data.frame(
    Y = rnorm(N),
    Z = rbinom(N, 1, .5),
    X = rnorm(N),
    J = sample(1:10, N, replace = T),
    W = runif(N)
  )


  ## Test functionality
  lm_robust(Y ~ Z, clusters = J, data = dat)

  lm_robust(Y ~ Z + X, clusters = J, data = dat)

  lm_robust(
    Y ~ Z + X,
    clusters = J,
    data = dat
  )


  lm_robust(
    Y ~ Z + X,
    clusters = J,
    se_type = "stata",
    data = dat,
    ci = T
  )

  expect_equivalent(
    as.matrix(
      tidy(
        lm_robust(
          Y ~ X + Z,
          clusters = J,
          ci = FALSE,
          data = dat
        )
      )[, c("p.value", "conf.low", "conf.high")]
    ),
    matrix(NA, nrow = 3, ncol = 3)
  )

  ## Test equality
  lm_interact <-
    lm_robust(
      Y ~ Z * X,
      clusters = J,
      data = dat
    )

  lm_interact_stata <-
    lm_robust(
      Y ~ Z * X,
      clusters = J,
      se_type = "stata",
      data = dat
    )

  lm_interact_simple <- lm(Y ~ Z * X, data = dat)

  bm_interact <-
    BMlmSE(
      lm_interact_simple,
      clustervar = as.factor(dat$J),
      IK = FALSE
    )

  bm_interact

  bm_interact_interval <-
    coef(lm_interact_simple)["Z:X"] +
    qt(0.975, df = bm_interact$dof["Z:X"]) * bm_interact$se["Z:X"] * c(-1, 1)

  bm_interact_stata_interval <-
    coef(lm_interact_simple)["Z:X"] +
    qt(0.975, df = length(unique(dat$J)) - 1) * bm_interact$se.Stata["Z:X"] * c(-1, 1)

  expect_equivalent(
    tidy(lm_interact)[4, c("std.error", "conf.low", "conf.high")],
    c(bm_interact$se["Z:X"], bm_interact_interval)
  )

  expect_equivalent(
    tidy(lm_interact_stata)[4, c("std.error", "conf.low", "conf.high")],
    c(bm_interact$se.Stata["Z:X"], bm_interact_stata_interval)
  )


  lm_full <-
    lm_robust(
      Y ~ Z + X,
      clusters = J,
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
  bm_full_lower <- coef(lm_full_simple) - bm_full_moe
  bm_full_upper <- coef(lm_full_simple) + bm_full_moe

  expect_equivalent(
    as.matrix(tidy(lm_full)[, c("std.error", "conf.low", "conf.high")]),
    cbind(bm_full$se, bm_full_lower, bm_full_upper)
  )

  ## Works with rank deficient case
  dat$X2 <- dat$X
  lmr_rd <- lm_robust(Y ~ X + Z + X2, data = dat, clusters = J, se_type = "stata")
  lmr_full <- lm_robust(Y ~ X + Z, data = dat, clusters = J, se_type = "stata")
  expect_identical(
    tidy(lmr_rd)[1:3, ],
    tidy(lmr_full)
  )

  lmr_rd_cr2 <- lm_robust(Y ~ X + Z + X2, data = dat, clusters = J, se_type = "CR2")
  lmr_full_cr2 <- lm_robust(Y ~ X + Z, data = dat, clusters = J, se_type = "CR2")
  expect_identical(
    tidy(lmr_rd_cr2)[1:3, ],
    tidy(lmr_full_cr2)
  )

  ## Test error handling
  expect_error(
    lm_robust(
      Y ~ Z,
      clusters = J,
      se_type = "HC2",
      data = dat
    ),
    "CR2"
  )

  expect_error(
    lm_robust(
      Y ~ Z,
      se_type = "CR2",
      data = dat
    ),
    "CR2"
  )

  # To easily do with and without weights
  test_lm_cluster_variance <- function(w) {
    # Test other estimators
    lm_cr0 <- lm_robust(Y ~ Z + X, data = dat, weights = w, clusters = J, se_type = "CR0")
    lm_stata <- lm_robust(Y ~ Z + X, data = dat, weights = w, clusters = J, se_type = "stata")
    lm_cr2 <- lm_robust(Y ~ Z + X, data = dat, weights = w, clusters = J, se_type = "CR2")

    # Stata is the same as CR0 but with finite sample
    expect_equivalent(
      lm_cr0$std.error ^ 2,
      lm_stata$std.error ^ 2 * (N - length(coef(lm_stata))) * (length(unique(dat$J)) - 1) / ((N - 1) * length(unique(dat$J)))
    )

    expect_false(all(lm_cr0$std.error == lm_stata$std.error))
    expect_false(all(lm_cr0$std.error == lm_cr2$std.error))
    expect_false(all(lm_stata$std.error == lm_cr2$std.error))
    expect_false(all(lm_stata$df == lm_cr2$df))

    expect_equivalent(
      lm_cr0$df,
      lm_stata$df
    )
  }

  # No weights first
  test_lm_cluster_variance(NULL)
  test_lm_cluster_variance(dat$W)

})

test_that("Clustered weighted SEs are correct", {
  lm_cr2 <- lm_robust(mpg ~ hp, data = mtcars, weights = wt, clusters = cyl, se_type = "CR2")
  lm_stata <- lm_robust(mpg ~ hp, data = mtcars, weights = wt, clusters = cyl, se_type = "stata")
  lm_cr0 <- lm_robust(mpg ~ hp, data = mtcars, weights = wt, clusters = cyl, se_type = "CR0")

  lm_o <- lm(mpg ~ hp, data = mtcars, weights = wt)

  expect_equivalent(
    vcov(lm_cr2),
    as.matrix(clubSandwich::vcovCR(lm_o, cluster = mtcars$cyl, type = "CR2"))
  )
  expect_equivalent(
    vcov(lm_cr0),
    as.matrix(clubSandwich::vcovCR(lm_o, cluster = mtcars$cyl, type = "CR0"))
  )
  expect_equivalent(
    vcov(lm_stata),
    as.matrix(clubSandwich::vcovCR(lm_o, cluster = mtcars$cyl, type = "CR1S"))
  )
})

test_that("lm cluster se with missingness", {
  dat <- data.frame(
    Y = rnorm(100),
    Z = rbinom(100, 1, .5),
    X = rnorm(100),
    J = sample(1:10, 100, replace = T),
    W = runif(100)
  )

  dat$X[23] <- NA
  dat$J[63] <- NA

  expect_warning(
    estimatr_cluster_out <- lm_robust(
      Y ~ Z + X,
      clusters = J,
      data = dat
    ),
    "missingness in the cluster"
  )

  estimatr_cluster_sub <- lm_robust(
    Y ~ Z + X,
    clusters = J,
    data = dat[-c(23, 63), ]
  )

  estimatr_cluster_out[["call"]] <- NULL
  estimatr_cluster_sub[["call"]] <- NULL
  expect_equal(
    estimatr_cluster_out,
    estimatr_cluster_sub
  )
})

test_that("lm works with quoted or unquoted vars and withor without factor clusters", {
  dat <- data.frame(
    Y = rnorm(100),
    Z = rbinom(100, 1, .5),
    X = rnorm(100),
    J = sample(1:10, 100, replace = T),
    W = runif(100)
  )

  lmr <- lm_robust(Y~Z, data = dat, weights = W)
  lmrq <- lm_robust(Y~Z, data = dat, weights = W)
  lmr[["call"]] <- NULL
  lmrq[["call"]] <- NULL
  expect_equal(
    lmr,
    lmrq
  )

  # works with char
  dat$J <- as.character(dat$J)

  lmrc <- lm_robust(Y~Z, data = dat, clusters = J)
  lmrcq <- lm_robust(Y~Z, data = dat, clusters = J)
  lmrc[["call"]] <- NULL
  lmrcq[["call"]] <- NULL
  expect_equal(
    lmrc,
    lmrcq
  )

  # works with num
  dat$J_num <- as.numeric(dat$J)

  lmrc_qnum <- lm_robust(Y~Z, data = dat, clusters = J_num)
  lmrc_qnum[["call"]] <- NULL
  expect_equal(
    lmrc,
    lmrc_qnum
  )


  # works with factor
  dat$J_fac <- as.factor(dat$J)
  expect_equivalent(
    rmcall(lm_robust(Y~Z, data = dat, clusters = J_fac)),
    rmcall(lm_robust(Y~Z, data = dat, clusters = J))
  )

  # works with being cast in the call
  lm_robust(Y~Z, data = dat, clusters = as.factor(J))
})

test_that("Clustered SEs work with clusters of size 1", {
  dat <- data.frame(
    Y = rnorm(100),
    X = rnorm(100),
    J = 1:100
  )

  lm_cr2 <- lm_robust(Y ~ X, data = dat, clusters = J)
  lm_stata <- lm_robust(Y ~ X, data = dat, clusters = J, se_type = "stata")
  lmo <- lm(Y ~ X, data = dat)

  bmo <-
    BMlmSE(
      lmo,
      clustervar = as.factor(dat$J),
      IK = FALSE
    )

  expect_equivalent(
    as.matrix(tidy(lm_cr2)[, c("estimate", "std.error", "df")]),
    cbind(coef(lmo), bmo$se, bmo$dof)
  )

  expect_equivalent(
    as.matrix(tidy(lm_stata)[, c("estimate", "std.error")]),
    cbind(coef(lmo), bmo$se.Stata)
  )
})

test_that("multiple outcomes", {

  lmo <- lm(cbind(mpg, hp) ~ wt, data = mtcars)
  lmro <- lm_robust(cbind(mpg, hp) ~ wt, data = mtcars, clusters = cyl)

  expect_equivalent(
    as.matrix(clubSandwich::vcovCR(lmo, cluster = mtcars$cyl, type = "CR2")),
    vcov(lmro)
  )

  expect_equivalent(
    as.matrix(clubSandwich::vcovCR(lmo, cluster = mtcars$cyl, type = "CR0")),
    vcov(lm_robust(cbind(mpg, hp) ~ wt, data = mtcars, clusters = cyl, se_type = "CR0"))
  )

  lmmpg <- lm_robust(mpg ~ wt, data = mtcars, clusters = cyl)
  lmhp <- lm_robust(hp ~ wt, data = mtcars, clusters = cyl)

  expect_equivalent(
    tidy(lmro)$df[1:2],
    lmmpg$df
  )

  expect_equivalent(
    tidy(lmro)$df[3:4],
    lmhp$df
  )

  expect_equivalent(lmro$r.squared[1], lmmpg$r.squared)
  expect_equivalent(lmro$r.squared[2], lmhp$r.squared)

  J <- length(unique(mtcars$cyl))
  n <- nrow(mtcars)
  r <- 2

  # Have to manually do correction because clubSandwich uses n*ny and r*ny in place of n and r in
  # stata correction
  expect_equivalent(
    as.matrix(clubSandwich::vcovCR(lmo, cluster = mtcars$cyl, type = "CR0")) *
      ((J * (n - 1)) / ((J - 1) * (n - r))) ,
    vcov(lm_robust(cbind(mpg, hp) ~ wt, data = mtcars, clusters = cyl, se_type = "stata"))
  )

})
