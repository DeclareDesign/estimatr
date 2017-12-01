context("lm lin")

test_that("Test LM Lin",{

  dat <- data.frame(Y = rnorm(100),
                   Z = rbinom(100, 1, .5),
                   X1 = rnorm(100),
                   X2 = rbinom(100, 1, .5),
                   cluster = sample(1:10, size = 100, replace = T))
  lm_lin_out <- lm_lin(Y ~ Z, data = dat, covariates = ~ X1 + X2)
  expect_error(
    lm_lin_out <- lm_lin(Y ~ Z, data = dat, covariates = ~ X1 + X2),
    NA
  )

  expect_error(
    lm_lin(Y ~ Z + X1,
           covariates = ~ X2,
           data = dat),
    'right-hand side'
  )

  dat2 <- dat
  dat2$Z <- rnorm(100)
  expect_error(
    lm_lin(Y ~ Z,
           covariates = ~ X1,
           data = dat2),
    'binary'
  )

  expect_error(
    lm_lin(Y ~ Z,
           covariates = Y ~ X1,
           data = dat),
    'right-hand sided equation'
  )

  expect_error(
    lm_lin(Y ~ cluster,
           ~ X1,
           data = dat),
    'no more than two levels'
  )

  # works with one covar
  expect_error(
    lm_lin(Y ~ Z,
           covariates = ~ X1,
           data = dat),
    NA
  )

  dat$X1_bar <- dat$X1 - mean(dat$X1)
  dat$X2_bar <- dat$X2 - mean(dat$X2)

  lm_rob_out <- lm_robust(Y ~ Z + Z*X1_bar + Z*X2_bar, data = dat)

  expect_identical(
    summary(lm_lin_out),
    summary(lm_rob_out)
  )

  expect_equivalent(
    lm_lin_out$est,
    lm(Y ~ Z + Z*X1_bar + Z*X2_bar, data = dat)$coefficients
  )


  ## Works with missing data in treatment
  dat$Z[23] <- NA
  dat$X1_bar <- dat$X1 - mean(dat$X1[-23])
  dat$X2_bar <- dat$X2 - mean(dat$X2[-23])

  expect_identical(
    summary(lm_lin(Y ~ Z,
           covariates = ~ X1 + X2,
           data = dat)),
    summary(lm_robust(Y ~ Z + Z*X1_bar + Z*X2_bar,
              data = dat))
  )

  ## Works with no intercept
  expect_identical(
    summary(lm_lin(Y ~ Z + 0,
           covariates = ~ X1 + X2,
           data = dat)),
    summary(lm_robust(Y ~ Z + Z*X1_bar + Z*X2_bar + 0,
              data = dat))
  )

  ## Test cluster passes through
  expect_identical(
    summary(lm_lin(Y ~ Z,
           covariates = ~ X1 + X2,
           data = dat,
           clusters = cluster)),
    summary(lm_robust(Y ~ Z + Z*X1_bar + Z*X2_bar,
              data = dat,
              clusters = cluster))
  )

  ## Test that it works with subset
  keep <- setdiff(which(dat$Y > 0), 23)
  dat$X1_bar <- dat$X1 - mean(dat$X1[keep])
  dat$X2_bar <- dat$X2 - mean(dat$X2[keep])

  expect_identical(
    summary(lm_lin(Y ~ Z,
           covariates = ~ X1 + X2,
           data = dat,
           clusters = cluster,
           subset = Y > 0)),
    summary(lm_robust(Y ~ Z + Z*X1_bar + Z*X2_bar,
              data = dat,
              clusters = cluster,
              subset = Y > 0))
  )

  # Works with factors
  dat <- data.frame(Y = rnorm(100),
                   treat = as.factor(rbinom(100, 1, .5)),
                   X1 = rnorm(100),
                   X2 = as.factor(rbinom(100, 1, .5)),
                   cluster = sample(1:10, size = 100, replace = T))

  dat$X1_bar <- dat$X1 - mean(dat$X1)
  dat$X21_bar <- as.numeric(dat$X2==1) - mean(dat$X2 == 1)

  expect_equivalent(
    summary(lm_lin(Y ~ treat,
           covariates = ~ X1 + X2,
           data = dat,
           clusters = cluster)),
    summary(lm_robust(Y ~ treat + treat*X1_bar + treat*X21_bar,
              data = dat,
              clusters = cluster))
  )

  ## Works with a factor with spaces in the name (often the case for clusters)
  dat$X2 <- as.factor(sample(c("This is a level", "Get on my level"),
                            size = 100,
                            replace = T))
  ## for lm_robust
  dat$X2_bar <-
    as.numeric(dat$X2=="This is a level") - mean(dat$X2 == "This is a level")

  ## Names will differ
  expect_equivalent(
    tidy(
      lm_lin(Y ~ treat,
             covariates = ~ X1 + X2,
             data = dat,
             clusters = cluster)
    )[, -1],
    tidy(
      lm_robust(Y ~ treat + treat*X1_bar + treat*X2_bar,
                data = dat,
                clusters = cluster)
    )[, -1]
  )

  ## Works with missingness on cluster
  dat$cluster[40] <- NA
  dat$X1_bar <- dat$X1 - mean(dat$X1[-40])
  dat$X2_bar <-
    as.numeric(dat$X2=="This is a level") - mean(dat$X2[-40] == "This is a level")

  expect_warning(
    lin_out <- lm_lin(Y ~ treat,
                      covariates = ~ X1 + X2,
                      data = dat,
                      clusters = cluster),
    'missingness in the cluster'
  )

  expect_warning(
    rob_out <- lm_robust(Y ~ treat + treat*X1_bar + treat*X2_bar,
                         data = dat,
                         clusters = cluster),
    'missingness in the cluster'
  )

  # drop coefficient name because names will differ
  expect_equivalent(
    tidy(lin_out)[, -1],
    tidy(rob_out)[, -1]
  )
})
