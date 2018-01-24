context("Estimator - lm_lin")

test_that("Test LM Lin", {
  dat <- data.frame(
    Y = rnorm(100),
    Z = rbinom(100, 1, .5),
    X1 = rnorm(100),
    X2 = rbinom(100, 1, .5),
    cluster = sample(1:10, size = 100, replace = T)
  )
  lm_lin_out <- lm_lin(Y ~ Z, covariates = ~ X1 + X2, data = dat)
  expect_error(
    lm_lin_out <- lm_lin(Y ~ Z, data = dat, covariates = ~ X1 + X2),
    NA
  )

  expect_error(
    lm_lin(
      Y ~ Z + X1,
      covariates = ~ X2,
      data = dat
    ),
    "one variable on the right-hand side"
  )

  dat2 <- dat
  dat2$Z <- rnorm(100)
  # For now we allow huge multi-valued treatments, but output is wonky
  # expect_error(
  #   lm_lin(Y ~ Z,
  #          covariates = ~ X1,
  #          data = dat2),
  #   'binary'
  # )

  expect_error(
    lm_lin(
      Y ~ Z,
      covariates = Y ~ X1,
      data = dat
    ),
    "right-sided formula"
  )

  # Now allows multi-valued treatments
  expect_error(
    lm_lin(
      Y ~ factor(cluster),
      ~ X1,
      data = dat
    ),
    NA
  )

  expect_error(
    lm_lin(Y ~ treat, dat$X1, data = dat),
    "must be specified as a formula"
  )

  expect_error(
    lm_lin(Y ~ treat, ~ 1, data = dat),
    "variable on the right-hand side"
  )


  # works with one covar
  expect_error(
    lm_lin(
      Y ~ Z,
      covariates = ~ X1,
      data = dat
    ),
    NA
  )

  dat$X1_bar <- dat$X1 - mean(dat$X1)
  dat$X2_bar <- dat$X2 - mean(dat$X2)

  lm_rob_out <- lm_robust(Y ~ Z + Z * X1_bar + Z * X2_bar, data = dat)

  expect_identical(
    tidy(lm_lin_out),
    tidy(lm_rob_out)
  )

  expect_equivalent(
    lm_lin_out$coefficients,
    lm(Y ~ Z + Z * X1_bar + Z * X2_bar, data = dat)$coefficients
  )


  ## Works with multi-valued treatments
  dat$Z_mult <- rep_len(1:3, length.out = 100)

  test_mult <- function(treatment_type, dat) {
    dat$Z_mult <-
      switch(treatment_type,
        int = rep_len(1:3, length.out = 100),
        num = rep_len(c(1.5, 2.5, 5), length.out = 100),
        char = rep_len(letters[1:4], length.out = 100),
        bin = rbinom(100, 1, 0.5)
      )

    mult_out <- lm_lin(
      Y ~ Z_mult,
      covariates = ~ X1 + X2,
      data = dat
    )
    fact_mult_out <- lm_lin(
      Y ~ factor(Z_mult),
      covariates = ~ X1 + X2,
      data = dat
    )
    noint_mult_out <- lm_lin(
      Y ~ Z_mult + 0,
      covariates = ~ X1 + X2,
      data = dat
    )
    noint_fact_mult_out <- lm_lin(
      Y ~ factor(Z_mult) + 0,
      covariates = ~ X1 + X2,
      data = dat
    )

    expect_identical(
      tidy(mult_out)[, -1],
      tidy(fact_mult_out)[, -1]
    )

    rob_fact_mult_out <- lm_robust(Y ~ factor(Z_mult) * X1_bar + factor(Z_mult) * X2_bar, data = dat)

    expect_identical(
      tidy(fact_mult_out),
      tidy(rob_fact_mult_out)
    )

    # Also works with no intercept!
    expect_identical(
      tidy(noint_mult_out)[, -1],
      tidy(noint_fact_mult_out)[, -1]
    )
  }

  test_mult("int", dat)
  test_mult("num", dat)
  test_mult("char", dat)
  # test_mult("bin", dat) this gives weird answers because it isn't really valid when not a factor!

  ## Works with missing data in treatment
  dat$Z[23] <- NA
  dat$X1_bar <- dat$X1 - mean(dat$X1[-23])
  dat$X2_bar <- dat$X2 - mean(dat$X2[-23])

  expect_equal(
    tidy(lm_lin(
      Y ~ Z,
      covariates = ~ X1 + X2,
      data = dat
    )),
    tidy(lm_robust(
      Y ~ Z + Z * X1_bar + Z * X2_bar,
      data = dat
    ))
  )

  ## Test cluster passes through
  expect_identical(
    tidy(lm_lin(
      Y ~ Z,
      covariates = ~ X1 + X2,
      data = dat,
      clusters = cluster
    )),
    tidy(lm_robust(
      Y ~ Z + Z * X1_bar + Z * X2_bar,
      data = dat,
      clusters = cluster
    ))
  )

  ## Test that it works with subset
  keep <- setdiff(which(dat$Y > 0), 23)
  dat$X1_bar <- dat$X1 - mean(dat$X1[keep])
  dat$X2_bar <- dat$X2 - mean(dat$X2[keep])

  expect_equal(
    tidy(lm_lin(
      Y ~ Z,
      covariates = ~ X1 + X2,
      data = dat,
      clusters = cluster,
      subset = Y > 0
    )),
    tidy(lm_robust(
      Y ~ Z + Z * X1_bar + Z * X2_bar,
      data = dat,
      clusters = cluster,
      subset = Y > 0
    ))
  )

  # Works with factors
  dat <- data.frame(
    Y = rnorm(100),
    treat = as.factor(rbinom(100, 1, .5)),
    X1 = rnorm(100),
    X2 = as.factor(rbinom(100, 1, .5)),
    cluster = sample(1:10, size = 100, replace = T)
  )

  dat$X1_bar <- dat$X1 - mean(dat$X1)
  dat$X21_bar <- as.numeric(dat$X2 == 1) - mean(dat$X2 == 1)

  expect_equivalent(
    tidy(lm_lin(
      Y ~ treat,
      covariates = ~ X1 + X2,
      data = dat,
      clusters = cluster
    )),
    tidy(lm_robust(
      Y ~ treat + treat * X1_bar + treat * X21_bar,
      data = dat,
      clusters = cluster
    ))
  )

  ## Works with a factor with spaces in the name (often the case for clusters)
  dat$X2 <- as.factor(sample(
    c("This is a level", "Get on my level"),
    size = 100,
    replace = T
  ))
  ## for lm_robust
  dat$X2_bar <-
    as.numeric(dat$X2 == "This is a level") - mean(dat$X2 == "This is a level")

  ## Names will differ
  expect_equivalent(
    tidy(
      lm_lin(
        Y ~ treat,
        covariates = ~ X1 + X2,
        data = dat,
        clusters = cluster
      )
    )[, -1],
    tidy(
      lm_robust(
        Y ~ treat + treat * X1_bar + treat * X2_bar,
        data = dat,
        clusters = cluster
      )
    )[, -1]
  )

  ## Works with missingness on cluster
  dat$cluster[40] <- NA
  dat$X1_bar <- dat$X1 - mean(dat$X1[-40])
  dat$X2_bar <-
    as.numeric(dat$X2 == "This is a level") - mean(dat$X2[-40] == "This is a level")

  expect_warning(
    lin_out <- lm_lin(
      Y ~ treat,
      covariates = ~ X1 + X2,
      data = dat,
      clusters = cluster
    ),
    "missingness in the cluster"
  )

  expect_warning(
    rob_out <- lm_robust(
      Y ~ treat + treat * X1_bar + treat * X2_bar,
      data = dat,
      clusters = cluster
    ),
    "missingness in the cluster"
  )

  # drop coefficient name because names will differ
  expect_equivalent(
    tidy(lin_out)[, -1],
    tidy(rob_out)[, -1]
  )

  # rank deficient cases
  dat$treat2 <- dat$treat
  dat$X1_2 <- dat$X1

  lm_lin(Y ~ treat, ~ treat2 + X1, data = dat) # somewhat odd behavior
  expect_equivalent(
    is.na(lm_lin(Y ~ treat, ~ X1_2 + X1, data = dat)$coefficients),
    c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE)
  )

  # Does lm_lin parse covariate names correctly?
  # Previously it dropped the factor(), now properly builds factor
  lmlo <- lm_lin(Y ~ treat, ~factor(cluster) + X1, data = dat)
  expect_equal(
    nrow(tidy(lmlo)),
    22
  )

  # works with a binary with no intercept (bit odd, though!)
  dat$z <- rbinom(nrow(dat), 1, 0.5)
  lmlo <- lm_lin(Y ~ z + 0, ~X1, data = dat)
  expect_equal(
    lmlo$coefficient_name,
    c("z", "X1_bar", "z:X1_bar")
  )
})
