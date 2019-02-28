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
    "must only have the treatment variable on the right-hand side of the formula"
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

  dat$X1_c <- dat$X1 - mean(dat$X1)
  dat$X2_c <- dat$X2 - mean(dat$X2)

  lm_rob_out <- lm_robust(Y ~ Z + Z * X1_c + Z * X2_c, data = dat)

  expect_equal(
    tidy(lm_lin_out),
    tidy(lm_rob_out)
  )

  expect_equivalent(
    coef(lm_lin_out),
    coef(lm(Y ~ Z + Z * X1_c + Z * X2_c, data = dat))
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

    expect_equal(
      tidy(mult_out)[, -1],
      tidy(fact_mult_out)[, -1]
    )

    rob_fact_mult_out <- lm_robust(Y ~ factor(Z_mult) * X1_c + factor(Z_mult) * X2_c, data = dat)

    expect_equal(
      tidy(fact_mult_out),
      tidy(rob_fact_mult_out)
    )

    # Also works with no intercept!
    expect_equal(
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
  dat$X1_c <- dat$X1 - mean(dat$X1[-23])
  dat$X2_c <- dat$X2 - mean(dat$X2[-23])

  expect_equal(
    tidy(lm_lin(
      Y ~ Z,
      covariates = ~ X1 + X2,
      data = dat
    )),
    tidy(lm_robust(
      Y ~ Z + Z * X1_c + Z * X2_c,
      data = dat
    ))
  )

  ## Test cluster passes through
  expect_equal(
    tidy(lm_lin(
      Y ~ Z,
      covariates = ~ X1 + X2,
      data = dat,
      clusters = cluster
    )),
    tidy(lm_robust(
      Y ~ Z + Z * X1_c + Z * X2_c,
      data = dat,
      clusters = cluster
    ))
  )

  ## Test that it works with subset
  keep <- setdiff(which(dat$Y > 0), 23)
  dat$X1_c <- dat$X1 - mean(dat$X1[keep])
  dat$X2_c <- dat$X2 - mean(dat$X2[keep])

  expect_equal(
    tidy(lm_lin(
      Y ~ Z,
      covariates = ~ X1 + X2,
      data = dat,
      clusters = cluster,
      subset = Y > 0
    )),
    tidy(lm_robust(
      Y ~ Z + Z * X1_c + Z * X2_c,
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

  dat$X1_c <- dat$X1 - mean(dat$X1)
  dat$X21_c <- as.numeric(dat$X2 == 1) - mean(dat$X2 == 1)

  expect_equivalent(
    tidy(lm_lin(
      Y ~ treat,
      covariates = ~ X1 + X2,
      data = dat,
      clusters = cluster
    )),
    tidy(lm_robust(
      Y ~ treat + treat * X1_c + treat * X21_c,
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
  dat$X2_c <-
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
        Y ~ treat + treat * X1_c + treat * X2_c,
        data = dat,
        clusters = cluster
      )
    )[, -1]
  )

  ## Works with missingness on cluster
  dat$cluster[40] <- NA
  dat$X1_c <- dat$X1 - mean(dat$X1[-40])
  dat$X2_c <-
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
      Y ~ treat + treat * X1_c + treat * X2_c,
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
    is.na(coef(lm_lin(Y ~ treat, ~ X1_2 + X1, data = dat))),
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
    lmlo$term,
    c("z", "X1_c", "z:X1_c")
  )

  # behaves correctly with "odd" covariates (see issue #283)
  lmlo <- lm_lin(Y ~ z, ~ is.na(X1), data = dat)
  expect_equal(
    lmlo$term,
    c("(Intercept)", "z", "(is.na(X1)TRUE)_c", "z:(is.na(X1)TRUE)_c")
  )
})

test_that("lm_lin same as sampling perspective", {

  # Unweighted matches sampling view
  lmo <- lm_lin(mpg ~ am, ~ hp, data = mtcars)
  m_hp <- mean(mtcars$hp)
  areg <- lm(mpg ~ hp, data = mtcars, subset = am == 1)
  breg <- lm(mpg ~ hp, data = mtcars, subset = am == 0)
  ate <-
    with(mtcars[mtcars$am == 1, ],
         mean(mpg) + (m_hp - mean(hp)) * coef(areg)[2]) -
    with(mtcars[mtcars$am == 0, ],
         mean(mpg) + (m_hp - mean(hp)) * coef(breg)[2])

  expect_equivalent(
    ate,
    coef(lmo)["am"]
  )
})

test_that("weighted lm_lin same as with one covar sampling view", {

  # Weighted matches (one covar)
  lmwo <- lm_lin(mpg ~ am, ~ hp, weights = wt, data = mtcars)
  hp_wmean <- weighted.mean(mtcars$hp, mtcars$wt)
  wareg <- lm(mpg ~ hp, data = mtcars, subset = am == 1, weights = wt)
  wbreg <- lm(mpg ~ hp, data = mtcars, subset = am == 0, weights = wt)
  wate <-
    with(mtcars[mtcars$am == 1, ],
         weighted.mean(mpg, wt) + (hp_wmean - weighted.mean(hp, wt)) * coef(wareg)[2]) -
    with(mtcars[mtcars$am == 0, ],
         weighted.mean(mpg, wt) + (hp_wmean - weighted.mean(hp, wt)) * coef(wbreg)[2])

  expect_equivalent(
    wate,
    coef(lmwo)["am"]
  )
})

test_that("weighted lm_lin same as with two covar sampling view", {

  # Weighted matches (two covars)
  lmw2o <- lm_lin(mpg ~ am, ~ hp + cyl, weights = wt, data = mtcars)
  hpcyl_wmean <- apply(mtcars[, c("hp", "cyl")], 2, weighted.mean, mtcars$wt)
  w2areg <- lm(mpg ~ hp + cyl, data = mtcars, subset = am == 1, weights = wt)
  w2breg <- lm(mpg ~ hp + cyl, data = mtcars, subset = am == 0, weights = wt)
  w2ate <-
    with(mtcars[mtcars$am == 1, ],
         weighted.mean(mpg, wt) + (hpcyl_wmean - apply(cbind(hp, cyl), 2, weighted.mean, wt)) %*% coef(w2areg)[2:3]) -
    with(mtcars[mtcars$am == 0, ],
         weighted.mean(mpg, wt) + (hpcyl_wmean - apply(cbind(hp, cyl), 2, weighted.mean, wt)) %*% coef(w2breg)[2:3])

  expect_equivalent(
    w2ate,
    coef(lmw2o)["am"]
  )
})

test_that("lm_lin properly renames trickily named variables", {

  # lm_lin should add parentheses around variables that have colons in their name
  # or that have parentheses in the name that are not in the first position
  lo <- lm_lin(mpg ~ am, ~ wt*cyl + log(wt), mtcars)

  expect_equal(
    lo$term,
    c("(Intercept)", "am", "wt_c", "cyl_c", "(log(wt))_c", "(wt:cyl)_c",
      "am:wt_c", "am:cyl_c", "am:(log(wt))_c", "am:(wt:cyl)_c")
  )
})

test_that("lm_lin works with multiple outcomes", {

  lmpg <- lm_lin(mpg ~ am, ~ cyl, mtcars)
  lwt <- lm_lin(wt ~ am, ~ cyl, mtcars)
  lboth <- lm_lin(cbind(mpg, wt) ~ am, ~ cyl, mtcars)

  expect_equivalent(
    tidy(lmpg),
    tidy(lboth)[1:4, ]
  )

  expect_equivalent(
    tidy(lwt),
    tidy(lboth)[5:8, ]
  )

  expect_equivalent(
    lboth$fstatistic[1:2],
    c(lmpg$fstatistic[1], lwt$fstatistic[1])
  )
})
