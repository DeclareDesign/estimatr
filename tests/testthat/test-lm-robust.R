context("Estimator - lm_robust, non-clustered")

test_that("lm robust se", {
  set.seed(42)
  N <- 40
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5), X = rnorm(N), W = runif(N))

  tidy(lm_robust(Y ~ Z, data = dat))

  lm_robust(Y ~ Z, se_type = "none", data = dat)

  lm_robust(Y ~ Z + X, data = dat)
  lm_robust(Y ~ Z * X, data = dat)

  expect_equivalent(
    coef(lm_robust(Y ~ 1, data = dat))[1],
    mean(dat$Y)
  )

  expect_error(
    lm_robust(Y ~ Z + X, data = dat, se_type = "not_a_real_one"),
    "`se_type` must be either 'HC0', 'HC1', 'stata', 'HC2', 'HC3',"
  )

  # Works with subset
  lmsub <- lm_robust(Y ~ Z + X, data = dat, subset = W > 0.5)
  lmbool <- lm_robust(Y ~ Z + X, data = dat[dat$W > 0.5, ])
  expect_equal(
    rmcall(lmsub),
    rmcall(lmbool)
  )

  lm_robust(Y ~ Z, weights = W, data = dat)

  # matches.
  # commarobust::commarobust(lm(Y ~ Z, weights = W, data = dat))

  # To easily do with and without weights
  test_lm_robust_variance <- function(w) {
    # Test other estimators
    lm_hc0 <- lm_robust(Y ~ Z + X, data = dat, weights = w, se_type = "HC0")
    lm_hc1 <- lm_robust(Y ~ Z + X, data = dat, weights = w, se_type = "HC1")
    lm_hc2 <- lm_robust(Y ~ Z + X, data = dat, weights = w, se_type = "HC2")
    lm_hc3 <- lm_robust(Y ~ Z + X, data = dat, weights = w, se_type = "HC3")
    lm_stata <- lm_robust(Y ~ Z + X, data = dat, weights = w, se_type = "stata")

    # Stata is the same as HC1
    expect_equal(
      rmcall(lm_hc1),
      rmcall(lm_stata)
    )

    expect_false(all(lm_hc0$std.error == lm_hc1$std.error))
    expect_false(all(lm_hc0$std.error == lm_hc2$std.error))
    expect_false(all(lm_hc0$std.error == lm_hc3$std.error))
    expect_false(all(lm_hc1$std.error == lm_hc2$std.error))
    expect_false(all(lm_hc1$std.error == lm_hc3$std.error))
    expect_false(all(lm_hc2$std.error == lm_hc3$std.error))

    expect_equivalent(
      lm_hc0$df,
      lm_hc1$df,
      lm_hc2$df,
      lm_hc3$df,
      lm_stata$df
    )

    expect_equivalent(
      lm_hc0$std.error ^ 2,
      lm_hc1$std.error ^ 2 * ((N - length(coef(lm_hc1))) / N)
    )
  }

  # No weights first
  test_lm_robust_variance(NULL)
  test_lm_robust_variance(dat$W)

  # works with formula in a variable (always worked)
  form <- Y ~ Z
  lm_form <- lm_robust(form, data = dat)

  # works with formula inside a function (didn't work before 0.4.0)
  f <- function(data) {
    form2 <- Y ~ Z
    return(lm_robust(form2, data = data))
  }
  lm_f_form <- f(dat)

  expect_equal(
    rmcall(lm_form),
    rmcall(lm_f_form)
  )

  # Drops unused levels appropriately
  dat$Z <- as.factor(sample(LETTERS[1:3], nrow(dat), replace = TRUE))
  lmall <- lm_robust(Y ~ Z, data = dat)
  lm1 <- lm_robust(Y ~ Z, data = dat[dat$Z %in% c("A", "B"), ])
  lm2 <- lm_robust(Y ~ Z, data = dat, subset = Z %in% c("A", "B"))

  expect_equal(
    rmcall(lm1),
    rmcall(lm2)
  )

  # pvals and cis diff because dof are diff
  expect_equal(
    tidy(lmall)[1:2, 1:3],
    tidy(lm1)[, 1:3]
  )

  # rlang works
  my_w_vec <- rlang::sym("W")
  expect_equal(
    tidy(lm_robust(Y ~ Z + X, data = dat, weights = !!my_w_vec, se_type = "HC2")),
    tidy(lm_robust(Y ~ Z + X, data = dat, weights = W, se_type = "HC2"))
  )

  my_dat <- rlang::sym("dat")
  expect_equal(
    tidy(lm_robust(Y ~ Z + X, data = !!my_dat, weights = W, se_type = "HC2")),
    tidy(lm_robust(Y ~ Z + X, data = dat, weights = W, se_type = "HC2"))
  )

  my_y <- rlang::sym("Y")
  expect_equal(
    tidy(lm_robust(!!my_y ~ Z + X, data = dat, weights = W, se_type = "HC2")),
    tidy(lm_robust(Y ~ Z + X, data = dat, weights = W, se_type = "HC2"))
  )

  my_formula <- Y ~ Z + X
  expect_equal(
    tidy(lm_robust(!!my_formula, data = dat, weights = W, se_type = "HC2")),
    tidy(lm_robust(Y ~ Z + X, data = dat, weights = W, se_type = "HC2"))
  )
})

test_that("lm robust F-tests are correct", {

  co <- lm_robust(mpg ~ hp + am, data = mtcars, se_type = "classical")
  caro <- car::linearHypothesis(co, c("hp = 0", "am = 0"), test = "F")
  carolm <- car::linearHypothesis(lm(mpg ~ hp + am, data = mtcars),
                                  c("hp = 0", "am = 0"),
                                  test = "F")
  expect_equivalent(
    co$fstatistic,
    c(caro$F[2], caro$Df[2], caro$Res.Df[2])
  )
  expect_equivalent(
    co$fstatistic,
    c(carolm$F[2], carolm$Df[2], carolm$Res.Df[2])
  )

  cow <- lm_robust(mpg ~ hp + am, data = mtcars, weights = wt, se_type = "classical")
  caro <- car::linearHypothesis(cow, c("hp = 0", "am = 0"), test = "F")
  expect_equivalent(
    cow$fstatistic,
    c(caro$F[2], caro$Df[2], caro$Res.Df[2])
  )

  hc0 <- lm_robust(mpg ~ hp + am, data = mtcars, se_type = "HC0")
  caro <- car::linearHypothesis(hc0, c("hp = 0", "am = 0"), test = "F")
  carolm <- car::linearHypothesis(lm(mpg ~ hp + am, data = mtcars),
                                  c("hp = 0", "am = 0"),
                                  test = "F",
                                  white.adjust = "hc0")
  expect_equivalent(
    hc0$fstatistic,
    c(caro$F[2], caro$Df[2], caro$Res.Df[2])
  )
  expect_equivalent(
    hc0$fstatistic,
    c(carolm$F[2], carolm$Df[2], carolm$Res.Df[2])
  )

  hc1 <- lm_robust(mpg ~ hp + am, data = mtcars, se_type = "HC1")
  caro <- car::linearHypothesis(hc1, c("hp = 0", "am = 0"), test = "F")
  carolm <- car::linearHypothesis(lm(mpg ~ hp + am, data = mtcars),
                                  c("hp = 0", "am = 0"),
                                  test = "F",
                                  white.adjust = "hc1")
  expect_equivalent(
    hc1$fstatistic,
    c(caro$F[2], caro$Df[2], caro$Res.Df[2])
  )
  expect_equivalent(
    hc1$fstatistic,
    c(carolm$F[2], carolm$Df[2], carolm$Res.Df[2])
  )

  hc1w <- lm_robust(mpg ~ hp + am, data = mtcars, weights = wt, se_type = "HC1")
  caro <- car::linearHypothesis(hc1w, c("hp = 0", "am = 0"), test = "F")
  expect_equivalent(
    hc1w$fstatistic,
    c(caro$F[2], caro$Df[2], caro$Res.Df[2])
  )

  hc2 <- lm_robust(mpg ~ hp + am, data = mtcars, se_type = "HC2")
  caro <- car::linearHypothesis(hc2, c("hp = 0", "am = 0"), test = "F")
  carolm <- car::linearHypothesis(lm(mpg ~ hp + am, data = mtcars),
                                  c("hp = 0", "am = 0"),
                                  test = "F",
                                  white.adjust = "hc2")
  expect_equivalent(
    hc2$fstatistic,
    c(caro$F[2], caro$Df[2], caro$Res.Df[2])
  )
  expect_equivalent(
    hc2$fstatistic,
    c(carolm$F[2], carolm$Df[2], carolm$Res.Df[2])
  )

  hc3 <- lm_robust(mpg ~ hp + am, data = mtcars, se_type = "HC3")
  caro <- car::linearHypothesis(hc3, c("hp = 0", "am = 0"), test = "F")
  carolm <- car::linearHypothesis(lm(mpg ~ hp + am, data = mtcars),
                                  c("hp = 0", "am = 0"),
                                  test = "F",
                                  white.adjust = "hc3")
  expect_equivalent(
    hc3$fstatistic,
    c(caro$F[2], caro$Df[2], caro$Res.Df[2])
  )
  expect_equivalent(
    hc3$fstatistic,
    c(carolm$F[2], carolm$Df[2], carolm$Res.Df[2])
  )

  cr0 <- lm_robust(mpg ~ hp + am, data = mtcars, clusters = carb, se_type = "CR0")
  caro <- clubSandwich::Wald_test(lm(mpg ~ hp + am, data = mtcars),
                                  cluster = mtcars$carb,
                                  c(FALSE, TRUE, TRUE),
                                  vcov = "CR0",
                                  test = "Naive-F")

  expect_equivalent(
    cr0$fstatistic[c(1, 3)],
    c(caro$Fstat, caro$df)
  )

  cr1s <- lm_robust(mpg ~ hp + am, data = mtcars, clusters = carb, se_type = "stata")
  caro <- clubSandwich::Wald_test(lm(mpg ~ hp + am, data = mtcars),
                                  cluster = mtcars$carb,
                                  c(FALSE, TRUE, TRUE),
                                  vcov = "CR1S",
                                  test = "Naive-F")

  expect_equivalent(
    cr1s$fstatistic[c(1, 3)],
    c(caro$Fstat, caro$df)
  )

  cr1sw <- lm_robust(mpg ~ hp + am, data = mtcars, clusters = carb, weights = wt, se_type = "stata")
  caro <- clubSandwich::Wald_test(lm(mpg ~ hp + am, weights = wt, data = mtcars),
                                  cluster = mtcars$carb,
                                  c(FALSE, TRUE, TRUE),
                                  vcov = "CR1S",
                                  test = "Naive-F")

  expect_equivalent(
    cr1sw$fstatistic[c(1, 3)],
    c(caro$Fstat, caro$df)
  )

  cr2 <- lm_robust(mpg ~ hp + am, data = mtcars, clusters = carb, se_type = "CR2")
  caro <- clubSandwich::Wald_test(lm(mpg ~ hp + am, data = mtcars),
                                  cluster = mtcars$carb,
                                  c(FALSE, TRUE, TRUE),
                                  vcov = "CR2",
                                  test = "Naive-F")

  expect_equivalent(
    cr2$fstatistic[c(1, 3)],
    c(caro$Fstat, caro$df)
  )

  cr2w <- lm_robust(mpg ~ hp + am, data = mtcars, clusters = carb, weights = wt, se_type = "CR2")
  caro <- clubSandwich::Wald_test(lm(mpg ~ hp + am, weights = wt, data = mtcars),
                                  cluster = mtcars$carb,
                                  c(FALSE, TRUE, TRUE),
                                  vcov = "CR2",
                                  test = "Naive-F")

  expect_equivalent(
    cr2w$fstatistic[c(1, 3)],
    c(caro$Fstat, caro$df)
  )

})

test_that("lm robust mlm gets right fstats", {
  cocyl <- lm_robust(cyl ~ hp + am, data = mtcars, se_type = "classical")
  compg <- lm_robust(mpg ~ hp + am, data = mtcars, se_type = "classical")
  co2 <- lm_robust(cbind(cyl, mpg) ~ hp + am, data = mtcars, se_type = "classical")

  expect_equivalent(
    co2$fstatistic[1:2],
    c(cocyl$fstatistic[1], compg$fstatistic[1])
  )

  hc2cyl <- lm_robust(cyl ~ hp + am, data = mtcars, se_type = "HC2")
  hc2mpg <- lm_robust(mpg ~ hp + am, data = mtcars, se_type = "HC2")
  hc22 <- lm_robust(cbind(cyl, mpg) ~ hp + am, data = mtcars, se_type = "HC2")

  expect_equivalent(
    hc22$fstatistic[1:2],
    c(hc2cyl$fstatistic[1], hc2mpg$fstatistic[1])
  )

  cr2cyl <- lm_robust(cyl ~ hp + am, data = mtcars, cluster = carb, se_type = "CR2")
  cr2mpg <- lm_robust(mpg ~ hp + am, data = mtcars, cluster = carb, se_type = "CR2")
  cr22 <- lm_robust(cbind(cyl, mpg) ~ hp + am, data = mtcars, cluster = carb, se_type = "CR2")

  expect_equivalent(
    cr22$fstatistic[1:2],
    c(cr2cyl$fstatistic[1], cr2mpg$fstatistic[1])
  )

  cowcyl <- lm_robust(cyl ~ hp + am, data = mtcars, weights = wt, se_type = "classical")
  cowmpg <- lm_robust(mpg ~ hp + am, data = mtcars, weights = wt, se_type = "classical")
  cow2 <- lm_robust(cbind(cyl, mpg) ~ hp + am, data = mtcars, weights = wt, se_type = "classical")

  expect_equivalent(
    cow2$fstatistic[1:2],
    c(cowcyl$fstatistic[1], cowmpg$fstatistic[1])
  )

  hc2wcyl <- lm_robust(cyl ~ hp + am, data = mtcars, weights = wt, se_type = "HC2")
  hc2wmpg <- lm_robust(mpg ~ hp + am, data = mtcars, weights = wt, se_type = "HC2")
  hc2w2 <- lm_robust(cbind(cyl, mpg) ~ hp + am, data = mtcars, weights = wt, se_type = "HC2")

  expect_equivalent(
    hc2w2$fstatistic[1:2],
    c(hc2wcyl$fstatistic[1], hc2wmpg$fstatistic[1])
  )

  cr2wcyl <- lm_robust(cyl ~ hp + am, data = mtcars, cluster = carb, weights = wt, se_type = "CR2")
  cr2wmpg <- lm_robust(mpg ~ hp + am, data = mtcars, cluster = carb, weights = wt, se_type = "CR2")
  cr2w2 <- lm_robust(cbind(cyl, mpg) ~ hp + am, data = mtcars, cluster = carb, weights = wt, se_type = "CR2")

  expect_equivalent(
    cr2w2$fstatistic[1:2],
    c(cr2wcyl$fstatistic[1], cr2wmpg$fstatistic[1])
  )

})


test_that("lm robust works with missingness", {
  dat <- data.frame(
    Y = rnorm(100),
    Z = rbinom(100, 1, .5),
    X = rnorm(100),
    W = runif(100)
  )

  dat$X[23] <- NA

  expect_equal(
    rmcall(lm_robust(Y ~ Z + X, data = dat)),
    rmcall(lm_robust(Y ~ Z + X, data = dat[-23, ]))
  )
  lm_robust(Y ~ Z + X, data = dat)
  lm_robust(Y ~ Z * X, data = dat)

  ## Outcome missingness
  dat$Y[35] <- NA

  estimatr_missout_out <- lm_robust(Y ~ Z + X, data = dat)

  lm_missout_out <- lm(Y ~ Z + X, data = dat)
  lm_missout_hc2 <- cbind(
    coef(lm_missout_out),
    sqrt(diag(sandwich::vcovHC(lm_missout_out, type = "HC2")))
  )

  expect_equivalent(
    as.matrix(tidy(estimatr_missout_out)[, c("estimate", "std.error")]),
    lm_missout_hc2
  )

  # nested DFs
  dat$Y2 <- matrix(dat$Y)
  expect_equivalent(
    tidy(lm_robust(Y ~ Z + X, data = dat))[, 1:6],
    tidy(lm_robust(Y2 ~ Z + X, data = dat))[, 1:6]
  )
})

test_that("lm_robust doesn't include aux variables when . is used", {
  n <- 10
  dat <- data.frame(y = rnorm(n), x = rnorm(n))

  # not in data.frame
  clust <- rep(1:5, each = 2)

  expect_equal(
    rmcall(lm_robust(y ~ ., clusters = clust, data = dat)),
    rmcall(lm_robust(y ~ x, clusters = clust, data = dat))
  )
})


test_that("lm robust works with weights", {
  N <- 100
  dat <- data.frame(
    Y = rnorm(N),
    Z = rbinom(N, 1, .5),
    X = rnorm(N),
    W = runif(N)
  )

  ## Make sure weighting works
  expect_error(
    estimatr_out <- lm_robust(Y ~ Z * X, weights = W, data = dat),
    NA
  )

  expect_true(
    any(grepl("Weighted", capture.output(summary(estimatr_out))))
  )

  # Compare to lm output
  lm_out <- lm(Y ~ Z * X, weights = W, data = dat)
  lmo_hc2 <- cbind(
    coef(lm_out),
    sqrt(diag(sandwich::vcovHC(lm_out, type = "HC2")))
  )

  expect_equivalent(
    as.matrix(tidy(estimatr_out)[, c("estimate", "std.error")]),
    lmo_hc2
  )

  ## Make sure weighting works with missingness
  dat$W[39] <- NA

  expect_warning(
    estimatr_miss_out <- lm_robust(Y ~ Z * X, weights = W, data = dat),
    "missing"
  )

  expect_equal(
    rmcall(estimatr_miss_out),
    rmcall(lm_robust(Y ~ Z * X, weights = W, data = dat[-39, ]))
  )

  # Compare to lm output
  lm_miss_out <- lm(Y ~ Z * X, weights = W, data = dat)
  lmo_miss_hc2 <- cbind(
    coef(lm_miss_out),
    sqrt(diag(sandwich::vcovHC(lm_miss_out, type = "HC2")))
  )

  expect_equivalent(
    as.matrix(tidy(estimatr_miss_out)[, c("estimate", "std.error")]),
    lmo_miss_hc2
  )

  expect_error(
    lm_robust(Y ~ Z, data = dat, weights = c(-0.5, runif(N - 1))),
    "`weights` must not be negative"
  )
})

test_that("lm_robust_fit adds column names", {
  n <- 10
  y <- rnorm(n)
  X <- matrix(rnorm(n * 3), ncol = 3)

  lm_o <- lm_robust_fit(
    y = y,
    X = X,
    weights = NULL,
    cluster = NULL,
    ci = TRUE,
    se_type = "classical",
    alpha = 0.05,
    return_vcov = TRUE,
    try_cholesky = TRUE,
    has_int = FALSE,
    iv_stage = list(0)
  )

  expect_equal(
    lm_o$term,
    c("X1", "X2", "X3")
  )
})

test_that("lm robust works with large data", {
  N <- 75000
  dat <- data.frame(
    Y = rbinom(N, 1, .5),
    X1 = rnorm(N),
    X2 = rnorm(N),
    X3 = rnorm(N)
  )
  expect_error(
    lm_robust(Y ~ X1 + X2 + X3, data = dat, se_type = "none"),
    NA
  )
})

set.seed(42)
N <- 100
dat <- data.frame(
  Y = rbinom(N, 1, .5),
  X1 = rnorm(N),
  X2 = rnorm(N),
  X3 = rnorm(N)
)

test_that("lm robust works with rank-deficient X", {

  dat$Z1 <- dat$X1
  sum_lm <- summary(lm(Y ~ X1 + X2 + Z1 + X3, data = dat))
  ## manually build vector of coefficients, can't extract from summary.lm
  out_sumlm <- matrix(NA, nrow = length(sum_lm$aliased), ncol = 2)
  j <- 1
  for (i in seq_along(sum_lm$aliased)) {
    if (!sum_lm$aliased[i]) {
      out_sumlm[i, ] <- coef(sum_lm)[j, 1:2]
      j <- j + 1
    }
  }

  ## order sometimes is different! Not stable order!
  # expect_equivalent(
  #   as.matrix(tidy(lm_robust(Y ~ X1 + X2 + Z1 + X3, data = dat, se_type = 'classical'))[, c('estimate', 'std.error')]),
  #   out_sumlm
  # )

  dat$Z1 <- dat$X1 + 5

  ## Not the same as LM! Different QR decompositions when dependency isn't just equivalency
  expect_equivalent(
    as.matrix(tidy(lm_robust(Y ~ X1 + X2 + Z1 + X3, data = dat, se_type = "classical"))[, c("estimate", "std.error")]),
    as.matrix(RcppEigen:::summary.fastLm(RcppEigen::fastLm(Y ~ X1 + X2 + Z1 + X3, data = dat))$coefficients[, 1:2])
  )

  # trigger cascade to QR from try_chol; set seed above because try_cholesky
  # sometimes will work!
  expect_equivalent(
    tidy(lm_robust(Y ~ X1 + X2 + Z1 + X3, data = dat)),
    tidy(lm_robust(Y ~ X1 + X2 + Z1 + X3, data = dat, try_cholesky = TRUE))
  )

  # Weighted rank deficient
  dat$w <- 1
  expect_equivalent(
    tidy(lm_robust(Y ~ X1 + X2 + Z1 + X3, data = dat)),
    tidy(lm_robust(Y ~ X1 + X2 + Z1 + X3, data = dat, weights = w))
  )

  expect_true(
    any(grepl(
      "not defined because the design matrix is rank deficient",
      capture.output(summary(lm_robust(Y ~ X1 + X2 + Z1 + X3, data = dat)))
    ))
  )
})

test_that("r squared is right", {
  lmo <- summary(lm(mpg ~ hp, mtcars))
  lmow <- summary(lm(mpg ~ hp, mtcars, weights = wt))
  lmon <- summary(lm(mpg ~ hp - 1, mtcars))
  lmown <- summary(lm(mpg ~ hp - 1, mtcars, weights = wt))

  lmro <- lm_robust(mpg ~ hp, mtcars)
  lmrow <- lm_robust(mpg ~ hp, mtcars, weights = wt)
  lmroclust <- lm_robust(mpg ~ hp, mtcars, clusters = carb)
  lmrowclust <- lm_robust(mpg ~ hp, mtcars, weights = wt, clusters = carb)
  lmron <- lm_robust(mpg ~ hp - 1, mtcars)
  lmrown <- lm_robust(mpg ~ hp - 1, mtcars, weights = wt)
  lmrclust <- lm_robust(mpg ~ hp - 1, mtcars, clusters = carb)
  lmrwclust <- lm_robust(mpg ~ hp - 1, mtcars, weights = wt, clusters = carb)

  # Use equivalent instead of equal because we change the name of the fstat value
  expect_equivalent(
    c(lmo$r.squared, lmo$adj.r.squared),
    c(lmro$r.squared, lmro$adj.r.squared)
  )

  expect_equivalent(
    c(lmow$r.squared, lmow$adj.r.squared),
    c(lmrow$r.squared, lmrow$adj.r.squared)
  )

  expect_equivalent(
    c(lmon$r.squared, lmon$adj.r.squared),
    c(lmron$r.squared, lmron$adj.r.squared)
  )

  expect_equivalent(
    c(lmown$r.squared, lmown$adj.r.squared),
    c(lmrown$r.squared, lmrown$adj.r.squared)
  )

  expect_equal(
    c(lmon$r.squared, lmon$adj.r.squared),
    c(lmrclust$r.squared, lmrclust$adj.r.squared)
  )

  expect_equal(
    c(lmown$r.squared, lmown$adj.r.squared),
    c(lmrwclust$r.squared, lmrwclust$adj.r.squared)
  )

  expect_equal(
    c(lmo$r.squared, lmo$adj.r.squared),
    c(lmroclust$r.squared, lmroclust$adj.r.squared)
  )

  expect_equal(
    c(lmow$r.squared, lmow$adj.r.squared),
    c(lmrowclust$r.squared, lmrowclust$adj.r.squared)
  )
})

test_that("multiple outcomes", {

  lmo <- lm(cbind(mpg, hp) ~ cyl, data = mtcars)
  lmro <- lm_robust(cbind(mpg, hp) ~ cyl, data = mtcars, se_type = "classical")
  mo <- tidy(lmro)

  expect_identical(
    mo$term,
    c("(Intercept)", "cyl", "(Intercept)", "cyl")
  )

  expect_equal(
    coef(lmro),
    coef(lmo)
  )

  expect_equal(
    vcov(lmo),
    vcov(lmro)
  )

  expect_equal(
    sandwich::vcovHC(lmo, type = "HC0"),
    vcov(lm_robust(cbind(mpg, hp) ~ cyl, data = mtcars, se_type = "HC0"))
  )
  expect_equal(
    sandwich::vcovHC(lmo, type = "HC1"),
    vcov(lm_robust(cbind(mpg, hp) ~ cyl, data = mtcars, se_type = "HC1"))
  )
  expect_equal(
    sandwich::vcovHC(lmo, type = "HC2"),
    vcov(lm_robust(cbind(mpg, hp) ~ cyl, data = mtcars, se_type = "HC2"))
  )
  expect_equal(
    sandwich::vcovHC(lmo, type = "HC3"),
    vcov(lm_robust(cbind(mpg, hp) ~ cyl, data = mtcars, se_type = "HC3"))
  )

  # with weights
  lmo <- lm(cbind(mpg, hp) ~ cyl, data = mtcars, weights = wt)
  lmro <- lm_robust(cbind(mpg, hp) ~ cyl, data = mtcars, weights = wt, se_type = "classical")
  mo <- tidy(lmro)

  expect_identical(
    mo$term,
    c("(Intercept)", "cyl", "(Intercept)", "cyl")
  )

  expect_equal(
    coef(lmro),
    coef(lmo)
  )

  expect_equivalent(
    sapply(summary(lmo)[[1]][c("r.squared", "adj.r.squared", "fstatistic")], `[`, 1),
    sapply(lmro[c("r.squared", "adj.r.squared", "fstatistic")], `[`, 1)
  )
  expect_equivalent(
    sapply(summary(lmo)[[2]][c("r.squared", "adj.r.squared", "fstatistic")], `[`, 1),
    sapply(lmro[c("r.squared", "adj.r.squared", "fstatistic")], `[`, 2)
  )

  # with missingness
  mtcarsmiss <- mtcars
  mtcarsmiss$hp[10] <- NA

  lmo <- lm(cbind(mpg, hp) ~ cyl, data = mtcarsmiss)
  lmro <- lm_robust(cbind(mpg, hp) ~ cyl, data = mtcarsmiss, se_type = "classical")

  expect_equivalent(
    do.call(rbind, lapply(summary(lmo), function(x) x$coefficients[, c(1, 2, 4)])),
    as.matrix(tidy(lmro)[, c("estimate", "std.error", "p.value")])
  )

  expect_equivalent(
    sapply(summary(lmo)[[1]][c("r.squared", "adj.r.squared", "fstatistic")], `[`, 1),
    sapply(lmro[c("r.squared", "adj.r.squared", "fstatistic")], `[`, 1)
  )
  expect_equivalent(
    sapply(summary(lmo)[[2]][c("r.squared", "adj.r.squared", "fstatistic")], `[`, 1),
    sapply(lmro[c("r.squared", "adj.r.squared", "fstatistic")], `[`, 2)
  )
})
