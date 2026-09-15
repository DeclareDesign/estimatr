library(estimatr)

# lm_robust's own behaviour: its defaults, how it reads a formula and missing
# values, what it returns, and the regressions filed against it by issue
# number or by the 2026-08-23 review's item.
#
# What is true of every estimator lives in the property files
# (test_invariance.R, test_equivalence.R, test_degenerate.R, test_errors.R,
# test_methods.R); correctness against other implementations lives in the
# test_vs_*.R files; absorbed fixed effects live in test_fixed_effects.R and
# test_fe_leverage.R.

lm_data <- function() {
  set.seed(42)
  n <- 200
  dat <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    z = rbinom(n, 1, 0.5),
    cl = rep(1:20, 10),
    block = rep(1:20, each = 10),
    w = runif(n, 0.5, 2)
  )
  dat
}

dat <- lm_data()
n <- nrow(dat)

# Strip what differs between two calls that fit the same model: the call itself,
# and the formula's environment.
rmcall <- function(x) {
  x$call <- NULL
  if (!is.null(x$terms)) attr(x$terms, ".Environment") <- NULL
  x
}

# ---- defaults ----

test_that("the default se_type is HC2, or CR2 with clusters", {
  expect_equal(lm_robust(y ~ x, data = dat)$se_type, "HC2")
  clustered <- lm_robust(y ~ x + z, data = dat, clusters = cl)
  expect_equal(clustered$se_type, "CR2")
  expect_equal(clustered$nclusters, 20L)
})

test_that("classical standard errors are lm()'s", {
  m0 <- lm_robust(y ~ x + z, data = dat, se_type = "classical")
  m1 <- lm(y ~ x + z, data = dat)
  expect_equal(unname(coef(m0)), unname(coef(m1)), tolerance = 1e-10)
  expect_equal(unname(m0$std.error), unname(summary(m1)$coef[, 2]), tolerance = 1e-10)
})

test_that("stata without clusters is HC1, which is HC0 scaled by n / (n - k)", {
  set.seed(42)
  N <- 40
  d <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, 0.5), X = rnorm(N))
  hc0 <- lm_robust(Y ~ Z + X, data = d, se_type = "HC0")
  hc1 <- lm_robust(Y ~ Z + X, data = d, se_type = "HC1")
  stata <- lm_robust(Y ~ Z + X, data = d, se_type = "stata")
  expect_equal(rmcall(hc1), rmcall(stata))
  k <- length(coef(hc0))
  expect_equal(hc0$std.error^2, hc1$std.error^2 * (N - k) / N, tolerance = 1e-12)
})

# ---- reading the formula and the data ----

test_that("an intercept-only model returns the mean", {
  set.seed(1)
  d <- data.frame(Y = rnorm(50))
  expect_equal(coef(lm_robust(Y ~ 1, data = d))[[1]], mean(d$Y))
})

test_that("a formula held in a variable, or built inside a function, is the same fit", {
  set.seed(1)
  d <- data.frame(Y = rnorm(40), Z = rbinom(40, 1, 0.5))
  form <- Y ~ Z
  inside <- function(dd) {
    form2 <- Y ~ Z
    lm_robust(form2, data = dd)
  }
  expect_equal(rmcall(lm_robust(form, data = d)), rmcall(inside(d)))
})

test_that("a . in the formula expands against the data, clusters included", {
  set.seed(42)
  clust <- rep(1:6, 10)
  d <- data.frame(y = rnorm(60), x = rnorm(60))
  expect_equal(rmcall(lm_robust(y ~ ., clusters = clust, data = d)),
               rmcall(lm_robust(y ~ x, clusters = clust, data = d)))
})

test_that("factor levels a subset removes are dropped", {
  set.seed(42)
  d <- data.frame(Y = rnorm(40), Z = factor(sample(LETTERS[1:3], 40, replace = TRUE)))
  expect_equal(rmcall(lm_robust(Y ~ Z, data = d[d$Z %in% c("A", "B"), ])),
               rmcall(lm_robust(Y ~ Z, data = d, subset = Z %in% c("A", "B"))))
})

test_that("a row with a missing value is the same fit as the row removed, field for field", {
  # test_equivalence.R holds the estimates of every estimator to this; here
  # every returned field of lm_robust is.
  set.seed(42)
  N <- 40
  d <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, 0.5), X = rnorm(N), W = runif(N),
                  cl = rep(1:8, 5))
  missing_y <- d
  missing_y$Y[5] <- NA
  expect_equal(rmcall(lm_robust(Y ~ Z + X, data = missing_y)),
               rmcall(lm_robust(Y ~ Z + X, data = d[-5, ])))
  missing_x <- d
  missing_x$X[23] <- NA
  expect_equal(rmcall(lm_robust(Y ~ Z + X, data = missing_x)),
               rmcall(lm_robust(Y ~ Z + X, data = d[-23, ])))
  weighted <- d
  weighted$Y[39] <- NA
  expect_equal(rmcall(lm_robust(Y ~ Z * X, weights = W, data = weighted)),
               rmcall(lm_robust(Y ~ Z * X, weights = W, data = d[-39, ])))
  expect_equal(rmcall(lm_robust(Y ~ X, clusters = cl, data = missing_y)),
               rmcall(lm_robust(Y ~ X, clusters = cl, data = d[-5, ])))
})

test_that("#421: an ordered factor cluster does not crash", {
  d <- dat
  d$cl_ord <- factor(dat$cl, ordered = TRUE)
  m <- lm_robust(y ~ x + z, data = d, clusters = cl_ord)
  expect_equal(m$nclusters, length(unique(dat$cl)))
  expect_equal(m$std.error, lm_robust(y ~ x + z, data = dat, clusters = cl)$std.error)
})

test_that("B13: a character cluster is coerced rather than handed to the C++", {
  set.seed(1)
  N <- 60
  d <- data.frame(y = rnorm(N), x = rnorm(N), cl = sample(6, N, TRUE))
  d$clch <- paste0("g", d$cl)
  for (se in c("CR0", "CR2", "stata")) {
    a <- lm_robust(y ~ x, clusters = clch, data = d, se_type = se)
    b <- lm_robust(y ~ x, clusters = cl, data = d, se_type = se)
    expect_equal(a$std.error, b$std.error, info = se)
    expect_equal(a$nclusters, 6L, info = se)
  }
})

test_that("B6: one cluster is refused rather than answered with a zero", {
  # Every cluster-robust estimator here divides by J - 1 somewhere except CR2,
  # whose Satterthwaite degrees of freedom never reach that guard, so a single
  # cluster produced standard errors of order 1e-17 and no warning.
  set.seed(1)
  N <- 60
  d <- data.frame(y = rnorm(N), x = rnorm(N), one = 1L)
  for (se in c("CR0", "CR2", "stata")) {
    expect_error(lm_robust(y ~ x, clusters = one, data = d, se_type = se),
                 "only one level", info = se)
  }
  expect_error(iv_robust(y ~ x | x, clusters = one, data = d), "only one level")
})

test_that("B5: zero weights are not observations", {
  # N was nrow(X), so ten zero weights among 100 rows gave df.residual 98
  # where lm() gives 88, and every classical, HC1 and stata standard error and
  # every p-value moved with it. A zero-weight row stays in residuals and
  # fitted.values, as lm() keeps it; it is only the counting that changes.
  set.seed(1)
  N_b5 <- 100
  d <- data.frame(x = rnorm(N_b5))
  d$y <- d$x + rnorm(N_b5)
  w_b5 <- runif(N_b5)
  w_b5[1:10] <- 0

  m <- lm_robust(y ~ x, data = d, weights = w_b5, se_type = "classical")
  l <- lm(y ~ x, data = d, weights = w_b5)
  expect_equal(m$df.residual, l$df.residual)
  expect_equal(m$nobs, nobs(l))
  expect_equal(m$std.error[["x"]],
               summary(l)$coefficients["x", "Std. Error"], tolerance = 1e-12)
  expect_equal(length(m$residuals), N_b5)

  # and the definition of a zero weight: the same answer as deleting the row
  keep <- w_b5 > 0
  for (se in c("classical", "HC0", "HC1", "HC2", "HC3")) {
    expect_equal(
      lm_robust(y ~ x, data = d, weights = w_b5, se_type = se)$std.error[["x"]],
      lm_robust(y ~ x, data = d[keep, ], weights = w_b5[keep],
                se_type = se)$std.error[["x"]],
      tolerance = 1e-12, info = se
    )
  }
})

test_that("a weighted R-squared stays in [0, 1]", {
  # Squaring the internal square-root weights a second time, which an early
  # version of 2.0 did, lets it leave the interval.
  set.seed(13)
  d <- data.frame(y = 1:10 + rnorm(10, 0, 0.01), x = 1:10, w = (1:10) / 10)
  m <- lm_robust(y ~ x, data = d, weights = w)
  expect_equal(m$r.squared, summary(lm(y ~ x, data = d, weights = w))$r.squared,
               tolerance = 1e-12)
})

test_that("B9: fitting does not advance the RNG", {
  # The hidden variable names were built with sample.int(), so every fit moved
  # the seed. This package lives inside DeclareDesign simulation loops, where
  # that changes what the next draw is.
  set.seed(1)
  N_b9 <- 50
  d <- data.frame(y = rnorm(N_b9), x = rnorm(N_b9), w = runif(N_b9),
                  bl = sample(5, N_b9, TRUE), cl = sample(8, N_b9, TRUE))

  # The fits themselves are beside the point here; only whether they move the
  # seed is.
  draw_after <- function(expr) {
    set.seed(99)
    suppressWarnings(force(expr))
    rnorm(1)
  }
  baseline <- draw_after(NULL)
  expect_equal(draw_after(lm_robust(y ~ x, data = d)), baseline)
  expect_equal(draw_after(lm_robust(y ~ x, weights = w, fixed_effects = ~ bl,
                                    clusters = cl, data = d)), baseline)
  expect_equal(draw_after(lm_lin(y ~ x, covariates = ~ w, data = d)), baseline)
})

test_that("B11: an offset() term is refused rather than ignored", {
  # It was parsed and never read, so the fit came back as though the term were
  # absent: a silently different model, not a refused one.
  set.seed(1)
  N_b11 <- 100
  d <- data.frame(x = rnorm(N_b11), off = rnorm(N_b11))
  d$y <- d$x + 2 * d$off + rnorm(N_b11)
  expect_error(lm_robust(y ~ x + offset(off), data = d), "not supported")

  # the rewrite the message names is exact
  expect_equal(coef(lm_robust(I(y - off) ~ x, data = d))[["x"]],
               coef(lm(y ~ x + offset(off), data = d))[["x"]])
})

# ---- residuals (estimatr #345) ----

test_that("#345: lm_robust returns residuals on the scale of the data", {
  m <- lm_robust(y ~ x + z, data = dat)
  expect_equal(length(m$residuals), n)
  expect_equal(unname(m$residuals), unname(residuals(lm(y ~ x + z, data = dat))),
               tolerance = 1e-12)
  # estimatr registers no residuals method: the default one reads the field.
  # Named explicitly, because clubSandwich registers a `residuals.lm_robust` of
  # its own, so merely loading that namespace anywhere in the suite changes
  # which method a bare `residuals(m)` would be testing.
  expect_equal(getS3method("residuals", "default")(m), m$residuals)
  expect_equal(unname(m$residuals + m$fitted.values), dat$y, tolerance = 1e-12)
})

test_that("#345: residuals with weights are on the unweighted scale", {
  m <- lm_robust(y ~ x + z, data = dat, weights = w)
  expect_equal(unname(m$residuals),
               unname(residuals(lm(y ~ x + z, data = dat, weights = dat$w))),
               tolerance = 1e-12)
  expect_equal(unname(m$residuals + m$fitted.values), dat$y, tolerance = 1e-12)
})

test_that("#345: residuals with clusters come back in the original row order", {
  # prep_data sorts rows by cluster internally; cl is interleaved here, so an
  # unsorted return would be silently misaligned with the data.
  m <- lm_robust(y ~ x + z, data = dat, clusters = cl)
  expect_equal(unname(m$residuals), unname(residuals(lm(y ~ x + z, data = dat))),
               tolerance = 1e-12)
  expect_equal(unname(m$residuals + m$fitted.values), dat$y, tolerance = 1e-12)
})

test_that("#345: residuals with fixed effects are the full-model residuals", {
  m <- lm_robust(y ~ x + z, data = dat, fixed_effects = ~ block)
  m_dummy <- lm(y ~ x + z + factor(block), data = dat)
  expect_equal(unname(m$residuals), unname(residuals(m_dummy)), tolerance = 1e-10)
  expect_equal(unname(m$residuals + m$fitted.values), dat$y, tolerance = 1e-10)
})

# ---- collinearity warnings (estimatr #411) ----

test_that("#411: dropped collinear regressors warn and name themselves", {
  d <- dat
  d$x_copy <- d$x
  expect_warning(lm_robust(y ~ x + x_copy, data = d), "collinear")
  expect_warning(lm_robust(y ~ x + x_copy, data = d), "x_copy")
  m <- suppressWarnings(lm_robust(y ~ x + x_copy, data = d))
  expect_true(is.na(m$coefficients[["x_copy"]]))
})

test_that("#411: no warning when the design matrix is full rank", {
  expect_no_warning(lm_robust(y ~ x + z, data = dat))
  expect_no_warning(lm_robust(y ~ x + z, data = dat, clusters = cl))
  expect_no_warning(lm_lin(y ~ z, covariates = ~ x, data = dat))
})

# ---- predict (estimatr #403, #404) ----

test_that("#403: predict() with no newdata returns the in-sample fit", {
  m <- lm_robust(y ~ x + z, data = dat)
  expect_equal(predict(m), m$fitted.values)
  expect_error(predict(m, se.fit = TRUE), "newdata")
})

test_that("#404: predict() works with fixed effects, with and without factors", {
  d <- dat
  d$f <- factor(rep(letters[1:4], length.out = n))
  m <- lm_robust(y ~ x + f, data = d, fixed_effects = ~ block)
  nd <- d[1:5, ]
  # equals the dummy-variable regression, which is the definition of correct
  expect_equal(unname(predict(m, nd)),
               unname(predict(lm(y ~ x + f + factor(block), data = d), nd)),
               tolerance = 1e-8)
  # and reproduces the in-sample fit
  expect_equal(unname(predict(m, d)), unname(m$fitted.values), tolerance = 1e-8)
})

test_that("#404: predict() rejects new FE levels and multi-way FE", {
  m <- lm_robust(y ~ x, data = dat, fixed_effects = ~ block)
  nd <- dat[1, ]
  nd$block <- 999
  expect_error(predict(m, nd), "new levels")
  m2 <- lm_robust(y ~ x, data = dat, fixed_effects = ~ block + cl, se_type = "HC1")
  expect_error(predict(m2, dat[1:5, ]), "fitted.values")
})

# ---- the exported fitting function (estimatr #269) ----

test_that("#269: lm_robust_fit accepts an integer X and an unnamed y", {
  fit <- lm_robust_fit(y = dat$y, X = matrix(as.integer(dat$z)), weights = NULL,
    cluster = NULL, ci = FALSE, se_type = "none", alpha = 0.05,
    return_vcov = FALSE, try_cholesky = FALSE, has_int = TRUE)
  tidied <- estimatr:::tidy.lm_robust(fit)
  expect_s3_class(tidied, "data.frame")
  expect_equal(nrow(tidied), 1L)
})
