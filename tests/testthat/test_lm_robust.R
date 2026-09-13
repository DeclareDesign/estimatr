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

# ---- collinearity (estimatr #351, #411) ----

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

test_that("#351: a constant regressor is detected as collinear, as in lm()", {
  d <- data.frame(y = rnorm(500), x = 1)
  m <- suppressWarnings(lm_robust(y ~ x, data = d))
  expect_true(is.na(m$coefficients[["x"]]))
  expect_equal(unname(is.na(coef(lm(y ~ x, data = d)))), unname(is.na(m$coefficients)))
})

# ---- leverage at or near 1 (estimatr #395) ----

test_that("#395: NaN standard errors from leverage-1 points are explained", {
  set.seed(7)
  N <- 50
  d <- data.frame(x = sample(1:40, N, TRUE), Z = sample(0:1, N, TRUE))
  d$Y <- 0.1 * d$Z + d$x + rnorm(N)
  # This design is rank deficient AND near-saturated, so several warnings fire
  # together: collinearity, leverage, and on some platforms a negative variance
  # diagonal. Collect them instead of nesting expect_warning(), which pins the
  # count as well as the content and so breaks whenever another one is added --
  # as it did on all five CI platforms and on none locally.
  ws <- character(0)
  withCallingHandlers(
    lm_lin(Y ~ Z, covariates = ~ as.factor(x), data = d),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("leverage", ws)))
  expect_true(any(grepl("collinear", ws)))
  expect_false(any(grepl("NaNs produced", ws, fixed = TRUE)))
  # classical SEs do not use leverage, so they stay finite
  m <- suppressWarnings(lm_lin(Y ~ Z, covariates = ~ as.factor(x), data = d,
                               se_type = "classical"))
  expect_false(is.nan(m$std.error[["Z"]]))
})

# The tests below pin what actually goes wrong at a high-leverage point,
# because the #395 test above pins only that a warning fires and the NEWS entry
# for it described a different failure than the one that occurs.
#
# A hat value is a projection diagonal and cannot exceed 1. Leverage of exactly
# 1 costs nothing in the answer: that observation has residual exactly 0, its
# meat contribution is a 0/0 the code resolves to 0, and HC2 returns a finite
# standard error (first test). The trouble is leverage computed marginally ABOVE
# 1, where 1 - h is a small negative number. HC2 then divides by it, half_meat
# takes the square root of the negative result, and every standard error in the
# fit is NaN however small the offending term. HC3 squares the denominator,
# which cancels the sign, so it returned a finite number carrying a spurious
# positive term and said nothing -- the worse of the two failures, because it
# looks like an answer.
#
# lm_variance() now sets the denominator to 0 wherever 1 - h <= 0, which sends
# both through the same isfinite trap that leverage of exactly 1 already used,
# and counts how many observations sit at or near leverage 1 so lm_robust() can
# warn on the condition rather than on a NaN. The count is tolerant, at
# sandwich's `h > 1 - sqrt(eps)`, while the clamp is the strict `1 - h <= 0`,
# because the two answer different questions. The clamp decides the number, and
# on an exactly saturated design the number is the same whichever side of 1 the
# rounding lands on, so a strict count would leave the warning to an ulp. The
# first test below pins the design where that is so.

test_that("#395: leverage of exactly 1 answers finitely, and says that it did", {
  # The single "c" observation is alone in its cell, so it is fitted exactly:
  # leverage is exactly 1 and the residual is exactly 0, and its contribution to
  # the meat is a 0/0 that is correctly resolved to 0 rather than to NaN. The
  # standard error is therefore built from 8 rows where the fit used 9, which is
  # what the warning reports. sandwich warns on this design too and on the same
  # grounds, returning NaN where estimatr drops the row and answers.
  d <- data.frame(
    g = factor(c("a", "a", "a", "a", "b", "b", "b", "b", "c")),
    Z = c(0, 1, 0, 1, 0, 1, 0, 1, 1),
    Y = c(1.0, 2.0, 1.5, 2.5, 3.0, 4.0, 3.5, 4.5, 9.0)
  )
  X <- model.matrix(Y ~ Z + g, d)
  qrx <- qr(X)
  h <- rowSums(qr.Q(qrx)^2)
  e <- as.vector(d$Y - X %*% qr.coef(qrx, d$Y))
  expect_equal(max(h), 1)          # exactly, not approximately
  expect_equal(min(1 - h), 0)
  expect_equal(e[9], 0)

  # The fit's own solver puts that same hat value at 1 + 2.2e-16 rather than at
  # the exact 1 qr() returns above, so this design is the one that decides
  # whether the warning is tolerant or is left to an ulp. Both halves are
  # pinned: the warning fires and names one observation, and the standard error
  # does not depend on which side of 1 the rounding landed.
  expect_warning(
    m <- lm_robust(Y ~ Z + g, data = d, se_type = "HC2"),
    "1 observation has a computed leverage at or near 1"
  )
  expect_false(is.nan(m$std.error[["Z"]]))
  expect_true(is.finite(m$std.error[["Z"]]))

  # Computed here rather than recorded. With observation 9 contributing 0 the
  # meat is the other 8 rows, and the HC2 standard error follows in closed form
  # from a bread and a hat vector that never go near estimatr's solver.
  M <- summary(lm(Y ~ Z + g, data = d))$cov.unscaled
  denom <- 1 - rowSums((X %*% M) * X)
  omega <- ifelse(denom <= 0, 0, e^2 / denom)
  se_hand <- sqrt(diag(M %*% (t(X) %*% (X * omega)) %*% M))
  expect_equal(m$std.error[["Z"]], se_hand[["Z"]])
})

test_that("#395: the leverage guard drops exactly the 1 - h < 0 observations", {
  # Whether a real fit puts a hat value above 1 is a rounding accident and
  # differs by platform, so the guard itself is pinned at the level of the
  # function that implements it, where the leverage is chosen rather than
  # observed. XtX_inv = 0.3 makes h = 0.3 * x^2, so the fourth observation sits
  # at 1.2 and the first three at 0.3.
  X <- matrix(c(1, 1, 1, 2), ncol = 1)
  M <- matrix(0.3, 1, 1)
  ei <- matrix(rep(1, 4), ncol = 1)
  h <- 0.3 * as.vector(X)^2
  expect_equal(h, c(0.3, 0.3, 0.3, 1.2))

  z_variance <- function(se_type) {
    estimatr:::lm_variance(
      X = X, Xunweighted = NULL, XtX_inv = M, ei = ei, weight_mean = 1,
      cluster = NULL, J = 0L, ci = TRUE, se_type = se_type,
      which_covs = TRUE, fe_rank = 0L, fe_leverage = NULL, n_eff = -1L
    )
  }

  # The guarded answers are the three well-behaved rows and nothing else.
  hc2 <- z_variance("HC2")
  hc3 <- z_variance("HC3")
  expect_equal(hc2[["n_leverage_near_one"]], 1L)
  expect_equal(hc3[["n_leverage_near_one"]], 1L)
  expect_equal(sqrt(hc2[["Vcov_hat"]][1, 1]), sqrt(0.09 * 3 / 0.7))
  expect_equal(sqrt(hc3[["Vcov_hat"]][1, 1]), sqrt(0.09 * 3 / 0.49))

  # What the same inputs gave before the guard, computed here rather than
  # recorded: HC2 was NaN throughout, and HC3 was finite, silent, and four times
  # too large, the offending row supplying 94% of the variance.
  expect_true(suppressWarnings(is.nan(sqrt(sum(0.09 * as.vector(X)^2 / (1 - h))))))
  expect_equal(sqrt(sum(0.09 * as.vector(X)^2 / (1 - h)^2)), 3.0904725, tolerance = 1e-7)
  expect_gt(3.0904725 / sqrt(hc3[["Vcov_hat"]][1, 1]), 4)

  # se_types that never read leverage are untouched and report no count.
  for (ty in c("HC0", "HC1", "classical")) {
    expect_equal(z_variance(ty)[["n_leverage_near_one"]], 0L,
                 label = paste(ty, "leverage count"))
  }
})

# These designs are rank deficient, so a collinearity warning fires on every fit
# regardless of se_type. Only the leverage warning is of interest here, so fits
# are run through this rather than through expect_silent().
fit_warnings <- function(expr) {
  ws <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = val, leverage_warning = any(grepl("leverage", ws, fixed = TRUE)))
}

# Leverage above 1 is a floating-point outcome, so it is asserted rather than
# assumed: if a platform's solver returns max(h) <= 1 for this design, the
# condition under test is absent and the guard has nothing to do.
lev_above_one <- function(fml, data) {
  X <- model.matrix(fml, data)
  fit <- estimatr:::lm_solver(X, matrix(data[[all.vars(fml)[1]]], ncol = 1), TRUE)
  keep <- !is.na(fit$beta_hat)
  Xk <- X[, keep, drop = FALSE]
  h <- rowSums((Xk %*% fit$XtX_inv) * Xk)
  list(any = any(1 - h < 0), max = max(h), h = h)
}

test_that("#395: HC2 and HC3 both warn, and neither returns NaN, above leverage 1", {
  set.seed(7)
  N <- 50
  d <- data.frame(x = sample(1:40, N, TRUE), Z = sample(0:1, N, TRUE))
  d$Y <- 0.1 * d$Z + d$x + rnorm(N)
  fml <- Y ~ Z * as.factor(x)

  lev <- lev_above_one(fml, d)
  skip_if_not(lev$any, "this platform's solver does not put any leverage above 1")
  expect_gt(lev$max, 1)

  # Before the guard HC2 was NaN here and HC3 was finite and silent. The
  # asymmetry between them was the defect; both now warn and both answer.
  for (ty in c("HC2", "HC3")) {
    m <- fit_warnings(lm_robust(fml, data = d, se_type = ty))
    expect_true(m$leverage_warning, label = paste(ty, "leverage warning"))
    expect_false(is.nan(m$fit$std.error[["Z"]]), label = paste(ty, "std.error is NaN"))
  }

  # The se_types that never touch leverage neither warn nor change.
  for (ty in c("HC0", "HC1", "classical")) {
    mm <- fit_warnings(lm_robust(fml, data = d, se_type = ty))
    expect_false(mm$leverage_warning, label = paste(ty, "leverage warning"))
  }
})

test_that("#395: the leverage warning counts the observations it dropped", {
  set.seed(7)
  N <- 50
  d <- data.frame(x = sample(1:40, N, TRUE), Z = sample(0:1, N, TRUE))
  d$Y <- 0.1 * d$Z + d$x + rnorm(N)
  fml <- Y ~ Z * as.factor(x)
  skip_if_not(lev_above_one(fml, d)$any, "no leverage above 1 on this platform")

  # Nesting expect_warning() would pin the number of warnings as well as their
  # content, which is what broke the first #395 test on all five CI platforms
  # and on none locally. Collect them and read the leverage one out.
  ws <- character(0)
  withCallingHandlers(
    lm_robust(fml, data = d, se_type = "HC2"),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("collinear", ws)))
  lev <- grep("computed leverage at or near 1", ws, value = TRUE)
  expect_length(lev, 1)

  # The number in the message is the count the C++ made, and the reference is
  # the same criterion applied to a hat vector computed outside it. This design
  # has more than one such observation, so the plural branch of the message is
  # exercised here and the singular branch in the exactly saturated test above.
  expect_match(lev, "^[0-9]+ observations have ")
  expect_equal(as.integer(sub(" .*$", "", lev)),
               sum(lev_above_one(fml, d)$h > 1 - sqrt(.Machine$double.eps)))
})

test_that("#395: CR2 has no analogous hole -- degeneracy goes through its clamp", {
  # CR2 never forms 1 - h_ii. It eigendecomposes
  # (I - H) - H' + Xo MUWTWUM Xo' per cluster and keeps 1/sqrt(eigenvalue) only
  # above 1e-12, so a degenerate cluster contributes 0 rather than a negative or
  # a sign-flipped term. Cluster 11 below is a singleton carrying its own factor
  # level, which is fitted exactly.
  set.seed(11)
  d <- data.frame(g = factor(c(rep("a", 10), rep("b", 10), "c")),
                  Z = c(rep(0:1, 5), rep(0:1, 5), 1))
  d$Y <- rnorm(21) + as.numeric(d$g)
  d$cl <- c(rep(1:5, 2), rep(6:10, 2), 11)

  X <- model.matrix(Y ~ Z + g, d)
  M <- solve(crossprod(X))
  MUWTWUM <- M %*% crossprod(X) %*% M
  block_eigen <- function(cl) {
    ix <- which(d$cl == cl)
    Xb <- X[ix, , drop = FALSE]
    H <- Xb %*% M %*% t(Xb)
    A <- (diag(length(ix)) - H) - t(H) + Xb %*% MUWTWUM %*% t(Xb)
    min(eigen(A, symmetric = TRUE)$values)
  }
  expect_lte(block_eigen(11), 1e-12)     # the clamp fires
  expect_gt(min(vapply(1:10, block_eigen, numeric(1))), 1e-12)

  m <- lm_robust(Y ~ Z + g, data = d, clusters = cl, se_type = "CR2")
  expect_true(is.finite(m$std.error[["Z"]]))
  expect_gt(m$std.error[["Z"]], 0)
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
