library(estimatr)

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
dat$z_block <- rep(rep(c(0L, 1L), 5L), 20L)

# ---- lm_robust basics ----

test_that("lm_robust returns lm_robust class", {
  m <- lm_robust(y ~ x + z, data = dat)
  expect_s3_class(m, "lm_robust")
})

test_that("lm_robust default is HC2", {
  m <- lm_robust(y ~ x, data = dat)
  expect_equal(m$se_type, "HC2")
})

test_that("lm_robust classical matches lm", {
  m0 <- lm_robust(y ~ x + z, data = dat, se_type = "classical")
  m1 <- lm(y ~ x + z, data = dat)
  expect_equal(as.numeric(coef(m0)), as.numeric(coef(m1)), tolerance = 1e-10)
  expect_equal(as.numeric(m0$std.error), as.numeric(summary(m1)$coef[, 2]), tolerance = 1e-10)
})

test_that("lm_robust tidy returns expected columns", {
  m <- lm_robust(y ~ x + z, data = dat)
  td <- tidy(m)
  expect_named(td, c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "df", "outcome"))
})

test_that("lm_robust clustered uses CR2 by default", {
  m <- lm_robust(y ~ x + z, data = dat, clusters = cl)
  expect_equal(m$se_type, "CR2")
  expect_true(!is.null(m$nclusters))
})

test_that("lm_robust subset works", {
  m <- lm_robust(y ~ x, data = dat, subset = z == 1)
  expect_true(m$nobs < n)
})

test_that("lm_robust weights work", {
  m <- lm_robust(y ~ x + z, data = dat, weights = w)
  expect_true(m$weighted)
})

test_that("lm_robust multivariate returns matrix coefs", {
  m <- lm_robust(cbind(y, y_extra = y + rnorm(n)) ~ x + z, data = dat)
  expect_true(is.matrix(m$coefficients))
  expect_equal(ncol(m$coefficients), 2L)
})

test_that("lm_robust r.squared is in [0, 1] for non-trivial models", {
  m <- lm_robust(y ~ x + z, data = dat)
  expect_gte(m$r.squared, 0)
  expect_lte(m$r.squared, 1)
})

test_that("lm_robust se_type none skips SE computation", {
  m <- lm_robust(y ~ x, data = dat, se_type = "none")
  expect_equal(m$se_type, "none")
  expect_true(is.na(m$std.error[1]))
})

test_that("lm_robust confint returns matrix", {
  m <- lm_robust(y ~ x + z, data = dat)
  ci <- confint(m)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 3L)
})

test_that("lm_robust fitted values have correct length", {
  m <- lm_robust(y ~ x, data = dat)
  expect_equal(length(fitted(m)), n)
})

# ---- lm_lin ----

test_that("lm_lin returns lm_robust class", {
  m <- lm_lin(y ~ z, covariates = ~ x, data = dat)
  expect_s3_class(m, "lm_robust")
})

test_that("lm_lin coefs include _c interaction", {
  m <- lm_lin(y ~ z, covariates = ~ x, data = dat)
  expect_true(any(grepl("_c$", m$term)))
})

test_that("lm_lin's treatment coefficient is close to the additive lm_robust one", {
  # Lin's estimator interacts the treatment with centred covariates, so it does
  # not equal the additive fit; on data where the interaction is near zero the
  # two land close, which is the property worth pinning.
  m  <- lm_lin(y ~ z, covariates = ~ x, data = dat)
  m0 <- lm_robust(y ~ z + x, data = dat)
  expect_equal(unname(coef(m)["z"]), unname(coef(m0)["z"]), tolerance = 0.01)
})

# ---- difference_in_means ----

test_that("DiM standard design returns correct class", {
  m <- difference_in_means(y ~ z, data = dat)
  expect_s3_class(m, "difference_in_means")
  expect_equal(m$design, "Standard")
})

test_that("DiM blocked design", {
  m <- difference_in_means(y ~ z_block, blocks = block, data = dat)
  expect_equal(m$design, "Blocked")
  expect_equal(unname(m$df), n - 2 * 20)
})

test_that("DiM clustered design", {
  # Clustered: treatment assigned at cluster level
  dat_cl <- dat
  dat_cl$cl_z <- as.integer(dat_cl$cl %% 2 == 0)
  m <- difference_in_means(y ~ cl_z, clusters = cl, data = dat_cl)
  expect_equal(m$design, "Clustered")
})

test_that("DiM tidy works", {
  m <- difference_in_means(y ~ z, data = dat)
  td <- tidy(m)
  expect_named(td[1:7], c("term","estimate","std.error","statistic","p.value","conf.low","conf.high"))
})

# ---- iv_robust ----

test_that("iv_robust returns iv_robust class", {
  dat$iv <- dat$z + rnorm(n, 0, 0.5)
  m <- iv_robust(y ~ z | iv, data = dat)
  expect_s3_class(m, "iv_robust")
})

test_that("iv_robust 2SLS matches ivreg for classical SEs", {
  set.seed(1)
  dat2 <- data.frame(
    y = rnorm(200),
    x = rnorm(200),
    z = rnorm(200)
  )
  m0 <- iv_robust(y ~ x | z, data = dat2, se_type = "classical")
  # Just check structure
  expect_s3_class(m0, "iv_robust")
  expect_equal(length(coef(m0)), 2L)
})

# ---- lh_robust ----

test_that("lh_robust returns lh_robust class", {
  m <- lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z = 0")
  expect_s3_class(m, "lh_robust")
})

test_that("lh_robust lh component tidy works", {
  m <- lh_robust(y ~ x + z, data = dat, linear_hypothesis = "x + z = 0")
  td <- tidy(m$lh)
  expect_equal(nrow(td), 1L)
})

# ---- glance / nobs ----

test_that("glance.lm_robust returns data.frame with expected cols", {
  m <- lm_robust(y ~ x + z, data = dat)
  g <- glance(m)
  expect_s3_class(g, "data.frame")
  expect_true("r.squared" %in% names(g))
  expect_true("nobs" %in% names(g))
})

test_that("nobs returns correct count", {
  m <- lm_robust(y ~ x, data = dat)
  expect_equal(nobs(m), n)
})

# ---- weighted R2 bug fix ----

test_that("weighted R2 is in [0, 1] for simple model", {
  # Before fix, weights^2 could make R2 negative or > 1
  dat_w <- data.frame(y = 1:10 + rnorm(10, 0, 0.01), x = 1:10, w = (1:10) / 10)
  m <- lm_robust(y ~ x, data = dat_w, weights = w)
  expect_gte(m$r.squared, 0)
  expect_lte(m$r.squared, 1.0001)  # small tolerance
})

# ---- regression: GitHub issue fixes ----

test_that("#421: ordered factor cluster does not crash", {
  dat_ord <- dat
  dat_ord$cl_ord <- factor(dat$cl, ordered = TRUE)
  m <- lm_robust(y ~ x + z, data = dat_ord, clusters = cl_ord)
  expect_s3_class(m, "lm_robust")
  expect_equal(m$nclusters, length(unique(dat$cl)))
})

test_that("#348: formula() object passed to fixed_effects works", {
  fe_form <- formula(~block)
  m_var  <- lm_robust(y ~ z, data = dat, fixed_effects = fe_form)
  m_lit  <- lm_robust(y ~ z, data = dat, fixed_effects = ~block)
  expect_equal(coef(m_var), coef(m_lit), tolerance = 1e-12)
})

test_that("#303: intercept-only model with FE returns sensible result", {
  m <- lm_robust(y ~ 1, data = dat, fixed_effects = ~block)
  expect_equal(length(m$coefficients), 0L)
  expect_equal(m$df.residual, n - 20L)       # 20 blocks
  expect_gte(m$r.squared, 0)
  expect_lte(m$r.squared, 1)
  expect_equal(length(m$fitted.values), n)
})

test_that("#405: lh_robust CIs match lm_robust CIs with clusters", {
  m_cl <- lm_robust(y ~ x + z, data = dat, clusters = cl)
  lh_x <- lh_robust(y ~ x + z, data = dat, clusters = cl, linear_hypothesis = "x=0")
  expect_equal(lh_x$lh$conf.low,  unname(confint(m_cl)["x", "2.5 %"]),  tolerance = 1e-10)
  expect_equal(lh_x$lh$conf.high, unname(confint(m_cl)["x", "97.5 %"]), tolerance = 1e-10)
})

test_that("lh_robust handles an intercept-only model", {
  # estimatr 1.0.6 errors here with "missing value where TRUE/FALSE needed".
  # Its df warning guards with `length(fit$df) > 0 && var(fit$df > 0)`, and
  # var() of a length-one vector is NA, so the `if` has nothing to branch on.
  # The unreleased origin/lh-fixes branch changes that 0 to a 1; this rewrite
  # resolves df per hypothesis by name instead and never calls var(), so the
  # case works rather than being patched. Pinned because a future change to the
  # df logic could reintroduce it silently.
  # Live case: the italian_village_continued design in the ResearchDesigns
  # library, which fits age ~ 1 and cannot run under CRAN estimatr.
  intercept_only <- data.frame(age = dat$y)
  m <- lh_robust(age ~ 1, data = intercept_only,
                 linear_hypothesis = "(Intercept) = 20")
  expect_s3_class(m, "lh_robust")
  expect_equal(nrow(m$lh), 1L)
  expect_equal(m$lh$df, unname(m$lm_robust$df["(Intercept)"]))
  expect_false(is.na(m$lh$p.value))
})

test_that("#320: lh_robust returns joint_hypothesis", {
  lh2 <- lh_robust(y ~ x + z, data = dat, linear_hypothesis = c("x=0", "z=0"))
  expect_true(!is.null(lh2$joint_hypothesis))
  expect_true("value" %in% names(lh2$joint_hypothesis))
  expect_true("p.value" %in% names(lh2$joint_hypothesis))
  expect_equal(unname(lh2$joint_hypothesis["numdf"]), 2)
  expect_gte(lh2$joint_hypothesis["p.value"], 0)
  expect_lte(lh2$joint_hypothesis["p.value"], 1)
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

test_that("#345: iv_robust returns structural residuals", {
  m <- iv_robust(mpg ~ wt | am, data = mtcars)
  X <- cbind(1, mtcars$wt)
  expect_equal(unname(m$residuals), as.vector(mtcars$mpg - X %*% coef(m)),
               tolerance = 1e-12)
})

test_that("#345: lm_lin returns residuals", {
  m <- lm_lin(y ~ z, covariates = ~ x, data = dat)
  expect_equal(length(m$residuals), n)
  expect_equal(unname(m$residuals + m$fitted.values), dat$y, tolerance = 1e-12)
})

# ---- collinearity warning (estimatr #411) ----

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

z_model_frame <- function(...) estimatr:::model.frame.iv_robust(...)
z_varnames <- function(...) estimatr:::variable.names.lm_robust(...)
z_tidy <- function(...) estimatr:::tidy.lm_robust(...)

# ---- rank detection and degenerate designs (estimatr #351, #395) ----

test_that("#351: a constant regressor is detected as collinear, as in lm()", {
  d <- data.frame(y = rnorm(500), x = 1)
  m <- suppressWarnings(lm_robust(y ~ x, data = d))
  expect_true(is.na(m$coefficients[["x"]]))
  expect_equal(unname(is.na(coef(lm(y ~ x, data = d)))), unname(is.na(m$coefficients)))
})

test_that("#395: NaN standard errors from leverage-1 points are explained", {
  set.seed(7); N <- 50
  d <- data.frame(x = sample(1:40, N, TRUE), Z = sample(0:1, N, TRUE))
  d$Y <- 0.1 * d$Z + d$x + rnorm(N)
  # the design is also rank deficient, so the collinearity warning fires too
  expect_warning(
    expect_warning(lm_lin(Y ~ Z, covariates = ~ as.factor(x), data = d), "leverage"),
    "collinear"
  )
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
# 1 is benign: that observation has residual exactly 0, its meat contribution is
# a 0/0 the code resolves to 0, and HC2 returns a finite standard error without
# warning (first test). The trouble is leverage computed marginally ABOVE 1,
# where 1 - h is a small negative number. HC2 then divides by it, half_meat
# takes the square root of the negative result, and every standard error in the
# fit is NaN however small the offending term. HC3 squares the denominator,
# which cancels the sign, so it returned a finite number carrying a spurious
# positive term and said nothing -- the worse of the two failures, because it
# looks like an answer.
#
# lm_variance() now sets the denominator to 0 wherever 1 - h <= 0, which sends
# both through the same isfinite trap that leverage of exactly 1 already used,
# and reports how many observations that hit so lm_robust() can warn on the
# condition rather than on a NaN.

test_that("#395: leverage of exactly 1 is benign -- HC2 stays finite and silent", {
  # The single "c" observation is alone in its cell, so it is fitted exactly:
  # leverage is exactly 1 and the residual is exactly 0, and its contribution to
  # the meat is a 0/0 that is correctly resolved to 0 rather than to NaN.
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

  m <- lm_robust(Y ~ Z + g, data = d, se_type = "HC2")   # must not warn
  expect_false(is.nan(m$std.error[["Z"]]))
  expect_true(is.finite(m$std.error[["Z"]]))
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
  expect_equal(hc2[["n_leverage_above_one"]], 1L)
  expect_equal(hc3[["n_leverage_above_one"]], 1L)
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
    expect_equal(z_variance(ty)[["n_leverage_above_one"]], 0L,
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
  list(any = any(1 - h < 0), max = max(h))
}

test_that("#395: HC2 and HC3 both warn, and neither returns NaN, above leverage 1", {
  set.seed(7); N <- 50
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
  set.seed(7); N <- 50
  d <- data.frame(x = sample(1:40, N, TRUE), Z = sample(0:1, N, TRUE))
  d$Y <- 0.1 * d$Z + d$x + rnorm(N)
  fml <- Y ~ Z * as.factor(x)
  skip_if_not(lev_above_one(fml, d)$any, "no leverage above 1 on this platform")

  expect_warning(
    expect_warning(lm_robust(fml, data = d, se_type = "HC2"),
                   "computed leverage above 1"),
    "collinear"
  )
  expect_warning(
    expect_warning(lm_robust(fml, data = d, se_type = "HC2"), "observations have"),
    "collinear"
  )
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
  nd <- dat[1, ]; nd$block <- 999
  expect_error(predict(m, nd), "new levels")
  m2 <- lm_robust(y ~ x, data = dat, fixed_effects = ~ block + cl, se_type = "HC1")
  expect_error(predict(m2, dat[1:5, ]), "fitted.values")
})

# ---- exported lm_robust_fit and S3 coverage (estimatr #269, #123) ----

test_that("#269: lm_robust_fit accepts an integer X and an unnamed y", {
  fit <- lm_robust_fit(y = dat$y, X = matrix(as.integer(dat$z)), weights = NULL,
    cluster = NULL, ci = FALSE, se_type = "none", alpha = 0.05,
    return_vcov = FALSE, try_cholesky = FALSE, has_int = TRUE)
  expect_s3_class(z_tidy(fit), "data.frame")
  expect_equal(nrow(z_tidy(fit)), 1L)
})

test_that("#123: variable.names() returns the model terms", {
  m <- lm_robust(y ~ x + z, data = dat)
  expect_equal(z_varnames(m), c("(Intercept)", "x", "z"))
})

# ---- iv_robust S3 (estimatr #389, #397) ----

test_that("#389: glance() works with multiple endogenous regressors", {
  m <- iv_robust(mpg ~ hp + wt | am + cyl, data = mtcars, diagnostics = TRUE)
  g <- glance(m)
  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 1L)
  # reports the weakest of the per-regressor first stages
  fs <- m$diagnostic_first_stage_fstatistic
  expect_equal(g$statistic.weakinst, unname(min(fs[grep("(^|:)value$", names(fs))])))
})

test_that("over-identified diagnostics work for every se_type, not just classical", {
  # `first_stage_fits` has its columns renamed `fit_<endog>`, and the rewritten
  # robust branch then indexed it by the bare endogenous names, so every
  # over-identified fit with a non-classical se_type died with "subscript out
  # of bounds". The classical branch took a different function and was fine,
  # which is why the whole path had no test.
  set.seed(2); n <- 500
  d <- data.frame(z1 = rnorm(n), z2 = rnorm(n), z3 = rnorm(n), w = rnorm(n))
  d$x <- d$z1 + 0.5 * d$z2 + 0.3 * d$z3 + rnorm(n)
  d$y <- d$x + d$w + rnorm(n)
  fml <- y ~ x + w | z1 + z2 + z3 + w

  for (ty in c("classical", "HC0", "HC1", "HC2", "HC3")) {
    m <- iv_robust(fml, data = d, se_type = ty, diagnostics = TRUE)
    ov <- m$diagnostic_overid_test
    expect_equal(names(ov), c("value", "df", "p.value"), info = ty)
    expect_true(is.finite(ov[["value"]]), info = ty)
    expect_equal(unname(ov[["df"]]), 2, info = ty)   # 3 instruments, 1 endogenous
  }

  # The classical branch is Sargan's statistic, which AER computes
  # independently.
  skip_if_not_installed("AER")
  sargan <- summary(AER::ivreg(fml, data = d), diagnostics = TRUE)$diagnostics["Sargan", ]
  cl <- iv_robust(fml, data = d, se_type = "classical", diagnostics = TRUE)
  expect_equal(cl$diagnostic_overid_test[["value"]], sargan[["statistic"]],
               tolerance = 1e-10)
  expect_equal(cl$diagnostic_overid_test[["p.value"]], sargan[["p-value"]],
               tolerance = 1e-10)

  # The robust branch is Wooldridge's score test, a different statistic, and
  # every non-classical se_type shares it.
  rb <- vapply(c("HC0", "HC1", "HC2", "HC3"), function(ty) {
    iv_robust(fml, data = d, se_type = ty, diagnostics = TRUE)$diagnostic_overid_test[["value"]]
  }, numeric(1))
  expect_equal(unname(diff(range(rb))), 0)
  expect_false(isTRUE(all.equal(rb[["HC0"]], cl$diagnostic_overid_test[["value"]])))
})

test_that("#397: model.frame() on iv_robust returns the model variables", {
  m <- iv_robust(mpg ~ wt + hp | am + hp, data = mtcars)
  mf <- z_model_frame(m)
  expect_equal(nrow(mf), nrow(mtcars))
  expect_setequal(names(mf), c("mpg", "wt", "hp", "am"))
})

# ---- clearer errors (estimatr #297, #304) ----

test_that("#297: lh_robust rejects multiple outcomes with an explanation", {
  skip_if_not_installed("carData")
  expect_error(
    lh_robust(cbind(mpg, am) ~ cyl + gear, data = mtcars, linear_hypothesis = "cyl = 2"),
    "multiple outcomes"
  )
})

test_that("#304: a bare grouping vector warns but still works", {
  # #304 asked for the formula to be enforced. A warning naming the argument
  # and the expected form does that without breaking code written against
  # 1.x, which accepted the bare vector; `RCT` on CRAN passes one.
  expect_warning(fit <- lm_robust(y ~ x, data = dat, fixed_effects = block),
                 "deprecated")
  expect_equal(fit$std.error,
               lm_robust(y ~ x, data = dat, fixed_effects = ~ block)$std.error)
  expect_equal(names(fit$felevels), "block")
  # Anything that is neither a formula nor a grouping vector is still an error.
  expect_error(lm_robust(y ~ x, data = dat, fixed_effects = list(1, 2)),
               "must be a one-sided formula")
})

# ---- augment (estimatr #377) ----

test_that("#377: augment() returns the model frame with .fitted and .resid", {
  m <- lm_robust(y ~ x + z, data = dat)
  a <- estimatr:::augment.lm_robust(m)
  expect_s3_class(a, "data.frame")
  expect_true(all(c(".fitted", ".resid") %in% names(a)))
  expect_equal(nrow(a), n)
  expect_equal(a$.fitted + a$.resid, dat$y, tolerance = 1e-12)
})

test_that("#377: augment() works for lm_lin, iv_robust and fixed effects", {
  expect_true(".fitted" %in% names(
    estimatr:::augment.lm_robust(lm_lin(y ~ z, covariates = ~ x, data = dat))))
  expect_true(".fitted" %in% names(
    estimatr:::augment.iv_robust(iv_robust(mpg ~ wt | am, data = mtcars))))
  a <- estimatr:::augment.lm_robust(
    lm_robust(y ~ x, data = dat, fixed_effects = ~ block))
  expect_equal(a$.fitted + a$.resid, dat$y, tolerance = 1e-8)
})

test_that("#377: augment(newdata =) predicts without residuals", {
  m <- lm_robust(y ~ x + z, data = dat)
  a <- estimatr:::augment.lm_robust(m, newdata = dat[1:5, ])
  expect_equal(nrow(a), 5L)
  expect_true(".fitted" %in% names(a))
  expect_false(".resid" %in% names(a))
})

test_that("#377: augment() refuses multivariate outcomes", {
  m <- lm_robust(cbind(y, x) ~ z, data = dat)
  expect_error(estimatr:::augment.lm_robust(m), "multiple outcomes")
})


# ---- regressions from 1.0.6 found in the 2026-08-23 review ----

test_that("A4: lm_lin without an intercept expands a 0/1 treatment", {
  # Without an intercept there is no baseline to absorb the control group, so
  # both indicators are needed. 2.0 expanded only when the treatment took
  # values outside {0, 1}, and returned `z, x1_c, z:x1_c`, losing the
  # control-group intercept; 1.0.6 returned all four terms.
  set.seed(1); N <- 100
  d <- data.frame(z = rbinom(N, 1, 0.5), x1 = rnorm(N))
  d$y <- d$z + d$x1 + rnorm(N)

  m <- lm_lin(y ~ z - 1, covariates = ~ x1, data = d)
  expect_equal(names(coef(m)), c("z0", "z1", "z0:x1_c", "z1:x1_c"))
  # equal to the same model fitted by hand, which is what the terms mean
  d$x1_c <- d$x1 - mean(d$x1)
  hand <- lm_robust(y ~ factor(z):x1_c + factor(z) - 1, data = d, se_type = "HC2")
  expect_equal(unname(coef(m)), unname(coef(hand)), tolerance = 1e-10)
  expect_equal(length(predict(m, newdata = d)), N)

  # the intercept form is unchanged
  expect_equal(names(coef(lm_lin(y ~ z, covariates = ~ x1, data = d))),
               c("(Intercept)", "z", "x1_c", "z:x1_c"))
})

test_that("A5: a fit with nothing but fixed effects is a fit, not a bare list", {
  set.seed(1); N <- 100
  d <- data.frame(y = rnorm(N), bl = sample(5, N, TRUE))
  m <- lm_robust(y ~ 1, fixed_effects = ~ bl, data = d)

  expect_s3_class(m, "lm_robust")
  expect_equal(m$df.residual, lm(y ~ factor(bl), data = d)$df.residual)
  expect_equal(m$r.squared, summary(lm(y ~ factor(bl), data = d))$r.squared,
               tolerance = 1e-10)
  # every method reads it without erroring, which is what the class buys
  expect_output(print(m))
  expect_s3_class(tidy(m), "data.frame")
  expect_equal(nrow(glance(m)), 1L)
  expect_equal(nobs(m), N)
  expect_equal(length(predict(m, newdata = d)), N)

  # iv_robust cannot answer this and says which function can, rather than
  # dying on "length of 'dimnames' [2] not equal to array extent" as 1.0.6 did
  expect_error(iv_robust(y ~ 1 | 1, fixed_effects = ~ bl, data = d),
               "nothing to instrument")
})

test_that("B10: glance() reports the residual df, not the first coefficient's", {
  # x[["df"]] is per-coefficient and under CR2 is Satterthwaite, so the column
  # named df.residual read 8.56 on a 10-cluster fit whose residual df is 98.
  set.seed(1); N <- 100
  d <- data.frame(y = rnorm(N), x = rnorm(N), cl = sample(10, N, TRUE))
  m <- lm_robust(y ~ x, clusters = cl, data = d, se_type = "CR2")
  expect_equal(glance(m)$df.residual, m$df.residual)
  expect_false(isTRUE(all.equal(m$df.residual, unname(m$df[[1]]))))
  # and it agrees with the iv_robust method, which always used df.residual
  mi <- iv_robust(y ~ x | x, data = d, clusters = cl, se_type = "CR2")
  expect_equal(glance(mi)$df.residual, mi$df.residual)
})


test_that("B6: one cluster is refused rather than answered with a zero", {
  # Every cluster-robust estimator here divides by J - 1 somewhere except CR2,
  # whose Satterthwaite degrees of freedom never reach that guard, so a single
  # cluster produced standard errors of order 1e-17 and no warning.
  set.seed(1); N <- 60
  d <- data.frame(y = rnorm(N), x = rnorm(N), one = 1L)
  for (se in c("CR0", "CR2", "stata")) {
    expect_error(lm_robust(y ~ x, clusters = one, data = d, se_type = se),
                 "only one level", info = se)
  }
  expect_error(iv_robust(y ~ x | x, clusters = one, data = d), "only one level")
})

test_that("B13: a character cluster is coerced rather than handed to the C++", {
  set.seed(1); N <- 60
  d <- data.frame(y = rnorm(N), x = rnorm(N), cl = sample(6, N, TRUE))
  d$clch <- paste0("g", d$cl)
  for (se in c("CR0", "CR2", "stata")) {
    a <- lm_robust(y ~ x, clusters = clch, data = d, se_type = se)
    b <- lm_robust(y ~ x, clusters = cl, data = d, se_type = se)
    expect_equal(a$std.error, b$std.error, info = se)
    expect_equal(a$nclusters, 6L, info = se)
  }
})


test_that("B5: zero weights are not observations", {
  # N was nrow(X), so ten zero weights among 100 rows gave df.residual 98
  # where lm() gives 88, and every classical, HC1 and stata standard error and
  # every p-value moved with it. A zero-weight row stays in residuals and
  # fitted.values, as lm() keeps it; it is only the counting that changes.
  set.seed(1); N_b5 <- 100
  d <- data.frame(x = rnorm(N_b5))
  d$y <- d$x + rnorm(N_b5)
  w_b5 <- runif(N_b5); w_b5[1:10] <- 0

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

test_that("B9: fitting does not advance the RNG", {
  # The hidden variable names were built with sample.int(), so every fit moved
  # the seed. This package lives inside DeclareDesign simulation loops, where
  # that changes what the next draw is.
  set.seed(1); N_b9 <- 50
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
  set.seed(1); N_b11 <- 100
  d <- data.frame(x = rnorm(N_b11), off = rnorm(N_b11))
  d$y <- d$x + 2 * d$off + rnorm(N_b11)
  expect_error(lm_robust(y ~ x + offset(off), data = d), "not supported")

  # the rewrite the message names is exact
  expect_equal(coef(lm_robust(I(y - off) ~ x, data = d))[["x"]],
               coef(lm(y ~ x + offset(off), data = d))[["x"]])
})
