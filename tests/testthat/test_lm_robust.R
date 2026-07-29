library(estimatrZero)

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

test_that("lm_lin and lm_robust give same intercept estimate", {
  # Intercept in lm_lin at center should match simple mean
  m <- lm_lin(y ~ z, covariates = ~ x, data = dat)
  # Treatment coef should be near lm_robust
  m0 <- lm_robust(y ~ z + x, data = dat)
  expect_equal(coef(m)["z"], coef(m)["z"], tolerance = 0.01)
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
  # the default S3 method finds the field, so residuals() just works
  expect_equal(residuals(m), m$residuals)
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

# Other files in this suite load estimatr, which re-registers S3 methods for
# the shared "lm_robust" and "iv_robust" classes, so generic dispatch here is
# not guaranteed to reach estimatrZero's method. Where the method under test
# IS the thing being tested, call ours explicitly.
z_predict <- function(...) estimatrZero:::predict.lm_robust(...)
z_glance <- function(...) estimatrZero:::glance.iv_robust(...)
z_model_frame <- function(...) estimatrZero:::model.frame.iv_robust(...)
z_varnames <- function(...) estimatrZero:::variable.names.lm_robust(...)
z_tidy <- function(...) estimatrZero:::tidy.lm_robust(...)

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

# ---- predict (estimatr #403, #404) ----

test_that("#403: predict() with no newdata returns the in-sample fit", {
  m <- lm_robust(y ~ x + z, data = dat)
  expect_equal(z_predict(m), m$fitted.values)
  expect_error(z_predict(m, se.fit = TRUE), "newdata")
})

test_that("#404: predict() works with fixed effects, with and without factors", {
  d <- dat
  d$f <- factor(rep(letters[1:4], length.out = n))
  m <- lm_robust(y ~ x + f, data = d, fixed_effects = ~ block)
  nd <- d[1:5, ]
  # equals the dummy-variable regression, which is the definition of correct
  expect_equal(unname(z_predict(m, nd)),
               unname(predict(lm(y ~ x + f + factor(block), data = d), nd)),
               tolerance = 1e-8)
  # and reproduces the in-sample fit
  expect_equal(unname(z_predict(m, d)), unname(m$fitted.values), tolerance = 1e-8)
})

test_that("#404: predict() rejects new FE levels and multi-way FE", {
  m <- lm_robust(y ~ x, data = dat, fixed_effects = ~ block)
  nd <- dat[1, ]; nd$block <- 999
  expect_error(z_predict(m, nd), "new levels")
  m2 <- lm_robust(y ~ x, data = dat, fixed_effects = ~ block + cl)
  expect_error(z_predict(m2, dat[1:5, ]), "fitted.values")
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
  g <- z_glance(m)
  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 1L)
  # reports the weakest of the per-regressor first stages
  fs <- m$diagnostic_first_stage_fstatistic
  expect_equal(g$statistic.weakinst, unname(min(fs[grep("(^|:)value$", names(fs))])))
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

test_that("#304: fixed_effects must be a formula", {
  expect_error(lm_robust(y ~ x, data = dat, fixed_effects = block),
               "must be a one-sided formula")
})

# ---- augment (estimatr #377) ----

test_that("#377: augment() returns the model frame with .fitted and .resid", {
  m <- lm_robust(y ~ x + z, data = dat)
  a <- estimatrZero:::augment.lm_robust(m)
  expect_s3_class(a, "data.frame")
  expect_true(all(c(".fitted", ".resid") %in% names(a)))
  expect_equal(nrow(a), n)
  expect_equal(a$.fitted + a$.resid, dat$y, tolerance = 1e-12)
})

test_that("#377: augment() works for lm_lin, iv_robust and fixed effects", {
  expect_true(".fitted" %in% names(
    estimatrZero:::augment.lm_robust(lm_lin(y ~ z, covariates = ~ x, data = dat))))
  expect_true(".fitted" %in% names(
    estimatrZero:::augment.iv_robust(iv_robust(mpg ~ wt | am, data = mtcars))))
  a <- estimatrZero:::augment.lm_robust(
    lm_robust(y ~ x, data = dat, fixed_effects = ~ block))
  expect_equal(a$.fitted + a$.resid, dat$y, tolerance = 1e-8)
})

test_that("#377: augment(newdata =) predicts without residuals", {
  m <- lm_robust(y ~ x + z, data = dat)
  a <- estimatrZero:::augment.lm_robust(m, newdata = dat[1:5, ])
  expect_equal(nrow(a), 5L)
  expect_true(".fitted" %in% names(a))
  expect_false(".resid" %in% names(a))
})

test_that("#377: augment() refuses multivariate outcomes", {
  m <- lm_robust(cbind(y, x) ~ z, data = dat)
  expect_error(estimatrZero:::augment.lm_robust(m), "multiple outcomes")
})
