library(estimatr)

# lm_lin's own behaviour: how it names and centres the design it builds, its
# defaults, predict() on the rebuilt design, and the regressions filed against
# it.
#
# lm_lin against Lin's specification written out by hand, across treatment and
# covariate kinds, weighted and clustered, is test_lm_lin_equivalence.R.

lin_data <- function() {
  set.seed(42)
  n <- 200
  d <- data.frame(y = rnorm(n), z = rbinom(n, 1, 0.5), x = rnorm(n), cl = rep(1:20, 10))
  d
}

d <- lin_data()

test_that("the design is the treatment, the centred covariates, and their interactions", {
  m <- lm_lin(y ~ z, covariates = ~ x, data = d)
  expect_equal(names(coef(m)), c("(Intercept)", "z", "x_c", "z:x_c"))
  expect_equal(m$se_type, "HC2")
  expect_equal(lm_lin(y ~ z, covariates = ~ x, data = d, clusters = cl)$se_type, "CR2")
})

test_that("covariates are centred on their mean, weighted when there are weights", {
  expect_equal(unname(lm_lin(y ~ z, covariates = ~ x, data = d)$scaled_center), mean(d$x),
               tolerance = 1e-10)

  set.seed(12)
  n <- 500
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  ps <- plogis(0.3 * x1 + 0.3 * x2)
  z <- rbinom(n, 1, ps)
  obs <- data.frame(y = 2 * z + 0.5 * x1 + 0.3 * x2 + rnorm(n), z = z, x1 = x1, x2 = x2,
                    w = ifelse(z == 1, 1 / ps, 1 / (1 - ps)))
  m <- lm_lin(y ~ z, covariates = ~ x1 + x2, data = obs, weights = w)
  expect_equal(unname(m$scaled_center), c(weighted.mean(x1, obs$w), weighted.mean(x2, obs$w)),
               tolerance = 1e-10)
})

test_that("a weighted R-squared with informative weights stays lm()'s", {
  set.seed(13)
  n <- 200
  x <- 1:n / n
  obs <- data.frame(y = x + rnorm(n, 0, 0.1), x = x, z = rbinom(n, 1, 0.5), w = x + 0.1)
  m <- lm_lin(y ~ z, covariates = ~ x, data = obs, weights = w)
  obs$x_c <- obs$x - weighted.mean(obs$x, obs$w)
  expect_equal(m$r.squared, summary(lm(y ~ z * x_c, data = obs, weights = w))$r.squared,
               tolerance = 1e-12)
})

test_that("#345: lm_lin returns residuals", {
  m <- lm_lin(y ~ z, covariates = ~ x, data = d)
  expect_equal(length(m$residuals), nrow(d))
  expect_equal(unname(m$residuals + m$fitted.values), d$y, tolerance = 1e-12)
})

test_that("A4: lm_lin without an intercept expands a 0/1 treatment", {
  # Without an intercept there is no baseline to absorb the control group, so
  # both indicators are needed. 2.0 expanded only when the treatment took
  # values outside {0, 1}, and returned `z, x1_c, z:x1_c`, losing the
  # control-group intercept; 1.0.6 returned all four terms.
  set.seed(1)
  N <- 100
  dd <- data.frame(z = rbinom(N, 1, 0.5), x1 = rnorm(N))
  dd$y <- dd$z + dd$x1 + rnorm(N)

  m <- lm_lin(y ~ z - 1, covariates = ~ x1, data = dd)
  expect_equal(names(coef(m)), c("z0", "z1", "z0:x1_c", "z1:x1_c"))
  # equal to the same model fitted by hand, which is what the terms mean
  dd$x1_c <- dd$x1 - mean(dd$x1)
  hand <- lm_robust(y ~ factor(z):x1_c + factor(z) - 1, data = dd, se_type = "HC2")
  expect_equal(unname(coef(m)), unname(coef(hand)), tolerance = 1e-10)
  expect_equal(length(predict(m, newdata = dd)), N)

  # the intercept form is unchanged
  expect_equal(names(coef(lm_lin(y ~ z, covariates = ~ x1, data = dd))),
               c("(Intercept)", "z", "x1_c", "z:x1_c"))
})

# ---- predict ----

# predict() rebuilt the lm_lin design by looking the treatment up by term
# label, so a factor `z` produced X[, "z"] when the design matrix holds `zb`
# and `zc`: subscript out of bounds for every treatment that is not a single
# binary column. The design is now rebuilt the way lm_lin() builds it and lined
# up against the coefficients by name.

lin_cases <- ref_lin_cases()
lin_dat <- ref_data_lin()

for (case in lin_cases) {
  test_that(paste0("predict on an lm_lin fit works: ", case$lbl), {
    fit <- lm_lin(reformulate(case$z, "y"), covariates = case$cov, data = lin_dat)
    p <- predict(fit, newdata = lin_dat)
    expect_length(p, nrow(lin_dat))
    expect_false(anyNA(p))
    # Predicting on the fitting data must reproduce the in-sample fit.
    expect_equal(unname(p), unname(fit$fitted.values))
  })
}

for (case in lin_cases) {
  test_that(paste0("predict agrees with estimatr: ", case$lbl), {
    f <- reformulate(case$z, "y")
    fit_z <- lm_lin(f, covariates = case$cov, data = lin_dat)
    p_z <- predict(fit_z, newdata = lin_dat)
    expect_equal(unname(p_z), ref(paste0("post_predict_", case$lbl)))
  })
}

test_that("predict on a subset of newdata uses the levels seen at fit time", {
  sub <- lin_dat[lin_dat$zf != "c", ]
  fit_z <- lm_lin(y ~ zf, covariates = ~ x, data = lin_dat)
  expect_equal(unname(predict(fit_z, newdata = sub)), ref("post_predict_subset"))
})

test_that("a numeric multi-valued treatment keeps its fit-time levels", {
  fit <- lm_lin(y ~ zn, covariates = ~ x, data = lin_dat)
  # Stored so predict() expands against these rather than against whatever
  # values happen to appear in newdata.
  expect_equal(unname(fit$treatment_vals), c(2, 5))
  expect_null(lm_lin(y ~ zf, covariates = ~ x, data = lin_dat)$treatment_vals)
})

test_that("newdata carrying only one treatment level still predicts", {
  # The stored terms keep the fit's factor levels, so the design rebuilds with
  # the unused indicator columns at zero rather than losing them. Worth pinning:
  # this is the case where a design rebuilt from `newdata` alone would silently
  # come out narrower than the coefficient vector.
  nd <- lin_dat
  nd$zf <- factor(rep("a", nrow(nd)), levels = "a")
  fit_z <- lm_lin(y ~ zf, covariates = ~ x, data = lin_dat)
  expect_equal(unname(predict(fit_z, newdata = nd)), ref("post_predict_one_level"))
})
