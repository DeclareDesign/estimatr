library(estimatrZero)

skip_if_not_installed("estimatr")

set.seed(42)
n  <- 200
dat <- data.frame(
  y  = rnorm(n),
  y2 = rnorm(n),
  x  = rnorm(n),
  z  = rbinom(n, 1, 0.5),
  bl = rep(1:20, each = 10),
  cl = rep(1:10, 20),
  iv = rnorm(n) + rnorm(n, 0, 0.3)
)

# ---- FWL equivalence ----

test_that("FE demeaning gives same coefs as dummy regression", {
  m_fe    <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  m_dummy <- lm_robust(y ~ z + x + factor(bl), data = dat)
  expect_equal(unname(coef(m_fe)["z"]), unname(coef(m_dummy)["z"]), tolerance = 1e-10)
  expect_equal(unname(coef(m_fe)["x"]), unname(coef(m_dummy)["x"]), tolerance = 1e-10)
})

test_that("FE residuals equal dummy regression residuals", {
  m_fe    <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  m_dummy <- lm_robust(y ~ z + x + factor(bl), data = dat)
  resid_fe    <- dat$y - m_fe$fitted.values
  resid_dummy <- dat$y - m_dummy$fitted.values
  expect_equal(resid_fe, resid_dummy, tolerance = 1e-10)
})

test_that("FE df.residual is n - k - (B-1) - 1", {
  m <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  # 2 covariates, 20 blocks (19 dummies + intercept absorbed) → 200 - 2 - 20 = 178
  expect_equal(m$df.residual, n - 2L - 20L)
})

# ---- defaults and se_type behaviour ----

test_that("FE default se_type is HC1/stata (no clusters)", {
  m <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl)
  expect_equal(m$se_type, "HC1")
})

test_that("FE default se_type is CR0 (with clusters)", {
  m <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl, clusters = cl)
  expect_equal(m$se_type, "CR0")
})

test_that("HC2 with FE errors with informative message", {
  expect_error(
    lm_robust(y ~ z, data = dat, fixed_effects = ~bl, se_type = "HC2"),
    "HC2.*fixed_effects|fixed_effects.*HC2"
  )
})

test_that("HC3 with FE errors", {
  expect_error(
    lm_robust(y ~ z, data = dat, fixed_effects = ~bl, se_type = "HC3"),
    "HC3"
  )
})

test_that("CR2 with FE errors", {
  expect_error(
    lm_robust(y ~ z, data = dat, fixed_effects = ~bl, clusters = cl, se_type = "CR2"),
    "CR2"
  )
})

test_that("HC1 with FE does not warn", {
  expect_no_warning(lm_robust(y ~ z, data = dat, fixed_effects = ~bl, se_type = "HC1"))
})

test_that("classical with FE does not warn", {
  expect_no_warning(lm_robust(y ~ z, data = dat, fixed_effects = ~bl, se_type = "classical"))
})

# ---- numerical identity with estimatr ----
#
# HC0/HC1/stata/classical/CR0 match estimatr exactly via FWL: these SE types
# use only residuals (and cluster memberships for CR0), not hat values, so the
# FWL demeaned-X model gives the same result as the full [X|FE dummies] model.
#
# HC2, HC3, CR2 are deliberately not supported with fixed_effects: they require
# hat values from the full [X|FE] design matrix which is O(J^2) memory and
# computation.  Users who need those SEs should use the dummy formulation:
#   lm_robust(y ~ x + factor(fe_var), se_type = "HC2")

test_that("FE HC1 coefs, SEs, df, and fitted values identical to estimatr", {
  m0 <- estimatr::lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl, se_type = "HC1")
  mz <- estimatrZero::lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl, se_type = "HC1")
  expect_equal(coef(mz),         coef(m0),         tolerance = 1e-12)
  expect_equal(mz$std.error,     m0$std.error,     tolerance = 1e-12)
  expect_equal(mz$df.residual,   m0$df.residual)
  expect_equal(mz$fitted.values, m0$fitted.values, tolerance = 1e-10)
})

test_that("FE HC0 coefs and SEs identical to estimatr", {
  m0 <- estimatr::lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl, se_type = "HC0")
  mz <- estimatrZero::lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl, se_type = "HC0")
  expect_equal(coef(mz),     coef(m0),     tolerance = 1e-12)
  expect_equal(mz$std.error, m0$std.error, tolerance = 1e-12)
})

test_that("FE stata cluster SEs identical to estimatr", {
  m0 <- estimatr::lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl,
                              clusters = cl, se_type = "stata")
  mz <- estimatrZero::lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl,
                                  clusters = cl, se_type = "stata")
  expect_equal(coef(mz),     coef(m0),     tolerance = 1e-12)
  expect_equal(mz$std.error, m0$std.error, tolerance = 1e-12)
})

test_that("FE default coefs and fitted values identical to estimatr HC1", {
  m0 <- estimatr::lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl, se_type = "HC1")
  mz <- estimatrZero::lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  expect_equal(coef(mz),         coef(m0),         tolerance = 1e-12)
  expect_equal(mz$fitted.values, m0$fitted.values, tolerance = 1e-10)
  expect_equal(mz$df.residual,   m0$df.residual)
})

test_that("FE weighted HC1 coefs and SEs identical to estimatr", {
  set.seed(99)
  dat_w <- dat
  dat_w$w <- runif(n, 0.5, 2)
  m0 <- estimatr::lm_robust(y ~ z + x, data = dat_w, fixed_effects = ~bl,
                              weights = w, se_type = "HC1")
  mz <- estimatrZero::lm_robust(y ~ z + x, data = dat_w, fixed_effects = ~bl,
                                  weights = w, se_type = "HC1")
  expect_equal(coef(mz),     coef(m0),     tolerance = 1e-10)
  expect_equal(mz$std.error, m0$std.error, tolerance = 1e-10)
})

# ---- return object ----

test_that("FE model has proj_r.squared and r.squared", {
  m <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  expect_true(!is.null(m$proj_r.squared))
  expect_true(!is.null(m$r.squared))
  expect_gte(m$r.squared, 0)
  expect_lte(m$r.squared, 1)
})

test_that("FE fes flag is TRUE", {
  m <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl)
  expect_true(m$fes)
})

test_that("no FE fes flag is FALSE", {
  m <- lm_robust(y ~ z, data = dat)
  expect_false(m$fes)
})

test_that("tidy works on FE model", {
  m  <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  td <- tidy(m)
  expect_s3_class(td, "data.frame")
  expect_true("z" %in% td$term)
})

test_that("summary works on FE model", {
  m <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  expect_no_error(summary(m))
})

# ---- multi-way FE ----

test_that("two-way FE converges and gives sensible results", {
  # Two-way FE: block + cluster
  m_2way <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl + cl)
  m_dum  <- lm_robust(y ~ z + factor(bl) + factor(cl), data = dat)
  expect_equal(unname(coef(m_2way)["z"]), unname(coef(m_dum)["z"]), tolerance = 1e-7)
})

# ---- iv_robust with FE ----

test_that("iv_robust with FE returns correct class", {
  m <- iv_robust(y ~ z | iv, data = dat, fixed_effects = ~bl)
  expect_s3_class(m, "iv_robust")
  expect_true(m$fes)
})

test_that("iv_robust with FE coefs match FWL manually", {
  # Demean y, z (endogenous), iv manually then run 2SLS
  y_dm  <- dat$y  - ave(dat$y,  dat$bl, FUN = mean)
  z_dm  <- dat$z  - ave(dat$z,  dat$bl, FUN = mean)
  iv_dm <- dat$iv - ave(dat$iv, dat$bl, FUN = mean)
  dat_dm <- data.frame(y=y_dm, z=z_dm, iv=iv_dm)

  m_fe  <- iv_robust(y ~ z | iv, data = dat, fixed_effects = ~bl)
  m_man <- iv_robust(y ~ z | iv, data = dat_dm)
  expect_equal(unname(coef(m_fe)["z"]), unname(coef(m_man)["z"]), tolerance = 1e-10)
})

test_that("iv_robust FE diagnostics are suppressed with warning", {
  expect_warning(
    iv_robust(y ~ z | iv, data = dat, fixed_effects = ~bl, diagnostics = TRUE),
    "Diagnostics"
  )
})
