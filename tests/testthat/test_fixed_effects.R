library(estimatr)

# Fixed effects are absorbed via Frisch-Waugh-Lovell (FWL) demeaning.
# This gives bit-identical coefficients and residuals to the dummy-variable
# formulation.
#
# SE type support with fixed_effects:
#
#   Supported (exact match with estimatr):
#     HC0, HC1, stata, classical       (no clusters)
#     CR0, stata                        (with clusters)
#
#   NOT supported — will error:
#     HC2, HC3                          (require full [X|FE] hat values)
#     CR2                               (same reason)
#
# To get HC2/HC3/CR2 SEs with FE, use the dummy formulation:
#   lm_robust(y ~ x + factor(blockID), se_type = "HC2")
#
# Why: HC2, HC3, and CR2 need hat values from H = [X,Z]([X,Z]'[X,Z])^{-1}[X,Z]'
# where Z is the full FE dummy matrix.  That is O(J^2) memory and O(J^3) compute
# for J FE levels — prohibitive for individual FE.  HC0/HC1/CR0 only need
# residuals and (for CR0) cluster membership, so FWL demeaning is sufficient.

dat <- ref_data_fe()
n <- nrow(dat)

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

test_that("one-way FE default se_type is HC2 (no clusters)", {
  # Was HC1: HC2 used to be unaffordable with absorbed FE. Under a single
  # factor the leverage identity makes it exact and cheap, so the default is
  # the package's ordinary one again, and it matches estimatr for this call.
  # See test_fe_leverage.R.
  m <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl)
  expect_equal(m$se_type, "HC2")
})

test_that("FE default se_type is CR0 (with clusters)", {
  m <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl, clusters = cl)
  expect_equal(m$se_type, "CR0")
})

test_that("HC2 and HC3 now work with one-way FE and match the dummy regression", {
  # These two used to assert an error. The restriction was real for two or more
  # FE factors and wrong for one; see test_fe_leverage.R for the identity.
  for (se in c("HC2", "HC3")) {
    fe  <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl, se_type = se)
    dum <- lm_robust(y ~ z + factor(bl), data = dat, se_type = se)
    expect_equal(unname(fe$std.error), unname(dum$std.error["z"]))
  }
})

test_that("HC2 and HC3 still error with two-way FE", {
  dat2 <- dat
  dat2$bl2 <- factor(rep(1:4, length.out = nrow(dat2)))
  expect_error(
    lm_robust(y ~ z, data = dat2, fixed_effects = ~ bl + bl2, se_type = "HC2"),
    "cannot be used with `fixed_effects`"
  )
  expect_error(
    lm_robust(y ~ z, data = dat2, fixed_effects = ~ bl + bl2, se_type = "HC3"),
    "cannot be used with `fixed_effects`"
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
  m0 <- ref("fe_HC1")
  mz <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl, se_type = "HC1")
  expect_equal(coef(mz),         coef(m0),         tolerance = 1e-12)
  expect_equal(mz$std.error,     m0$std.error,     tolerance = 1e-12)
  expect_equal(mz$df.residual,   m0$df.residual)
  expect_equal(mz$fitted.values, m0$fitted.values, tolerance = 1e-10)
})

test_that("FE HC0 coefs and SEs identical to estimatr", {
  m0 <- ref("fe_HC0")
  mz <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl, se_type = "HC0")
  expect_equal(coef(mz),     coef(m0),     tolerance = 1e-12)
  expect_equal(mz$std.error, m0$std.error, tolerance = 1e-12)
})

test_that("FE stata cluster SEs identical to estimatr", {
  m0 <- ref("fe_cl_stata")
  mz <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl,
                  clusters = cl, se_type = "stata")
  expect_equal(coef(mz),     coef(m0),     tolerance = 1e-12)
  expect_equal(mz$std.error, m0$std.error, tolerance = 1e-12)
})

test_that("FE default coefs and fitted values identical to estimatr HC1", {
  m0 <- ref("fe_HC1")
  mz <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  expect_equal(coef(mz),         coef(m0),         tolerance = 1e-12)
  expect_equal(mz$fitted.values, m0$fitted.values, tolerance = 1e-10)
  expect_equal(mz$df.residual,   m0$df.residual)
})

test_that("FE weighted HC1 coefs and SEs identical to estimatr", {
  dat_w <- ref_data_fe_weighted()
  m0 <- ref("fe_w_HC1")
  mz <- lm_robust(y ~ z + x, data = dat_w, fixed_effects = ~bl,
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
