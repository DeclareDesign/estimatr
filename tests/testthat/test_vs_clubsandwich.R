library(estimatr)

# CR2 against clubSandwich.
#
# This is the highest-value comparison in the package. CR2 is the default
# cluster-robust variance here, it is the estimator this package is most often
# cited for, and clubSandwich is the only independent implementation of it in
# R. It is also intricate enough that an error would not be obvious from
# inspection: the adjustment matrices are per-cluster matrix square roots, and
# the Satterthwaite degrees of freedom are a separate calculation on top.
#
# Live rather than recorded, for the reasons given in test_vs_sandwich.R.

skip_if_not_installed("clubSandwich")

d <- ext_data_ols()

club_se <- function(fit, coefs) {
  sqrt(diag(as.matrix(fit)))[coefs]
}

# ---- CR2 variance ----

test_that("CR2 matches clubSandwich::vcovCR, balanced and unbalanced clusters", {
  m <- lm(y ~ x + z, data = d)
  for (cv in c("cl", "clu")) {
    fit <- lm_robust(y ~ x + z, data = d, clusters = d[[cv]], se_type = "CR2")
    target <- clubSandwich::vcovCR(m, cluster = d[[cv]], type = "CR2")
    expect_equal(
      unname(fit$vcov), unname(as.matrix(target)),
      tolerance = LIVE_TOL, label = paste0("CR2 vcov, clusters = ", cv)
    )
  }
})

# Under weights the two estimators follow different conventions, which
# `?lm_robust` states and which is invisible at the call site: weighted CR2 is
# built against a working model with identity covariance and weighted HC2
# against precision weights. `inverse_var` is passed explicitly on both sides
# below rather than left to clubSandwich's default, so that a change in that
# default fails the test instead of silently asserting the other convention.
test_that("weighted CR2 matches clubSandwich::vcovCR at inverse_var = FALSE", {
  m <- lm(y ~ x + z, data = d, weights = w)
  fit <- lm_robust(y ~ x + z, data = d, clusters = cl, weights = w, se_type = "CR2")
  target <- clubSandwich::vcovCR(m, cluster = d$cl, type = "CR2",
                                 inverse_var = FALSE)
  expect_equal(unname(fit$vcov), unname(as.matrix(target)), tolerance = LIVE_TOL)
})

test_that("weighted HC2 matches clubSandwich::vcovCR at inverse_var = TRUE", {
  m <- lm(y ~ x + z, data = d, weights = w)
  fit <- lm_robust(y ~ x + z, data = d, weights = w, se_type = "HC2")
  # One cluster per observation is the heteroskedasticity-consistent case.
  target <- clubSandwich::vcovCR(m, cluster = seq_len(nrow(d)), type = "CR2",
                                 inverse_var = TRUE)
  expect_equal(unname(fit$vcov), unname(as.matrix(target)), tolerance = LIVE_TOL)
})

test_that("CR0 matches clubSandwich::vcovCR type CR0", {
  m <- lm(y ~ x + z, data = d)
  fit <- lm_robust(y ~ x + z, data = d, clusters = cl, se_type = "CR0")
  target <- clubSandwich::vcovCR(m, cluster = d$cl, type = "CR0")
  expect_equal(unname(fit$vcov), unname(as.matrix(target)), tolerance = LIVE_TOL)
})

# ---- Satterthwaite degrees of freedom ----
#
# The degrees of freedom are what make a CR2 interval different from a CR0
# interval with a t quantile, and they vary by coefficient. Comparing them
# separately from the variance means a failure says which of the two is wrong.

test_that("CR2 degrees of freedom match clubSandwich Satterthwaite", {
  m <- lm(y ~ x + z, data = d)
  for (cv in c("cl", "clu")) {
    fit <- lm_robust(y ~ x + z, data = d, clusters = d[[cv]], se_type = "CR2")
    ct <- clubSandwich::coef_test(
      m,
      vcov = clubSandwich::vcovCR(m, cluster = d[[cv]], type = "CR2"),
      test = "Satterthwaite"
    )
    expect_equal(
      unname(fit$df), ct$df_Satt,
      tolerance = LIVE_TOL, label = paste0("CR2 df, clusters = ", cv)
    )
  }
})

test_that("CR2 degrees of freedom differ across coefficients", {
  # Guards the test above: if this package returned a single residual degrees
  # of freedom for every coefficient, and clubSandwich did too, the comparison
  # would pass while both were wrong. They are supposed to differ.
  fit <- lm_robust(y ~ x + z, data = d, clusters = clu, se_type = "CR2")
  expect_gt(diff(range(fit$df)), 0.1)
  expect_true(all(fit$df < length(unique(d$clu))))
})

# ---- CR2 with absorbed fixed effects ----
#
# The absorbed path computes CR2 without ever forming the dummy design.
# clubSandwich has no absorption, so the reference is the dummy expansion,
# computed here by a package that shares no code with this one. This is the
# cluster-robust counterpart of the leverage identity checked in
# test_fe_leverage.R.

test_that("CR2 with absorbed fixed effects matches clubSandwich on the dummy expansion", {
  dfe <- ext_data_fe()
  fit <- lm_robust(y ~ x + z, fixed_effects = ~ g, clusters = cl,
                   data = dfe, se_type = "CR2")
  m <- lm(y ~ x + z + factor(g), data = dfe)
  target <- clubSandwich::vcovCR(m, cluster = dfe$cl, type = "CR2")

  expect_equal(unname(fit$std.error), unname(club_se(target, c("x", "z"))),
               tolerance = LIVE_TOL)

  ct <- clubSandwich::coef_test(m, vcov = target, test = "Satterthwaite")
  expect_equal(unname(fit$df),
               ct$df_Satt[match(c("x", "z"), ct$Coef)],
               tolerance = LIVE_TOL)
})

# ---- instrumental variables ----

test_that("iv_robust CR2 matches clubSandwich on AER::ivreg", {
  skip_if_not_installed("AER")
  set.seed(42)
  n <- 300
  di <- data.frame(x = rnorm(n), inst = rnorm(n), cl = rep(1:15, each = 20))
  di$en <- di$inst + rnorm(n, 0, 0.5)
  di$y <- di$en + rnorm(n)

  fit <- iv_robust(y ~ en + x | inst + x, data = di, clusters = cl, se_type = "CR2")
  mi <- AER::ivreg(y ~ en + x | inst + x, data = di)
  target <- clubSandwich::vcovCR(mi, cluster = di$cl, type = "CR2")

  expect_equal(unname(fit$coefficients), unname(coef(mi)), tolerance = LIVE_TOL)
  expect_equal(unname(fit$std.error), unname(sqrt(diag(as.matrix(target)))),
               tolerance = LIVE_TOL)

  ct <- clubSandwich::coef_test(mi, vcov = target, test = "Satterthwaite")
  expect_equal(unname(fit$df), ct$df_Satt, tolerance = LIVE_TOL)
})

# ---- a hostile design ----
#
# ext_data_hard() in helper-external.R is built to strain the arithmetic, and
# on it clubSandwich cannot be held to LIVE_TOL; the reasons, and the tolerance
# it can be held to, are set out beside HARD_CLUB_TOL. The tight reference on
# this design is CR2 written from its definition, and clubSandwich is the loose
# one.

dh <- ext_data_hard()
hard_formula <- y ~ z + income + share + age + age2

test_that("CR2 on the hostile design matches CR2 written from its definition", {
  X <- model.matrix(hard_formula, dh)
  for (tc in c(FALSE, TRUE)) {
    fit <- lm_robust(hard_formula, data = dh, clusters = cl, se_type = "CR2",
                     try_cholesky = tc)
    expect_equal(unname(fit$vcov), cr2_by_definition(X, dh$y, dh$cl),
                 tolerance = LIVE_TOL, label = paste("try_cholesky =", tc))
    fit <- lm_robust(hard_formula, data = dh, clusters = cl, weights = w,
                     se_type = "CR2", try_cholesky = tc)
    expect_equal(unname(fit$vcov), cr2_by_definition(X, dh$y, dh$cl, dh$w),
                 tolerance = LIVE_TOL, label = paste("weighted, try_cholesky =", tc))
  }
})

test_that("CR2 and its degrees of freedom match clubSandwich on the hostile design", {
  fits <- list(
    unweighted = list(m = lm(hard_formula, data = dh), w = NULL),
    weighted = list(m = lm(hard_formula, data = dh, weights = w), w = dh$w)
  )
  for (nm in names(fits)) {
    m <- fits[[nm]]$m
    target <- clubSandwich::vcovCR(m, cluster = dh$cl, type = "CR2", inverse_var = FALSE)
    ct <- clubSandwich::coef_test(m, vcov = target, test = "Satterthwaite")
    for (tc in c(FALSE, TRUE)) {
      fit <- lm_robust(hard_formula, data = dh, clusters = cl, weights = fits[[nm]]$w,
                       se_type = "CR2", try_cholesky = tc)
      label <- paste0(nm, ", try_cholesky = ", tc)
      expect_equal(unname(fit$vcov), unname(as.matrix(target)),
                   tolerance = HARD_CLUB_TOL, label = paste(label, "vcov"))
      expect_equal(unname(fit$df), ct$df_Satt,
                   tolerance = HARD_CLUB_TOL, label = paste(label, "df"))
    }
  }

  skip_if_not_installed("AER")
  iv_formula <- y ~ en + z + income + share | inst + z + income + share
  mi <- AER::ivreg(iv_formula, data = dh)
  target <- clubSandwich::vcovCR(mi, cluster = dh$cl, type = "CR2")
  ct <- clubSandwich::coef_test(mi, vcov = target, test = "Satterthwaite")
  fit <- iv_robust(iv_formula, data = dh, clusters = cl, se_type = "CR2")
  expect_equal(unname(fit$vcov), unname(as.matrix(target)), tolerance = HARD_CLUB_TOL)
  expect_equal(unname(fit$df), ct$df_Satt, tolerance = HARD_CLUB_TOL)
})
