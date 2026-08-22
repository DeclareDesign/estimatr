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

test_that("weighted CR2 matches clubSandwich::vcovCR", {
  m <- lm(y ~ x + z, data = d, weights = w)
  fit <- lm_robust(y ~ x + z, data = d, clusters = cl, weights = w, se_type = "CR2")
  target <- clubSandwich::vcovCR(m, cluster = d$cl, type = "CR2")
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
