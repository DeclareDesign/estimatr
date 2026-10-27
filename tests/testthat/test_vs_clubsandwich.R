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

# A cluster with more rows than 2p (p without weights) has its adjustment
# computed in the span of its design columns rather than from the full
# n_c-by-n_c matrix (#276, #318). The sizes put clusters on both sides of that
# switch, and `z` is constant within cluster, so the columns spanning a cluster
# are collinear.
test_that("CR2 agrees with its definition on both sides of the low-rank switch", {
  set.seed(343)
  sizes <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 30, 200)
  d_lr <- data.frame(cl = rep(seq_along(sizes), sizes))
  n_lr <- nrow(d_lr)
  d_lr$z <- rbinom(length(sizes), 1, 0.5)[d_lr$cl]
  d_lr$x <- rnorm(n_lr)
  d_lr$w <- rexp(n_lr)
  d_lr$y <- d_lr$z + d_lr$x + rnorm(length(sizes))[d_lr$cl] + rnorm(n_lr)
  X <- model.matrix(~ z + x, d_lr)
  fit <- lm_robust(y ~ z + x, data = d_lr, clusters = cl, se_type = "CR2")
  expect_equal(unname(fit$vcov), cr2_by_definition(X, d_lr$y, d_lr$cl),
               tolerance = LIVE_TOL)
  fit <- lm_robust(y ~ z + x, data = d_lr, clusters = cl, weights = w,
                   se_type = "CR2")
  expect_equal(unname(fit$vcov), cr2_by_definition(X, d_lr$y, d_lr$cl, d_lr$w),
               tolerance = LIVE_TOL)
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

# ---- weighted 2SLS, weighted lm_lin, multivariate outcomes, hypotheses ----
#
# CR2 and its degrees of freedom on the paths no comparison above reached
# (review C4). Measured gaps are at most 1.1e-13. Every weighted cell names
# inverse_var = FALSE, the working model estimatr's weighted CR2 uses, so a
# change to clubSandwich's default cannot move what this package is held to.

test_that("weighted iv_robust CR2 matches clubSandwich on AER::ivreg", {
  skip_if_not_installed("AER")
  set.seed(42)
  n_iv <- 300
  di <- data.frame(x = rnorm(n_iv), inst = rnorm(n_iv), cl = rep(1:15, each = 20),
                   w = runif(n_iv, 0.5, 2))
  di$en <- di$inst + rnorm(n_iv, 0, 0.5)
  di$y <- di$en + rnorm(n_iv)

  fit <- iv_robust(y ~ en + x | inst + x, data = di, weights = w, clusters = cl, se_type = "CR2")
  mi <- AER::ivreg(y ~ en + x | inst + x, data = di, weights = w)
  target <- clubSandwich::vcovCR(mi, cluster = di$cl, type = "CR2", inverse_var = FALSE)
  expect_equal(unname(fit$vcov), unname(as.matrix(target)), tolerance = LIVE_TOL)
  ct <- clubSandwich::coef_test(mi, vcov = target, test = "Satterthwaite")
  expect_equal(unname(fit$df), ct$df_Satt, tolerance = LIVE_TOL)
})

test_that("weighted lm_lin CR2 matches clubSandwich on lm with the covariate centred on its weighted mean", {
  centred <- d
  centred$x_c <- d$x - weighted.mean(d$x, d$w)
  reference <- lm(y ~ z * x_c, data = centred, weights = w)
  target <- clubSandwich::vcovCR(reference, cluster = centred$cl, type = "CR2", inverse_var = FALSE)
  fit <- lm_lin(y ~ z, covariates = ~ x, data = d, weights = w, clusters = cl, se_type = "CR2")
  expect_equal(unname(fit$vcov), unname(as.matrix(target)), tolerance = LIVE_TOL)
  ct <- clubSandwich::coef_test(reference, vcov = target, test = "Satterthwaite")
  expect_equal(unname(fit$df), ct$df_Satt, tolerance = LIVE_TOL)
})

test_that("a multivariate lm_robust CR2 matches clubSandwich on the mlm fit", {
  multi <- d
  set.seed(4)
  multi$y2 <- multi$x - multi$z + rnorm(nrow(multi))
  reference <- lm(cbind(y, y2) ~ x + z, data = multi)
  target <- clubSandwich::vcovCR(reference, cluster = multi$cl, type = "CR2")
  fit <- lm_robust(cbind(y, y2) ~ x + z, data = multi, clusters = cl, se_type = "CR2")
  expect_equal(unname(fit$vcov), unname(as.matrix(target)), tolerance = LIVE_TOL)
  ct <- clubSandwich::coef_test(reference, vcov = target, test = "Satterthwaite")
  # coef_test names a multivariate fit's coefficients "<outcome>:<term>", as
  # the variance matrix does; matched by name rather than trusted by position.
  expect_equal(unname(as.vector(fit$df)),
               ct$df_Satt[match(rownames(fit$vcov), ct$Coef)], tolerance = LIVE_TOL)
})

test_that("lh_robust CR2 matches clubSandwich's contrast, degrees of freedom included", {
  # A combination of coefficients has its own Satterthwaite degrees of
  # freedom. lh_robust used to give it the smallest per-coefficient df, which
  # its comment called conservative and which is not a bound: on this data that
  # was 16.10 for the first hypothesis, where its own is 15.75.
  hypotheses <- c("z + 2*x = 0.5", "(Intercept) - z = 0")
  contrasts <- rbind(c(0, 2, 1), c(1, 0, -1))
  rhs <- c(0.5, 0)

  m <- lm(y ~ x + z, data = d)
  expected <- clubSandwich::linear_contrast(
    m, vcov = clubSandwich::vcovCR(m, cluster = d$cl, type = "CR2"),
    contrasts = contrasts, test = "Satterthwaite"
  )
  lh <- lh_robust(y ~ x + z, data = d, clusters = cl, se_type = "CR2",
                  linear_hypothesis = hypotheses)
  expect_equal(unname(lh$lh$coefficients), expected$Est - rhs, tolerance = LIVE_TOL)
  expect_equal(unname(lh$lh$std.error), expected$SE, tolerance = LIVE_TOL)
  expect_equal(unname(lh$lh$df), expected$df, tolerance = LIVE_TOL)

  # Weighted, at the working model estimatr's weighted CR2 uses.
  mw <- lm(y ~ x + z, data = d, weights = w)
  expected_w <- clubSandwich::linear_contrast(
    mw, vcov = clubSandwich::vcovCR(mw, cluster = d$cl, type = "CR2", inverse_var = FALSE),
    contrasts = contrasts, test = "Satterthwaite"
  )
  lh_w <- lh_robust(y ~ x + z, data = d, clusters = cl, weights = w, se_type = "CR2",
                    linear_hypothesis = hypotheses)
  expect_equal(unname(lh_w$lh$std.error), expected_w$SE, tolerance = LIVE_TOL)
  expect_equal(unname(lh_w$lh$df), expected_w$df, tolerance = LIVE_TOL)

  # Absorbed fixed effects, where CR2 expands the dummies behind the
  # coefficients, so the hypothesis has to be lined up with the coefficients
  # and not with the dummy columns.
  dfe <- ext_data_fe()
  mfe <- lm(y ~ x + z + factor(g), data = dfe)
  contrast_fe <- matrix(c(0, 2, 1, rep(0, nlevels(dfe$g) - 1)), nrow = 1)
  expected_fe <- clubSandwich::linear_contrast(
    mfe, vcov = clubSandwich::vcovCR(mfe, cluster = dfe$cl, type = "CR2"),
    contrasts = contrast_fe, test = "Satterthwaite"
  )
  lh_fe <- lh_robust(y ~ x + z, data = dfe, fixed_effects = ~ g, clusters = cl,
                     se_type = "CR2", linear_hypothesis = "z + 2*x = 0")
  expect_equal(unname(lh_fe$lh$std.error), expected_fe$SE, tolerance = LIVE_TOL)
  expect_equal(unname(lh_fe$lh$df), expected_fe$df, tolerance = LIVE_TOL)
})

# ---- a singleton cluster's own dummy ----
#
# A dummy for a one-member cluster is the CR2 analogue of the full-leverage
# observation that HC2 and HC3 return NA for. That cluster's (I - H) block is
# singular, the pseudo-inverse zeroes the direction, and the dummy's standard
# error is assembled from clusters that carry no information about it.
# clubSandwich takes the same pseudo-inverse and returns the same number, so
# this is the estimator as published rather than an artefact here, and it is
# pinned as such: a change that made CR2 return NA on this design, or a
# different finite number, fails this test and reopens the question on purpose.

test_that("CR2 on a singleton cluster's own dummy matches clubSandwich, and is finite", {
  set.seed(343)
  n_cl <- 30
  ds <- data.frame(cl = c(rep(seq_len(n_cl), each = 3), n_cl + 1))
  ds$x <- rnorm(nrow(ds))
  ds$d <- as.numeric(ds$cl == n_cl + 1)
  ds$y <- 0.5 * ds$x + rnorm(nrow(ds))
  fit <- lm_robust(y ~ x + d, data = ds, clusters = cl, se_type = "CR2")
  m <- lm(y ~ x + d, data = ds)
  target <- clubSandwich::vcovCR(m, cluster = ds$cl, type = "CR2")
  expect_false(anyNA(fit$std.error))
  expect_equal(unname(fit$vcov), unname(as.matrix(target)), tolerance = LIVE_TOL)
  ct <- clubSandwich::coef_test(m, vcov = target, test = "Satterthwaite")
  expect_equal(unname(fit$df), ct$df_Satt, tolerance = LIVE_TOL)
})

# ---- CR3 ----
#
# CR3 adjusts by the generalized inverse of each cluster's block of I - H where
# CR2 takes its inverse square root. For OLS and WLS it is also the sum of
# squared leave-one-cluster-out changes in the coefficients, which gives a
# reference that does not share clubSandwich's arithmetic and that still exists
# where clubSandwich's solve() fails.

loco_sum <- function(formula, data, cluster, w = rep(1, nrow(data))) {
  data$loco_w <- w
  b <- coef(lm(formula, data, weights = loco_w))
  dev <- lapply(unique(cluster), function(g) {
    b_g <- coef(lm(formula, data[cluster != g, ], weights = loco_w))
    replace(b - b_g, is.na(b_g), NA)
  })
  Reduce(`+`, lapply(dev, tcrossprod))
}

test_that("CR3 and its degrees of freedom match clubSandwich, weighted and not", {
  for (cv in c("cl", "clu")) {
    fit <- lm_robust(y ~ x + z, data = d, clusters = d[[cv]], se_type = "CR3")
    m <- lm(y ~ x + z, data = d)
    expect_equal(unname(fit$vcov),
                 unname(as.matrix(clubSandwich::vcovCR(m, cluster = d[[cv]], type = "CR3"))),
                 tolerance = LIVE_TOL, label = paste0("CR3 vcov, clusters = ", cv))
    ct <- clubSandwich::coef_test(m, vcov = "CR3", cluster = d[[cv]], test = "Satterthwaite")
    expect_equal(unname(fit$df), ct$df_Satt, tolerance = LIVE_TOL,
                 label = paste0("CR3 df, clusters = ", cv))
  }
  fit <- lm_robust(y ~ x + z, data = d, clusters = cl, weights = w, se_type = "CR3")
  m <- lm(y ~ x + z, data = d, weights = w)
  expect_equal(unname(fit$vcov),
               unname(as.matrix(clubSandwich::vcovCR(m, cluster = d$cl, type = "CR3", inverse_var = FALSE))),
               tolerance = LIVE_TOL)
  ct <- clubSandwich::coef_test(m, vcov = "CR3", cluster = d$cl, test = "Satterthwaite", inverse_var = FALSE)
  expect_equal(unname(fit$df), ct$df_Satt, tolerance = LIVE_TOL)
})

test_that("CR3 is the leave-one-cluster-out sum, under weights and absorbed fixed effects", {
  expect_equal(
    unname(lm_robust(y ~ x + z, data = d, clusters = clu, se_type = "CR3")$vcov),
    unname(loco_sum(y ~ x + z, d, d$clu)), tolerance = LIVE_TOL
  )
  expect_equal(
    unname(lm_robust(y ~ x + z, data = d, clusters = cl, weights = w, se_type = "CR3")$vcov),
    unname(loco_sum(y ~ x + z, d, d$cl, d$w)), tolerance = LIVE_TOL
  )
  d$fe <- factor(rep(1:4, length.out = nrow(d)))
  expect_equal(
    unname(lm_robust(y ~ x + z, data = d, clusters = clu, fixed_effects = ~fe, se_type = "CR3")$vcov),
    unname(loco_sum(y ~ x + z + fe, d, d$clu)[2:3, 2:3]), tolerance = LIVE_TOL
  )
})

# A one-row cluster carrying its own dummy leaves I - H singular in that block:
# clubSandwich stops with "system is computationally singular", while dropping
# the cluster still identifies the slope, and the generalized inverse agrees.
test_that("CR3 on a cluster that identifies its own dummy is the leave-one-cluster-out sum", {
  set.seed(343)
  ds <- data.frame(cl = rep(1:12, c(1, rep(5, 11))), x1 = rnorm(56))
  ds$y <- ds$x1 + rnorm(12)[ds$cl] + rnorm(56)
  fit <- lm_robust(y ~ x1 + factor(cl), data = ds, clusters = cl, se_type = "CR3")
  expect_equal(fit$vcov["x1", "x1"],
               loco_sum(y ~ x1 + factor(cl), ds, ds$cl)[2, 2],
               tolerance = LIVE_TOL)
})

test_that("iv_robust and lh_robust CR3 match clubSandwich", {
  skip_if_not_installed("AER")
  set.seed(343)
  di <- data.frame(cl = rep(1:15, each = 8), x1 = rnorm(120), zz = rbinom(120, 1, 0.5))
  di$d <- rbinom(120, 1, plogis(di$zz + di$x1))
  di$y <- di$d + di$x1 + rnorm(15)[di$cl] + rnorm(120)
  fit <- iv_robust(y ~ d + x1 | zz + x1, data = di, clusters = cl, se_type = "CR3")
  target <- clubSandwich::vcovCR(AER::ivreg(y ~ d + x1 | zz + x1, data = di), cluster = di$cl, type = "CR3")
  expect_equal(unname(fit$vcov), unname(as.matrix(target)), tolerance = LIVE_TOL)

  lh <- lh_robust(y ~ x + z, data = d, clusters = cl, se_type = "CR3", linear_hypothesis = "x + z = 0")
  lc <- clubSandwich::linear_contrast(lm(y ~ x + z, data = d), vcov = "CR3", cluster = d$cl,
                                      contrasts = matrix(c(0, 1, 1), nrow = 1), test = "Satterthwaite")
  expect_equal(lh$lh$std.error, lc$SE, tolerance = LIVE_TOL)
  expect_equal(lh$lh$df, lc$df, tolerance = LIVE_TOL)
})
