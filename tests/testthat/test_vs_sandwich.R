library(estimatr)

# Heteroskedasticity-consistent variance against sandwich.
#
# sandwich is the reference implementation of HC0 through HC3 in R and shares
# no code with this package. These comparisons are live rather than recorded:
# both sides are computed in the same session on the same BLAS, so they can be
# held to a far tighter tolerance than a recording allows, and sandwich is
# stable enough that a live call is not a source of spurious failure.
#
# The surface is every se_type this package offers that sandwich also computes,
# crossed with weighted and unweighted, because weighting is where the
# definition of the hat value becomes contestable and where this package parts
# company with Stata (see test_vs_stata.R).

skip_if_not_installed("sandwich")

d <- ext_data_ols()

# `vcov` rather than `std.error`: the full matrix is what downstream inference
# uses, and comparing only the diagonal would miss an error in the covariances
# that lh_robust and the multi-parameter tests depend on.
expect_vcov_equal <- function(fit, target, label) {
  expect_equal(
    unname(fit$vcov), unname(as.matrix(target)),
    tolerance = LIVE_TOL, label = label
  )
}

# ---- HC0 through HC3, no weights ----

test_that("lm_robust HC0-HC3 match sandwich::vcovHC", {
  m <- lm(y ~ x + z, data = d)
  for (ty in c("HC0", "HC1", "HC2", "HC3")) {
    fit <- lm_robust(y ~ x + z, data = d, se_type = ty)
    expect_vcov_equal(fit, sandwich::vcovHC(m, type = ty), ty)
  }
})

test_that("lm_robust classical matches stats::vcov", {
  m <- lm(y ~ x + z, data = d)
  fit <- lm_robust(y ~ x + z, data = d, se_type = "classical")
  expect_vcov_equal(fit, vcov(m), "classical")
})

# ---- HC0 through HC3, with weights ----
#
# This package computes hat values from the weighted design, which is what
# sandwich and stats do. Stata answers differently, and that difference is
# pinned in test_vs_stata.R rather than left as a gap.

test_that("weighted lm_robust HC0-HC3 match sandwich::vcovHC", {
  m <- lm(y ~ x + z, data = d, weights = w)
  for (ty in c("HC0", "HC1", "HC2", "HC3")) {
    fit <- lm_robust(y ~ x + z, data = d, weights = w, se_type = ty)
    expect_vcov_equal(fit, sandwich::vcovHC(m, type = ty), paste0(ty, ", weighted"))
  }
})

test_that("weighted lm_robust classical matches stats::vcov", {
  m <- lm(y ~ x + z, data = d, weights = w)
  fit <- lm_robust(y ~ x + z, data = d, weights = w, se_type = "classical")
  expect_vcov_equal(fit, vcov(m), "classical, weighted")
})

# ---- cluster-robust ----
#
# sandwich spells the two conventions this package calls CR0 and "stata" as
# arguments to vcovCL: `cadjust` controls the G/(G-1) factor and `type`
# controls the residual adjustment. Naming both explicitly here is the point of
# the test, since the difference between them is a factor most users never see
# and would not detect if it were wrong.

test_that("CR0 matches sandwich::vcovCL with no cluster adjustment", {
  m <- lm(y ~ x + z, data = d)
  for (cv in c("cl", "clu")) {
    fit <- lm_robust(y ~ x + z, data = d, clusters = d[[cv]], se_type = "CR0")
    expect_vcov_equal(
      fit,
      sandwich::vcovCL(m, cluster = d[[cv]], type = "HC0", cadjust = FALSE),
      paste0("CR0, clusters = ", cv)
    )
  }
})

test_that("se_type 'stata' matches sandwich::vcovCL with the Stata corrections", {
  m <- lm(y ~ x + z, data = d)
  for (cv in c("cl", "clu")) {
    fit <- lm_robust(y ~ x + z, data = d, clusters = d[[cv]], se_type = "stata")
    expect_vcov_equal(
      fit,
      sandwich::vcovCL(m, cluster = d[[cv]], type = "HC1", cadjust = TRUE),
      paste0("stata, clusters = ", cv)
    )
  }
})

test_that("weighted cluster-robust matches sandwich::vcovCL", {
  m <- lm(y ~ x + z, data = d, weights = w)
  fit_cr0 <- lm_robust(y ~ x + z, data = d, clusters = cl, weights = w, se_type = "CR0")
  expect_vcov_equal(
    fit_cr0,
    sandwich::vcovCL(m, cluster = d$cl, type = "HC0", cadjust = FALSE),
    "CR0, weighted"
  )
  fit_st <- lm_robust(y ~ x + z, data = d, clusters = cl, weights = w, se_type = "stata")
  expect_vcov_equal(
    fit_st,
    sandwich::vcovCL(m, cluster = d$cl, type = "HC1", cadjust = TRUE),
    "stata, weighted"
  )
})

# ---- the corrections are not interchangeable ----
#
# Every assertion above would still pass if this package and sandwich both
# ignored `cadjust`. This one fails in that case, and so establishes that the
# tests above are discriminating between the conventions rather than agreeing
# on a constant.

test_that("the CR0 and Stata cluster corrections differ as expected", {
  m <- lm(y ~ x + z, data = d)
  cr0 <- lm_robust(y ~ x + z, data = d, clusters = cl, se_type = "CR0")
  sta <- lm_robust(y ~ x + z, data = d, clusters = cl, se_type = "stata")

  # G/(G-1) * (N-1)/(N-K), applied to the variance.
  g <- length(unique(d$cl))
  n <- nrow(d)
  k <- length(coef(m))
  expect_equal(
    unname(diag(sta$vcov)),
    unname(diag(cr0$vcov)) * (g / (g - 1)) * ((n - 1) / (n - k)),
    tolerance = LIVE_TOL
  )
  expect_gt(max(abs(diag(sta$vcov) - diag(cr0$vcov))), 0)
})

# ---- instrumental variables ----

iv_test_data <- function() {
  set.seed(42)
  n <- 300
  di <- data.frame(x = rnorm(n), inst = rnorm(n), cl = rep(1:15, each = 20))
  di$en <- di$inst + rnorm(n, 0, 0.5)
  di$y <- di$en + rnorm(n)
  di
}

test_that("iv_robust HC0, HC1 and CR0 match sandwich on AER::ivreg", {
  skip_if_not_installed("AER")
  di <- iv_test_data()
  mi <- AER::ivreg(y ~ en + x | inst + x, data = di)

  for (ty in c("HC0", "HC1")) {
    fit <- iv_robust(y ~ en + x | inst + x, data = di, se_type = ty)
    expect_equal(unname(fit$vcov), unname(sandwich::vcovHC(mi, type = ty)),
                 tolerance = LIVE_TOL, label = paste0("iv ", ty))
  }
  fit <- iv_robust(y ~ en + x | inst + x, data = di, clusters = cl, se_type = "CR0")
  expect_equal(
    unname(fit$vcov),
    unname(sandwich::vcovCL(mi, cluster = di$cl, type = "HC0", cadjust = FALSE)),
    tolerance = LIVE_TOL
  )
})

# Two-stage least squares admits two different leverage values, and which one
# you get depends on which 2SLS implementation sandwich is handed. Writing
# X-hat for the instrumented design, this package uses the leverage of the
# second-stage regression,
#
#     h_i = xhat_i' (Xhat'Xhat)^-1 xhat_i,
#
# the diagonal of an orthogonal projection. The other candidate is the
# derivative of the fitted value with respect to y_i,
#
#     h_i = xhat_i' (Xhat'Xhat)^-1 x_i,
#
# the diagonal of H* = X(X'ZZ'Z^-1 Z'X)^-1 X'Z(Z'Z)^-1 Z', the matrix taking y
# to the fitted values. Both sum to K.
#
# This is not an open question in the literature, and this package is on the
# settled side of it. Belsley, Kuh and Welsch (1980) treat 2SLS diagnostics,
# observe that H* is idempotent but asymmetric, and recommend the second-stage
# hatvalues precisely because the diagonal of an asymmetric matrix is not a
# leverage. Fox, Kleiber and Zeileis implement that recommendation as the
# default in the ivreg package, whose vignette states that the diagonal
# elements of H* "can't be treated as summary measures of leverage, that is,
# as hatvalues". Zeileis is also the author of sandwich.
#
# So the divergence is not with sandwich. It is with AER::ivreg, whose
# hatvalues method predates ivreg and returns diag(H*): the very quantity
# ivreg's authors deprecate. sandwich itself has no leverage convention of its
# own; it calls hatvalues() on whatever fit it is given. Handed a fit from the
# maintained implementation, it agrees with this package exactly, which is what
# the first test below asserts.
#
# The difference against AER's convention is not small, and an earlier version
# of this comment said it was. It reaches 8.6% at HC2 and 18.5% at HC3 on mtcars,
# where both are defined; see the vignette. Constructing both by hand is what
# makes the second test a test rather than a note.

test_that("iv_robust HC2 and HC3 match sandwich on an ivreg::ivreg fit", {
  skip_if_not_installed("ivreg")
  di <- iv_test_data()
  mi <- ivreg::ivreg(y ~ en + x | inst + x, data = di)

  for (type in c("HC2", "HC3")) {
    fit <- iv_robust(y ~ en + x | inst + x, data = di, se_type = type)
    expect_equal(
      unname(fit$std.error),
      unname(sqrt(diag(sandwich::vcovHC(mi, type = type)))),
      tolerance = LIVE_TOL
    )
  }

  # And the agreement is at the level of the leverage itself, not a
  # coincidence of the sandwich filling.
  X <- model.matrix(~ en + x, di)
  Z <- model.matrix(~ inst + x, di)
  xhat <- Z %*% solve(crossprod(Z), crossprod(Z, X))
  h_second_stage <- rowSums((xhat %*% solve(crossprod(xhat))) * xhat)
  expect_equal(unname(hatvalues(mi)), unname(h_second_stage), tolerance = LIVE_TOL)
})

test_that("iv_robust uses second-stage leverage, not the influence diagonal", {
  di <- iv_test_data()

  X <- model.matrix(~ en + x, di)
  Z <- model.matrix(~ inst + x, di)
  xhat <- Z %*% solve(crossprod(Z), crossprod(Z, X))
  bread <- solve(crossprod(xhat))
  resid <- as.vector(di$y - X %*% (bread %*% crossprod(xhat, di$y)))

  h_second_stage <- rowSums((xhat %*% bread) * xhat)
  h_influence <- rowSums((xhat %*% bread) * X)
  expect_equal(sum(h_second_stage), ncol(X))
  expect_equal(sum(h_influence), ncol(X))

  meat_vcov <- function(h, power) {
    adj <- resid / (1 - h)^power
    bread %*% crossprod(xhat * adj) %*% bread
  }

  for (spec in list(c("HC2", "0.5"), c("HC3", "1"))) {
    ty <- spec[1]
    power <- as.numeric(spec[2])
    fit <- iv_robust(y ~ en + x | inst + x, data = di, se_type = ty)

    expect_equal(unname(fit$vcov), unname(meat_vcov(h_second_stage, power)),
                 tolerance = LIVE_TOL, label = paste0("iv ", ty, " second-stage leverage"))
    # The two are genuinely different, bounded from both sides so that a change
    # in the size of the gap is a failure rather than a silent drift. On this
    # data the gap is small; on mtcars it reaches 8.6% at HC2 and 18.5% at HC3.
    rel <- max(abs(sqrt(diag(meat_vcov(h_second_stage, power))) -
                     sqrt(diag(meat_vcov(h_influence, power)))) /
                 sqrt(diag(fit$vcov)))
    expect_gt(rel, 1e-6)
    expect_lt(rel, 1e-2)
  }
})

# Why the 2SLS leverage choice is unforced.
#
# HC2 exists to make the variance estimator unbiased under homoskedasticity.
# For least squares that works exactly, because the hat matrix is a symmetric
# projection and E[e_i^2] = sigma^2 (1 - h_ii). For 2SLS the corresponding
# matrix H = X (Xhat'Xhat)^-1 Xhat' is idempotent but NOT symmetric, so
#
#     E[e_i^2] / sigma^2 = [(I - H)(I - H)']_ii = 1 - 2 h_ii + [H H']_ii,
#
# which is not 1 - h_ii for either candidate leverage, and which can exceed 1.
# So neither this package's convention nor sandwich's is the unbiased one:
# both are approximations to a quantity neither computes, and the HC2
# derivation simply does not extend uniquely to 2SLS.
#
# A Monte Carlo under homoskedasticity puts the two within a few tenths of a
# percent of each other in bias, on either side depending on the coefficient,
# against the five to eight percent that HC2 is there to remove in the first
# place. There is therefore no third implementation that could settle the
# question, because it is not a question about correctness.
#
# This test exists so that the asymmetry is on the record, and a later change
# toward sandwich's convention is made deliberately rather than on the
# assumption that sandwich is canonical.

test_that("the 2SLS hat matrix is idempotent but not symmetric", {
  di <- iv_test_data()
  X <- model.matrix(~ en + x, di)
  Z <- model.matrix(~ inst + x, di)
  xhat <- Z %*% solve(crossprod(Z), crossprod(Z, X))
  bread <- solve(crossprod(xhat))
  hat <- X %*% bread %*% t(xhat)

  expect_equal(hat %*% hat, hat, tolerance = 1e-10)
  expect_gt(max(abs(hat - t(hat))), 0.01)

  # The exact homoskedastic expectation, and the two approximations to it.
  resid_maker <- diag(nrow(X)) - hat
  exact <- diag(tcrossprod(resid_maker))
  h_second_stage <- rowSums((xhat %*% bread) * xhat)
  h_influence <- rowSums((xhat %*% bread) * X)

  expect_equal(exact, 1 - 2 * h_influence + diag(tcrossprod(hat)), tolerance = 1e-10)

  # Neither convention reproduces it, and the exact value is not even of the
  # form 1 - h for any leverage bounded by one.
  expect_gt(max(abs(exact - (1 - h_second_stage))), 1e-4)
  expect_gt(max(abs(exact - (1 - h_influence))), 1e-4)
  expect_gt(max(exact), 1)
})

# Why this package keeps its convention.
#
# Second-stage leverage is the diagonal of an orthogonal projection onto the
# column space of Xhat, so it lies in [0, 1] by construction and 1 - h is never
# negative. diag(H*) carries no such bound, because H* is idempotent but not
# symmetric: on mtcars it is already negative for one observation, and when the
# first stage is weak it routinely exceeds 1, at which point sqrt(1 - h) is
# undefined and the HC2 standard errors are NaN. sandwich says as much itself,
# warning that "HC2 covariances are numerically unstable for hat values close
# to 1".
#
# That boundedness argument is Belsley, Kuh and Welsch's, not this package's,
# and it is why the maintained implementation adopted it. The test below is
# kept because it fixes the regime where the two conventions visibly separate.

test_that("second-stage leverage stays in [0, 1] where influence leverage does not", {
  # A deliberately weak first stage: the instrument explains almost none of the
  # endogenous regressor. This is the regime the two conventions separate in.
  set.seed(20260822)
  n <- 40
  di <- data.frame(x = rnorm(n), inst = rnorm(n))
  di$en <- 0.12 * di$inst + rnorm(n)
  di$y <- 1 + di$en + 0.5 * di$x + rnorm(n)

  X <- model.matrix(~ en + x, di)
  Z <- model.matrix(~ inst + x, di)
  xhat <- Z %*% solve(crossprod(Z), crossprod(Z, X))
  bread <- solve(crossprod(xhat))
  h_second_stage <- rowSums((xhat %*% bread) * xhat)
  h_influence <- rowSums((xhat %*% bread) * X)

  expect_true(all(h_second_stage >= 0 & h_second_stage <= 1))
  expect_gt(max(h_influence), 1)

  # Both still sum to the number of parameters, so the bound is the only thing
  # separating them.
  expect_equal(sum(h_second_stage), ncol(X))
  expect_equal(sum(h_influence), ncol(X))

  # iv_robust returns finite standard errors on this data. sandwich cannot.
  fit <- iv_robust(y ~ en + x | inst + x, data = di, se_type = "HC2")
  expect_true(all(is.finite(fit$std.error)))

  # Under the influence convention HC2 is not merely inaccurate here, it is
  # undefined: 1 - h goes negative, so its square root does not exist. Asserted
  # by hand rather than through a package, since which convention a given 2SLS
  # object dispatches to depends on which packages are loaded.
  expect_true(any(1 - h_influence < 0))
  expect_true(any(is.nan(suppressWarnings(sqrt(1 - h_influence)))))
})

# The size of the gap, pinned on data outside this package's control.
#
# The vignette quotes these figures, so they are measured here rather than
# asserted there. mtcars is the example the ivreg documentation itself uses,
# and both conventions are well defined on it: the influence diagonal reaches
# -0.0101 but never exceeds 1, so nothing is NaN and the comparison is between
# two computable answers rather than between an answer and a failure.

test_that("the two conventions differ by 8.6% at HC2 and 18.5% at HC3 on mtcars", {
  X <- model.matrix(~ hp + am, mtcars)
  Z <- model.matrix(~ wt + gear, mtcars)
  xhat <- Z %*% solve(crossprod(Z), crossprod(Z, X))
  bread <- solve(crossprod(xhat))

  h_second_stage <- rowSums((xhat %*% bread) * xhat)
  h_influence <- rowSums((xhat %*% bread) * X)

  # Both are well defined here, and the influence diagonal goes negative.
  expect_true(all(h_second_stage >= 0 & h_second_stage <= 1))
  expect_lt(min(h_influence), 0)
  expect_lt(max(h_influence), 1)

  resid <- as.vector(mtcars$mpg - X %*% (bread %*% crossprod(xhat, mtcars$mpg)))
  se_of <- function(h, power) {
    meat <- crossprod(xhat, (resid^2 / (1 - h)^power) * xhat)
    sqrt(diag(bread %*% meat %*% bread))
  }

  # HC2 divides the squared residual by (1 - h), HC3 by (1 - h)^2. Bounds are
  # set from the measured spread across the three coefficients, 7.4-8.6% and
  # 16.5-18.5%, so a drift in either direction fails.
  for (spec in list(list(p = 1, lo = 0.07, hi = 0.09),
                    list(p = 2, lo = 0.16, hi = 0.19))) {
    ratio <- 1 - se_of(h_second_stage, spec$p) / se_of(h_influence, spec$p)
    expect_true(all(ratio > spec$lo & ratio < spec$hi))
  }
})
