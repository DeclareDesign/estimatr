library(estimatr)

# Rank detection must not depend on the units a column is measured in.
# stats::lm() gets this from LINPACK dqrdc2, which compares each column's
# remaining norm against its own original norm. Eigen's setThreshold()
# compares every pivot against the largest pivot in the matrix, so before
# lm_solver() and getMeatXtX() normalized their columns, one regressor in
# large units pushed the others under the threshold and they were dropped as
# collinear on a full-rank design: a wrong answer with a warning, not an
# error. The property below is stronger than any fixed design, and it is what
# the rest of the suite does not check.

set.seed(42)
n <- 200
dat <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  x3 = rnorm(n),
  z = rnorm(n),
  cl = rep(1:20, each = 10)
)
dat$y <- 1 + dat$x1 + dat$x2 + dat$x3 + rnorm(n)

powers <- 0:12

# Refit with x2 multiplied by 10^k, then undo the scaling on the x2 row. Every
# coefficient and standard error must come back to what the unscaled fit gave.
expect_scale_invariant <- function(fitter) {
  base <- fitter(dat)
  for (k in powers) {
    scaled <- dat
    scaled$x2 <- scaled$x2 * 10^k
    fit <- fitter(scaled)

    coefs <- coef(fit)
    ses <- fit$std.error
    coefs["x2"] <- coefs["x2"] * 10^k
    ses["x2"] <- ses["x2"] * 10^k

    expect_false(anyNA(coefs), label = paste("no coefficient dropped at k =", k))
    expect_equal(coefs, coef(base), tolerance = 1e-8,
                 info = paste("coefficients at k =", k))
    expect_equal(ses, base$std.error, tolerance = 1e-8,
                 info = paste("standard errors at k =", k))
  }
}

test_that("lm_robust is invariant to column scaling", {
  expect_scale_invariant(function(d) lm_robust(y ~ x1 + x2 + x3, data = d))
})

test_that("clustered lm_robust is invariant to column scaling", {
  expect_scale_invariant(
    function(d) lm_robust(y ~ x1 + x2 + x3, clusters = cl, data = d)
  )
})

# getMeatXtX() is the path HC2, HC3, and CR2 read the hat values off, and it
# carries the same threshold as lm_solver(). A fix applied to only one of the
# two leaves the variance read off a rank the coefficients were not fitted at.
test_that("every hat-value se_type is invariant to column scaling", {
  for (se_type in c("HC0", "HC1", "HC2", "HC3", "classical")) {
    local({
      this_type <- se_type
      expect_scale_invariant(
        function(d) lm_robust(y ~ x1 + x2 + x3, se_type = this_type, data = d)
      )
    })
  }
  expect_scale_invariant(
    function(d) lm_robust(y ~ x1 + x2 + x3, clusters = cl,
                          se_type = "CR0", data = d)
  )
})

test_that("iv_robust is invariant to column scaling", {
  expect_scale_invariant(
    function(d) iv_robust(y ~ x1 + x2 + x3 | x1 + z + x3, data = d)
  )
  expect_scale_invariant(
    function(d) iv_robust(y ~ x1 + x2 + x3 | x1 + z + x3,
                          clusters = cl, data = d)
  )
})

# The Cholesky path is a second rank determination, and Eigen's LLT reports
# success on a numerically singular Gram matrix, so info() alone never caught
# a rank-deficient design. Normalizing the columns makes each L_ii the
# column's own residual norm, which is dqrdc2's test, and the path falls back
# to the QR below the same 1e-7.
test_that("try_cholesky reaches the same answer as the QR path", {
  expect_scale_invariant(
    function(d) lm_robust(y ~ x1 + x2 + x3, data = d, try_cholesky = TRUE)
  )

  qr_fit <- lm_robust(y ~ x1 + x2 + x3, data = dat, try_cholesky = FALSE)
  ch_fit <- lm_robust(y ~ x1 + x2 + x3, data = dat, try_cholesky = TRUE)
  expect_equal(coef(ch_fit), coef(qr_fit), tolerance = 1e-10)
  expect_equal(ch_fit$std.error, qr_fit$std.error, tolerance = 1e-10)
})

test_that("try_cholesky still drops an exactly collinear column", {
  dup <- dat
  dup$x2 <- dup$x1
  fit <- suppressWarnings(
    lm_robust(y ~ x1 + x2 + x3, data = dup, try_cholesky = TRUE)
  )
  expect_equal(sum(is.na(coef(fit))), 1L)
  expect_equal(unname(coef(fit)),
               unname(coef(lm(y ~ x1 + x2 + x3, data = dup))),
               tolerance = 1e-10)

  const <- data.frame(x = rep(1, n), z = rnorm(n), y = rnorm(n))
  cfit <- suppressWarnings(lm_robust(y ~ x + z, data = const, try_cholesky = TRUE))
  expect_equal(sum(is.na(coef(cfit))), 1L)
})

test_that("a full-rank design in large units agrees with lm", {
  d <- dat
  d$x2 <- d$x2 * 1e9
  fit <- lm_robust(y ~ x1 + x2 + x3, data = d, se_type = "classical")
  ref <- lm(y ~ x1 + x2 + x3, data = d)
  expect_false(anyNA(coef(fit)))
  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-10)
  expect_equal(unname(fit$std.error),
               unname(summary(ref)$coefficients[, 2]), tolerance = 1e-10)
})

# The normalization must not cost the rank deficiency the 1e-7 threshold was
# added to catch (#351, #395): an exactly collinear column can otherwise
# survive as a pivot of order 1e-14 and produce a coefficient of order 1e11
# where lm() gives NA.
test_that("exactly collinear columns are still dropped", {
  const <- data.frame(x = rep(1, n), z = rnorm(n), y = rnorm(n))
  expect_warning(lm_robust(y ~ x + z, data = const), "collinear")
  fit_const <- suppressWarnings(lm_robust(y ~ x + z, data = const))
  expect_equal(sum(is.na(coef(fit_const))), 1L)
  expect_equal(sum(is.na(coef(fit_const))),
               sum(is.na(coef(lm(y ~ x + z, data = const)))))

  dup <- dat
  dup$x2 <- dup$x1
  expect_warning(lm_robust(y ~ x1 + x2 + x3, data = dup), "collinear")
  fit_dup <- suppressWarnings(lm_robust(y ~ x1 + x2 + x3, data = dup))
  expect_equal(sum(is.na(coef(fit_dup))), 1L)
  expect_equal(sum(is.na(coef(fit_dup))),
               sum(is.na(coef(lm(y ~ x1 + x2 + x3, data = dup)))))
})

test_that("a rescaled collinear column is still dropped", {
  dup <- dat
  dup$x2 <- dup$x1 * 1e9
  expect_warning(lm_robust(y ~ x1 + x2 + x3, data = dup), "collinear")
  fit <- suppressWarnings(lm_robust(y ~ x1 + x2 + x3, data = dup))
  expect_equal(sum(is.na(coef(fit))), 1L)
  expect_equal(sum(is.na(coef(fit))),
               sum(is.na(coef(lm(y ~ x1 + x2 + x3, data = dup)))))
})

test_that("an all-zero column does not divide by its norm", {
  zeroed <- dat
  zeroed$x2 <- 0
  expect_warning(lm_robust(y ~ x1 + x2 + x3, data = zeroed), "collinear")
  fit <- suppressWarnings(lm_robust(y ~ x1 + x2 + x3, data = zeroed))
  expect_true(is.na(coef(fit)[["x2"]]))
  expect_false(anyNA(coef(fit)[c("(Intercept)", "x1", "x3")]))
  expect_equal(unname(coef(fit)),
               unname(coef(lm(y ~ x1 + x2 + x3, data = zeroed))),
               tolerance = 1e-10)
})

# Multiplicity two. getMeatXtX() removes the QR's tossed columns by in-place
# left shift in pivot order, which is correct only in descending order, and
# with one dropped column there is nothing to order. Every probe recorded as
# verification of that code until now was rank deficient by exactly one, so
# none of them tested the ordering at all.
#
# The reference has to be independent rather than the other path: at exact
# rank deficiency the Cholesky path falls back to the QR, so the two paths
# run the same code and their agreement proves nothing. Making x2 and x3
# exact copies of x1 fixes the comparison instead: whichever column survives
# the pivot carries x1's column, so the surviving estimate and its standard
# error must equal the reduced fit's, for every se_type that reads the meat.
test_that("rank deficiency of two keeps the reduced fit's answers", {
  dup2 <- dat
  dup2$x2 <- dup2$x1
  dup2$x3 <- dup2$x1

  ref_lm <- lm(y ~ x1 + x2 + x3, data = dup2)
  expect_equal(sum(is.na(coef(ref_lm))), 2L)

  for (this_type in c("classical", "HC0", "HC1", "HC2", "HC3")) {
    full <- suppressWarnings(
      lm_robust(y ~ x1 + x2 + x3, data = dup2, se_type = this_type)
    )
    reduced <- lm_robust(y ~ x1, data = dup2, se_type = this_type)

    expect_equal(sum(is.na(coef(full))), 2L,
                 info = paste("dropped count,", this_type))
    kept <- !is.na(coef(full))
    expect_equal(unname(coef(full)[kept]), unname(coef(reduced)),
                 tolerance = 1e-10, info = paste("coefficients,", this_type))
    expect_equal(unname(full$std.error[kept]), unname(reduced$std.error),
                 tolerance = 1e-10, info = paste("standard errors,", this_type))
  }

  for (this_type in c("CR0", "CR2")) {
    full <- suppressWarnings(
      lm_robust(y ~ x1 + x2 + x3, data = dup2, clusters = cl,
                se_type = this_type)
    )
    reduced <- lm_robust(y ~ x1, data = dup2, clusters = cl,
                         se_type = this_type)
    kept <- !is.na(coef(full))
    expect_equal(sum(is.na(coef(full))), 2L,
                 info = paste("dropped count,", this_type))
    expect_equal(unname(full$std.error[kept]), unname(reduced$std.error),
                 tolerance = 1e-10, info = paste("standard errors,", this_type))
  }

  # Fitted values are unique whichever columns the pivot keeps, so they are
  # the one thing that can be held against lm() without pinning B14.
  full <- suppressWarnings(lm_robust(y ~ x1 + x2 + x3, data = dup2))
  expect_equal(unname(fitted(full)), unname(fitted(ref_lm)), tolerance = 1e-10)
})

# A deficiency of two built from two different causes, so the tossed columns
# are not adjacent in the original ordering.
test_that("a constant column beside a duplicate is dropped as two", {
  mixed <- dat
  mixed$x2 <- 1
  mixed$x3 <- mixed$x1
  fit <- suppressWarnings(lm_robust(y ~ x1 + x2 + x3, data = mixed))
  expect_equal(sum(is.na(coef(fit))), 2L)
  expect_equal(sum(is.na(coef(fit))),
               sum(is.na(coef(lm(y ~ x1 + x2 + x3, data = mixed)))))
  expect_equal(unname(fitted(fit)),
               unname(fitted(lm(y ~ x1 + x2 + x3, data = mixed))),
               tolerance = 1e-10)
})

# The rank decision across the whole approach to singularity, rather than at
# one design. Every other rank test here is either exactly deficient or
# comfortably full rank; the interesting region is between them, and it is
# where a retuned threshold would first show up. The two paths and lm() must
# agree about where rank is lost at every step, which is the property that
# makes try_cholesky safe to expose: the argument selects the arithmetic and
# never the answer.
#
# Condition indices are computed on the column-scaled design, as Belsley, Kuh
# and Welsch do, because that is the quantity the normalization in front of
# both decompositions makes relevant. The raw condition number of a design in
# mixed units is large for a reason that does not affect the fit.
test_that("both paths track lm's rank decision as collinearity approaches", {
  set.seed(99)
  m <- 400
  for (eps in 10^-(1:12)) {
    a <- rnorm(m)
    d <- data.frame(x1 = a, x2 = rnorm(m), x3 = a + eps * rnorm(m))
    d$y <- rnorm(m)

    ref <- lm(y ~ x1 + x2 + x3, data = d)
    qr_fit <- suppressWarnings(
      lm_robust(y ~ x1 + x2 + x3, data = d, try_cholesky = FALSE)
    )
    ch_fit <- suppressWarnings(
      lm_robust(y ~ x1 + x2 + x3, data = d, try_cholesky = TRUE)
    )

    lab <- paste("eps =", format(eps))
    expect_equal(sum(is.na(coef(qr_fit))), sum(is.na(coef(ref))),
                 info = paste("QR path against lm at", lab))
    expect_equal(sum(is.na(coef(ch_fit))), sum(is.na(coef(ref))),
                 info = paste("Cholesky path against lm at", lab))

    # Well inside the safe zone the two paths must agree to far more digits
    # than anything reportable. Measured on this design, the worst relative
    # gap over eps >= 1e-3 is 4.2e-10, at a scaled condition index of 1,800,
    # so 1e-8 keeps about 24x of headroom for the platforms where the linear
    # algebra differs. The band covers every design this package is built
    # for by a wide margin: a blocked trial with covariates measures about
    # 38. Do not widen it to eps = 1e-4 without re-measuring, since the gap
    # there is 6.2e-8 and the assertion would be tighter than the arithmetic.
    if (eps >= 1e-3) {
      expect_equal(coef(ch_fit), coef(qr_fit), tolerance = 1e-8,
                   info = paste("coefficients at", lab))
      expect_equal(ch_fit$std.error, qr_fit$std.error, tolerance = 1e-8,
                   info = paste("standard errors at", lab))
    }
  }
})

# The fallback is shared code, but each estimator reaches it through its own
# fit call, and until now only unclustered, unweighted lm_robust() exercised
# it. A wiring mistake at one entry point would be invisible everywhere else.
test_that("the Cholesky path is wired correctly at every entry point", {
  dup <- dat
  dup$x2 <- dup$x1
  dup$w <- runif(n, 0.5, 2)
  dup$Z <- rep(c(0, 1), length.out = n)

  same_both_ways <- function(fitter, label) {
    q <- suppressWarnings(fitter(FALSE))
    ch <- suppressWarnings(fitter(TRUE))
    expect_equal(coef(ch), coef(q), tolerance = 1e-10,
                 info = paste("coefficients,", label))
    expect_equal(ch$std.error, q$std.error, tolerance = 1e-10,
                 info = paste("standard errors,", label))
    expect_equal(sum(is.na(coef(ch))), 1L,
                 info = paste("dropped count,", label))
  }

  same_both_ways(function(tc) lm_robust(y ~ x1 + x2 + x3, data = dup,
                                        clusters = cl, se_type = "CR2",
                                        try_cholesky = tc), "clustered CR2")
  same_both_ways(function(tc) lm_robust(y ~ x1 + x2 + x3, data = dup,
                                        weights = w, try_cholesky = tc),
                 "weighted")
  same_both_ways(function(tc) iv_robust(y ~ x1 + x2 + x3 | x1 + x2 + z,
                                        data = dup, try_cholesky = tc),
                 "iv_robust")

  # lm_lin builds its own centered interactions, so the design it hands the
  # solver is not the one written in the formula.
  lin_q <- lm_lin(y ~ Z, ~ x1 + x3, data = dup, try_cholesky = FALSE)
  lin_ch <- lm_lin(y ~ Z, ~ x1 + x3, data = dup, try_cholesky = TRUE)
  expect_equal(coef(lin_ch), coef(lin_q), tolerance = 1e-10)
  expect_equal(lin_ch$std.error, lin_q$std.error, tolerance = 1e-10)

  # A full-rank fixed-effects fit never reaches the fallback, so it checks
  # that the fast path is right where it actually runs.
  fe_q <- lm_robust(y ~ x1 + x3, data = dat, fixed_effects = ~ cl,
                    try_cholesky = FALSE)
  fe_ch <- lm_robust(y ~ x1 + x3, data = dat, fixed_effects = ~ cl,
                     try_cholesky = TRUE)
  expect_equal(coef(fe_ch), coef(fe_q), tolerance = 1e-10)
  expect_equal(fe_ch$std.error, fe_q$std.error, tolerance = 1e-10)
})
