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
