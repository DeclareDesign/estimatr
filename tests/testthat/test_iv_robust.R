library(estimatr)

# iv_robust's own behaviour: the 2SLS estimate, its residuals, its diagnostics,
# and the regressions filed against it.
#
# HC0 through HC3, CR0 and CR2 against sandwich and clubSandwich are in
# test_vs_sandwich.R and test_vs_clubsandwich.R; Stata's ivregress in
# test_vs_stata.R; reparameterisations and respelled instruments in
# test_invariance.R; underidentified models in test_degenerate.R.

iv_data <- function() {
  set.seed(42)
  n <- 200
  d <- data.frame(
    x = rnorm(n),
    w = rnorm(n),
    inst = rnorm(n),
    inst2 = rnorm(n),
    wt = runif(n, 0.5, 2)
  )
  d$en <- d$inst + 0.5 * d$inst2 + rnorm(n)
  d$y <- 1 + d$en + d$w + rnorm(n)
  d
}

d <- iv_data()

# Measured gaps against AER below are at most 4e-14.
IV_TOL <- 1e-10

test_that("2SLS matches AER::ivreg, with an exogenous covariate and with weights", {
  # Classical, where AER's own vcov() is the same estimator. The robust
  # variances are compared through sandwich in test_vs_sandwich.R.
  skip_if_not_installed("AER")
  fit <- iv_robust(y ~ en + w | inst + w, data = d, se_type = "classical")
  ref <- AER::ivreg(y ~ en + w | inst + w, data = d)
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = IV_TOL)
  expect_equal(unname(fit$vcov), unname(vcov(ref)), tolerance = IV_TOL)

  fit <- iv_robust(y ~ en + w | inst + w, data = d, se_type = "classical", weights = wt)
  ref <- AER::ivreg(y ~ en + w | inst + w, data = d, weights = wt)
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = IV_TOL)
  expect_equal(unname(fit$vcov), unname(vcov(ref)), tolerance = IV_TOL)
})

test_that("just identified with one instrument, the slope is the ratio of covariances", {
  fit <- iv_robust(y ~ en | inst, data = d)
  expect_equal(fit$coefficients[["en"]], cov(d$y, d$inst) / cov(d$en, d$inst),
               tolerance = IV_TOL)
})

test_that("#345: iv_robust returns structural residuals", {
  m <- iv_robust(mpg ~ wt | am, data = mtcars)
  X <- cbind(1, mtcars$wt)
  expect_equal(unname(m$residuals), as.vector(mtcars$mpg - X %*% coef(m)),
               tolerance = 1e-12)
})

# ---- diagnostics ----

test_that("classical diagnostics match AER: weak instruments, Wu-Hausman, Sargan", {
  skip_if_not_installed("AER")
  fit <- iv_robust(y ~ en + w | inst + inst2 + w, data = d, se_type = "classical",
                   diagnostics = TRUE)
  ref <- summary(AER::ivreg(y ~ en + w | inst + inst2 + w, data = d),
                 diagnostics = TRUE)$diagnostics

  weak <- fit$diagnostic_first_stage_fstatistic
  expect_equal(weak[["value"]], ref["Weak instruments", "statistic"], tolerance = IV_TOL)
  expect_equal(weak[["nomdf"]], ref["Weak instruments", "df1"])
  expect_equal(weak[["dendf"]], ref["Weak instruments", "df2"])
  expect_equal(weak[["p.value"]], ref["Weak instruments", "p-value"], tolerance = IV_TOL)

  wu <- fit$diagnostic_endogeneity_test
  expect_equal(wu[["value"]], ref["Wu-Hausman", "statistic"], tolerance = IV_TOL)
  expect_equal(wu[["numdf"]], ref["Wu-Hausman", "df1"])
  expect_equal(wu[["dendf"]], ref["Wu-Hausman", "df2"])
  expect_equal(wu[["p.value"]], ref["Wu-Hausman", "p-value"], tolerance = IV_TOL)

  overid <- fit$diagnostic_overid_test
  expect_equal(overid[["value"]], ref["Sargan", "statistic"], tolerance = IV_TOL)
  expect_equal(overid[["df"]], ref["Sargan", "df1"])
  expect_equal(overid[["p.value"]], ref["Sargan", "p-value"], tolerance = IV_TOL)
})

test_that("over-identified diagnostics work for every se_type, not just classical", {
  # `first_stage_fits` has its columns renamed `fit_<endog>`, and the rewritten
  # robust branch then indexed it by the bare endogenous names, so every
  # over-identified fit with a non-classical se_type died with "subscript out
  # of bounds". The classical branch took a different function and was fine,
  # which is why the whole path had no test.
  set.seed(2)
  n <- 500
  dd <- data.frame(z1 = rnorm(n), z2 = rnorm(n), z3 = rnorm(n), w = rnorm(n))
  dd$x <- dd$z1 + 0.5 * dd$z2 + 0.3 * dd$z3 + rnorm(n)
  dd$y <- dd$x + dd$w + rnorm(n)
  fml <- y ~ x + w | z1 + z2 + z3 + w

  for (ty in c("classical", "HC0", "HC1", "HC2", "HC3")) {
    m <- iv_robust(fml, data = dd, se_type = ty, diagnostics = TRUE)
    ov <- m$diagnostic_overid_test
    expect_equal(names(ov), c("value", "df", "p.value"), info = ty)
    expect_true(is.finite(ov[["value"]]), info = ty)
    expect_equal(unname(ov[["df"]]), 2, info = ty)   # 3 instruments, 1 endogenous
  }

  # The robust branch is Wooldridge's score test, a different statistic from
  # Sargan's, and every non-classical se_type shares it.
  classical <- iv_robust(fml, data = dd, se_type = "classical", diagnostics = TRUE)
  rb <- vapply(c("HC0", "HC1", "HC2", "HC3"), function(ty) {
    iv_robust(fml, data = dd, se_type = ty, diagnostics = TRUE)$diagnostic_overid_test[["value"]]
  }, numeric(1))
  expect_equal(unname(diff(range(rb))), 0)
  expect_false(isTRUE(all.equal(rb[["HC0"]], classical$diagnostic_overid_test[["value"]])))
})

test_that("#389: glance() works with multiple endogenous regressors", {
  m <- iv_robust(mpg ~ hp + wt | am + cyl, data = mtcars, diagnostics = TRUE)
  g <- glance(m)
  expect_equal(nrow(g), 1L)
  # reports the weakest of the per-regressor first stages
  fs <- m$diagnostic_first_stage_fstatistic
  expect_equal(g$statistic.weakinst, unname(min(fs[grep("(^|:)value$", names(fs))])))
})

test_that("#397: model.frame() on iv_robust returns the model variables", {
  # Called directly: loading clubSandwich registers a model.frame method for
  # this class too, so bare dispatch depends on what the suite loaded first.
  m <- iv_robust(mpg ~ wt + hp | am + hp, data = mtcars)
  mf <- estimatr:::model.frame.iv_robust(m)
  expect_equal(nrow(mf), nrow(mtcars))
  expect_setequal(names(mf), c("mpg", "wt", "hp", "am"))
})
