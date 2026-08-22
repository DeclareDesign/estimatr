library(estimatr)

# Against Stata, from output frozen in 2019.
#
# Stata is the only reference here that is not an R package, so it is the only
# one that can corroborate the `se_type = "stata"` conventions, the small-sample
# corrections in `areg`, and the `ivregress 2sls` diagnostics. The output was
# produced by the do-files in data-raw/stata/ and has been carried unchanged
# since; the data is `mtcars`, so nothing about the comparison can drift.
#
# The fixtures were removed from the package in 2019 (commit 6219de7) and are
# restored here. They cost nothing to keep and are irreplaceable: regenerating
# them needs a Stata licence, and the divergences they pin are ones no R
# reference can show.
#
# Tolerances are derived per value from the number of digits Stata printed;
# see stata_rel_tol() in helper-external.R.

d <- ext_data_stata()

# Stata orders coefficients with the covariates first and `_cons` last; R puts
# the intercept first. Every comparison below names its coefficient rather than
# relying on position, because the position mapping is exactly the kind of
# thing that is wrong silently.

# ---- reg ----

lm_fits <- list(
  classical  = lm_robust(mpg ~ hp, data = d, se_type = "classical"),
  HC1        = lm_robust(mpg ~ hp, data = d, se_type = "HC1"),
  HC2        = lm_robust(mpg ~ hp, data = d, se_type = "HC2"),
  HC3        = lm_robust(mpg ~ hp, data = d, se_type = "HC3"),
  stata_cl   = lm_robust(mpg ~ hp, clusters = cyl, data = d, se_type = "stata"),
  classicalw = lm_robust(mpg ~ hp, data = d, weights = w, se_type = "classical"),
  HC1w       = lm_robust(mpg ~ hp, data = d, weights = w, se_type = "HC1"),
  HC2w       = lm_robust(mpg ~ hp, data = d, weights = w, se_type = "HC2"),
  HC3w       = lm_robust(mpg ~ hp, data = d, weights = w, se_type = "HC3"),
  stata_clw  = lm_robust(mpg ~ hp, clusters = cyl, weights = w, data = d, se_type = "stata")
)

# The two configurations where this package deliberately does not reproduce
# Stata. They are asserted below rather than dropped.
WLS_HAT_ROWS <- c("HC2w", "HC3w")

test_that("lm_robust reproduces Stata's reg", {
  st <- read_stata_fixture("stata-ests.txt", c("model", "v_hp", "v_cons", "df", "F"))
  expect_identical(st$model, names(lm_fits))

  for (i in seq_len(nrow(st))) {
    nm <- st$model[i]
    if (nm %in% WLS_HAT_ROWS) next
    fit <- lm_fits[[nm]]
    expect_equal_stata(fit$vcov["hp", "hp"], st$v_hp[i], paste0(nm, ": V(hp)"))
    expect_equal_stata(fit$vcov["(Intercept)", "(Intercept)"], st$v_cons[i],
                       paste0(nm, ": V(_cons)"))
    expect_equal_stata(fit$fstatistic[1], st$F[i], paste0(nm, ": F"))
    expect_equal(unname(fit$df[1]), as.numeric(st$df[i]), label = paste0(nm, ": df"))
  }
})

# ---- the weighted-leverage divergence ----
#
# With weights, Stata's vce(hc2) and vce(hc3) build the hat matrix differently
# from R. This package follows R: test_vs_sandwich.R asserts that the same two
# fits agree with sandwich::vcovHC to machine precision, so the divergence is a
# choice about whose convention to follow rather than an error in either.
#
# Bounding the difference from both sides is the point. A bare `expect_true(not
# equal)` would still pass if the two answers drifted apart by a factor of ten,
# and the reason to keep this test at all is to notice if that happens.

test_that("weighted HC2 and HC3 differ from Stata by the documented amount", {
  st <- read_stata_fixture("stata-ests.txt", c("model", "v_hp", "v_cons", "df", "F"))
  for (nm in WLS_HAT_ROWS) {
    i <- match(nm, st$model)
    fit <- lm_fits[[nm]]
    rel <- abs(fit$vcov["hp", "hp"] - as.numeric(st$v_hp[i])) / abs(as.numeric(st$v_hp[i]))
    expect_gt(rel, stata_rel_tol(st$v_hp[i]))
    expect_lt(rel, 0.02)
  }
})

test_that("the weighted-leverage divergence is confined to HC2 and HC3", {
  # Weighted classical, HC1 and clustered fits are in the main table above and
  # match Stata. Stating it here as well means the scope of the divergence is
  # itself pinned: if it ever spread to HC1, this fails with a clear name.
  st <- read_stata_fixture("stata-ests.txt", c("model", "v_hp", "v_cons", "df", "F"))
  for (nm in c("classicalw", "HC1w", "stata_clw")) {
    i <- match(nm, st$model)
    expect_equal_stata(lm_fits[[nm]]$vcov["hp", "hp"], st$v_hp[i], nm)
  }
})

# ---- areg ----

test_that("lm_robust with fixed_effects reproduces Stata's areg", {
  st <- read_stata_fixture("stata-fe-ests.txt", c("model", "v_hp", "F"))
  fe_fits <- list(
    classical  = lm_robust(mpg ~ hp, fixed_effects = ~ carb, data = d, se_type = "classical"),
    HC1        = lm_robust(mpg ~ hp, fixed_effects = ~ carb, data = d, se_type = "HC1"),
    stata_cl   = lm_robust(mpg ~ hp, fixed_effects = ~ carb, clusters = cyl, data = d,
                           se_type = "stata"),
    classicalw = lm_robust(mpg ~ hp, fixed_effects = ~ carb, data = d, weights = w,
                           se_type = "classical"),
    HC1w       = lm_robust(mpg ~ hp, fixed_effects = ~ carb, data = d, weights = w,
                           se_type = "HC1"),
    stata_clw  = lm_robust(mpg ~ hp, fixed_effects = ~ carb, clusters = cyl, weights = w,
                           data = d, se_type = "stata")
  )
  expect_identical(st$model, names(fe_fits))

  for (i in seq_len(nrow(st))) {
    nm <- st$model[i]
    expect_equal_stata(fe_fits[[nm]]$vcov["hp", "hp"], st$v_hp[i], paste0("areg ", nm))
  }

  # Stata reports an F for areg and this package returns none for any fit with
  # absorbed effects. That is inherited from estimatr 1.0.6, which also returned
  # NULL here, so it is recorded as current behaviour rather than compared.
  expect_null(fe_fits$HC1$fstatistic)
})

# ---- ivregress 2sls ----

test_that("iv_robust reproduces Stata's ivregress 2sls", {
  st <- read_stata_fixture(
    "stata-iv-ests.txt",
    c("model", "v_hp", "v_am", "v_cons", "F", "r2", "r2a", "rmse")
  )
  iv_fits <- list(
    classical   = iv_robust(mpg ~ hp + am | wt + gear, data = d, se_type = "classical"),
    rob         = iv_robust(mpg ~ hp + am | wt + gear, data = d, se_type = "HC1"),
    cl          = iv_robust(mpg ~ hp + am | wt + gear, clusters = cyl, data = d,
                            se_type = "stata"),
    classical_w = iv_robust(mpg ~ hp + am | wt + gear, data = d, weights = w,
                            se_type = "classical"),
    rob_w       = iv_robust(mpg ~ hp + am | wt + gear, data = d, weights = w,
                            se_type = "HC1"),
    cl_w        = iv_robust(mpg ~ hp + am | wt + gear, clusters = cyl, weights = w,
                            data = d, se_type = "stata")
  )
  expect_identical(st$model, names(iv_fits))

  for (i in seq_len(nrow(st))) {
    nm <- st$model[i]
    fit <- iv_fits[[nm]]
    expect_equal_stata(fit$vcov["hp", "hp"], st$v_hp[i], paste0(nm, ": V(hp)"))
    expect_equal_stata(fit$vcov["am", "am"], st$v_am[i], paste0(nm, ": V(am)"))
    expect_equal_stata(fit$vcov["(Intercept)", "(Intercept)"], st$v_cons[i],
                       paste0(nm, ": V(_cons)"))
    expect_equal_stata(fit$fstatistic[1], st$F[i], paste0(nm, ": F"))
    expect_equal_stata(fit$r.squared, st$r2[i], paste0(nm, ": r2"))
    expect_equal_stata(fit$adj.r.squared, st$r2a[i], paste0(nm, ": adj r2"))
    # Unweighted RMSE only; see the weighted case below.
    if (!grepl("_w$", nm)) {
      expect_equal_stata(sqrt(fit$res_var), st$rmse[i], paste0(nm, ": rmse"))
    }
  }
})

test_that("weighted 2SLS root MSE follows AER::ivreg rather than Stata", {
  # Stata rescales the weights when forming the root mean squared error for
  # aweighted 2SLS; R's ivreg does not, and this package follows ivreg. As with
  # the leverage case above, the divergence is pinned from both sides: equal to
  # the R reference, and different from Stata by a bounded amount.
  skip_if_not_installed("AER")
  st <- read_stata_fixture(
    "stata-iv-ests.txt",
    c("model", "v_hp", "v_am", "v_cons", "F", "r2", "r2a", "rmse")
  )
  fit <- iv_robust(mpg ~ hp + am | wt + gear, data = d, weights = w, se_type = "classical")
  aer <- AER::ivreg(mpg ~ hp + am | wt + gear, data = d, weights = w)

  expect_equal(sqrt(fit$res_var), aer$sigma, tolerance = LIVE_TOL)

  i <- match("classical_w", st$model)
  rel <- abs(sqrt(fit$res_var) - as.numeric(st$rmse[i])) / as.numeric(st$rmse[i])
  expect_gt(rel, stata_rel_tol(st$rmse[i]))
  expect_lt(rel, 0.25)
})
