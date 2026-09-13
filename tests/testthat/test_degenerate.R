library(estimatr)

# Designs at the edge of what can be estimated, and what each estimator must
# do there.
#
# Three answers are acceptable at a degenerate design, in order of preference:
# the answer lm() or the reduced model gives, an error that says what is wrong,
# or an answer with a warning that says why it cannot be trusted. What is not
# acceptable is a clean-looking fit, and every refusal below exists because
# 1.0.6 or an earlier 2.0 returned one.
#
# Rank deficiency in lm_robust itself, at multiplicity one and two and across
# the approach to singularity, is in test_rank_scale_invariance.R; leverage at
# or near 1 is in test_lm_robust.R under #395; small and singleton blocks are
# in test_blocked_variance.R; singleton fixed-effect groups are in
# test_fe_leverage.R.
#
# A rank-deficient fit and its reduced model run the same arithmetic on the
# same columns, so they agree to rounding: every comparison below holds at
# 1e-15, and DEG_TOL leaves the cross-platform room the other files leave.
DEG_TOL <- 1e-10

deg_data <- function() {
  set.seed(3)
  n <- 80
  d <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    z = rep(0:1, 40),
    cl = rep(1:16, 5),
    g = factor(rep(1:8, each = 10)),
    w = runif(n, 0.5, 2),
    inst = rnorm(n)
  )
  d$en <- d$inst + rnorm(n)
  d$y <- 1 + d$x1 + d$en + rnorm(n)
  d$dup <- d$x1
  d$dup2 <- 2 * d$x1
  d$dup3 <- d$x2 - d$x1
  d$inst_copy <- d$inst
  d$const <- 5
  d$g_level <- as.numeric(d$g) / 3
  d
}

d <- deg_data()
n <- nrow(d)

# The rank-deficient fit's surviving coefficients against the reduced fit.
expect_reduces_to <- function(full, reduced, label) {
  kept <- names(reduced$coefficients)
  expect_true(all(is.na(full$coefficients[setdiff(names(full$coefficients), kept)])),
              label = paste(label, "dropped terms are NA"))
  expect_equal(full$coefficients[kept], reduced$coefficients, tolerance = DEG_TOL,
               label = paste(label, "coefficients"))
  expect_equal(full$std.error[kept], reduced$std.error, tolerance = DEG_TOL,
               label = paste(label, "standard errors"))
  expect_equal(full$df[kept], reduced$df, tolerance = DEG_TOL, label = paste(label, "df"))
}

# ---- collinearity outside lm_robust's own design ----

test_that("iv_robust with a collinear exogenous regressor is the fit without it", {
  for (st in c("HC2", "classical", "CR2")) {
    cluster_ids <- if (st == "CR2") d$cl else NULL
    expect_warning(
      full <- iv_robust(y ~ en + x1 + dup | inst + x1 + dup, data = d, se_type = st,
                        clusters = cluster_ids),
      "collinear"
    )
    reduced <- iv_robust(y ~ en + x1 | inst + x1, data = d, se_type = st,
                         clusters = cluster_ids)
    expect_reduces_to(full, reduced, st)
  }
})

test_that("a redundant instrument changes nothing, diagnostics included", {
  copied <- iv_robust(y ~ en + x1 | inst + inst_copy + x1, data = d, diagnostics = TRUE)
  single <- iv_robust(y ~ en + x1 | inst + x1, data = d, diagnostics = TRUE)
  expect_reduces_to(copied, single, "copied instrument")
  expect_equal(copied$diagnostic_first_stage_fstatistic,
               single$diagnostic_first_stage_fstatistic, tolerance = DEG_TOL)
  expect_equal(copied$diagnostic_endogeneity_test,
               single$diagnostic_endogeneity_test, tolerance = DEG_TOL)
  # Just identified once the copy is set aside, so there is nothing to test.
  expect_true(is.na(copied$diagnostic_overid_test[["value"]]))
})

test_that("lm_lin with a duplicated or constant covariate is lm_lin without it", {
  reduced <- lm_lin(y ~ z, covariates = ~ x1, data = d)
  expect_warning(
    duplicated <- lm_lin(y ~ z, covariates = ~ x1 + dup, data = d),
    "dup_c, z:dup_c"
  )
  expect_reduces_to(duplicated, reduced, "duplicated covariate")
  expect_warning(
    constant <- lm_lin(y ~ z, covariates = ~ x1 + const, data = d),
    "const_c, z:const_c"
  )
  expect_reduces_to(constant, reduced, "constant covariate")
})

test_that("a regressor the fixed effects absorb is NA, and the rest are unchanged", {
  cases <- list(
    HC2 = list(se_type = "HC2", clusters = NULL),
    CR2 = list(se_type = "CR2", clusters = d$cl)
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    expect_warning(
      full <- lm_robust(y ~ x1 + g_level, data = d, fixed_effects = ~ g,
                        se_type = cs$se_type, clusters = cs$clusters),
      "g_level"
    )
    reduced <- lm_robust(y ~ x1, data = d, fixed_effects = ~ g,
                         se_type = cs$se_type, clusters = cs$clusters)
    expect_reduces_to(full, reduced, nm)
  }
  # The same NA the dummy expansion gives, with the dummies ahead of the
  # regressor as absorption effectively places them; written the other way
  # round, lm() drops a dummy instead.
  expect_true(is.na(coef(lm(y ~ x1 + g + g_level, data = d))[["g_level"]]))
})

test_that("weighted, clustered designs short of full rank by two are the reduced fit", {
  for (st in c("CR2", "stata", "CR0")) {
    expect_warning(
      full <- lm_robust(y ~ x1 + x2 + dup2 + dup3, data = d, weights = w,
                        clusters = cl, se_type = st),
      "dup2, dup3"
    )
    reduced <- lm_robust(y ~ x1 + x2, data = d, weights = w, clusters = cl, se_type = st)
    expect_reduces_to(full, reduced, st)
  }
})

# ---- no estimate exists ----

test_that("an underidentified instrumental-variables model is refused", {
  # Too few instruments by count: 1.0.6, and 2.0 until now, warned, dropped
  # the intercept, and reported hp and cyl.
  expect_error(
    iv_robust(mpg ~ hp + cyl | am, data = mtcars),
    "do not identify every regressor"
  )
  # Enough instruments by count, but one is a copy of the exogenous regressor,
  # so the endogenous regressor has none. Both versions dropped `en` and
  # returned the rest.
  expect_error(
    iv_robust(y ~ en + x1 | x1 + dup, data = d),
    "do not identify every regressor"
  )
  expect_error(
    iv_robust(y ~ en + x1 | x1 + dup, data = d, clusters = cl, se_type = "CR2"),
    "do not identify every regressor"
  )
})

test_that("a fit with no observations left is refused", {
  no_outcome <- d
  no_outcome$y <- NA
  expect_error(lm_robust(y ~ x1, data = no_outcome), "No observations are left")
  expect_error(lm_robust(y ~ x1, data = d, subset = x1 > 100), "No observations are left")
  expect_error(iv_robust(y ~ en | inst, data = no_outcome), "No observations are left")
  expect_error(lm_lin(y ~ z, covariates = ~ x1, data = no_outcome), "No observations are left")
  expect_error(difference_in_means(y ~ z, data = no_outcome), "No observations are left")
  expect_error(difference_in_means(y ~ z, data = d, blocks = g, subset = x1 > 100),
               "No observations are left")
})

test_that("a difference in means with a single unit in an arm is refused in the right words", {
  d1 <- d
  d1$treated_once <- as.integer(seq_len(n) == 5)
  err <- expect_error(difference_in_means(y ~ treated_once, data = d1), "at least two units")
  expect_false(grepl("every block", conditionMessage(err)))

  d1$all_treated <- 1L
  expect_error(difference_in_means(y ~ all_treated, data = d1), "more than one value")
})

# ---- estimates exist, and inference does not ----

test_that("a saturated design returns lm()'s coefficients and says inference failed", {
  saturated <- d[1:3, ]
  for (st in c("classical", "HC1", "HC2")) {
    warnings <- character()
    fit <- withCallingHandlers(
      lm_robust(y ~ x1 + x2, data = saturated, se_type = st),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_equal(unname(fit$coefficients), unname(coef(lm(y ~ x1 + x2, data = saturated))),
                 tolerance = DEG_TOL, label = st)
    expect_true(any(grepl("degrees of freedom have been estimated as negative or zero", warnings)),
                label = paste(st, "warns about the degrees of freedom"))
    expect_true(all(is.na(fit$p.value)), label = paste(st, "p-values"))
  }
})

test_that("more coefficients than observations drops as many as lm() does", {
  tiny <- d[1:2, ]
  warnings <- character()
  fit <- withCallingHandlers(
    lm_robust(y ~ x1 + x2, data = tiny, se_type = "classical"),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(sum(is.na(fit$coefficients)), sum(is.na(coef(lm(y ~ x1 + x2, data = tiny)))))
  expect_true(any(grepl("collinear", warnings)))
  expect_true(any(grepl("degrees of freedom have been estimated as negative or zero", warnings)))
})

test_that("an outcome the regressors fit exactly has vanishing standard errors", {
  d1 <- d
  d1$y_exact <- 1 + 2 * d$x1 - d$x2
  d1$y_const <- 3
  for (st in c("classical", "HC2", "CR2")) {
    cluster_ids <- if (st == "CR2") d1$cl else NULL
    exact <- lm_robust(y_exact ~ x1 + x2, data = d1, se_type = st, clusters = cluster_ids)
    expect_equal(unname(exact$coefficients), c(1, 2, -1), tolerance = DEG_TOL, label = st)
    expect_lt(max(exact$std.error), 1e-12)

    constant <- lm_robust(y_const ~ x1 + x2, data = d1, se_type = st, clusters = cluster_ids)
    expect_equal(unname(constant$coefficients), c(3, 0, 0), tolerance = DEG_TOL, label = st)
    expect_lt(max(constant$std.error), 1e-12)
  }
})
