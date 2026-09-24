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
# Rank deficiency in lm_robust's own design, at multiplicity one and two and
# across the approach to singularity, and leverage at or near 1 (#395), are the
# last sections. Small and singleton blocks are in test_blocked_variance.R;
# singleton fixed-effect groups are in test_fe_leverage.R.
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

# The rank-deficient fit against the reduced fit, over the whole return surface.
#
# This compared coefficients, standard errors and degrees of freedom, and
# nothing else, which is how b96ab8b shipped in 2.0.0: summary()$fstatistic
# came back NA whenever the dropped column was not the last one, under every
# robust se_type and in lm_robust(), lm_lin() and iv_robust() alike, and no
# call site here was looking at a model-level field. test_return_surface.R
# reads every field but its data is full rank, so the defect sat in the cell
# where the two suites cross. The walk below closes it for every existing call
# site at once, and for any estimator added later.
#
# `k` counts the columns asked for where `rank` counts the ones kept, so it
# differs from the reduced fit by construction; the rest are inputs, or carry
# environments.
REDUCES_SKIP <- c("call", "terms", "formula", "terms_regressors", "outcome",
                  "weights", "model", "k")

expect_reduces_to <- function(full, reduced, label) {
  kept <- names(reduced$coefficients)
  expect_true(all(is.na(full$coefficients[setdiff(names(full$coefficients), kept)])),
              label = paste(label, "dropped terms are NA"))
  expect_equal(full$coefficients[kept], reduced$coefficients, tolerance = DEG_TOL,
               label = paste(label, "coefficients"))
  expect_equal(full$std.error[kept], reduced$std.error, tolerance = DEG_TOL,
               label = paste(label, "standard errors"))
  expect_equal(full$df[kept], reduced$df, tolerance = DEG_TOL, label = paste(label, "df"))

  # Everything else the fit returns. A field is coefficient-indexed if it is as
  # long as the full coefficient vector, in which case it is subset to the kept
  # terms; otherwise it has to match outright.
  at <- match(kept, full$term)
  n_checked <- 0L
  for (f in setdiff(names(reduced), REDUCES_SKIP)) {
    a <- full[[f]]
    b <- reduced[[f]]
    if (is.null(a) || is.null(b)) next
    if (!is.numeric(a) && !is.character(a) && !is.logical(a)) next
    lab <- paste(label, f)
    if (is.matrix(a) || is.matrix(b)) {
      if (!identical(dim(a), dim(b))) next
      expect_equal(unname(a), unname(b), tolerance = DEG_TOL, label = lab)
    } else if (length(a) == length(full$coefficients) && length(b) == length(kept)) {
      expect_equal(unname(a[at]), unname(b), tolerance = DEG_TOL, label = lab)
    } else if (length(a) == length(b)) {
      expect_equal(unname(a), unname(b), tolerance = DEG_TOL, label = lab)
    } else {
      next
    }
    n_checked <- n_checked + 1L
  }
  # A walk that silently compared nothing would be a green empty test.
  expect_gt(n_checked, 0)

  sf <- summary(full)
  sr <- summary(reduced)
  expect_false(anyNA(sf$fstatistic), label = paste(label, "F statistic is not NA"))
  expect_equal(unname(sf$fstatistic), unname(sr$fstatistic), tolerance = DEG_TOL,
               label = paste(label, "F statistic"))
}

# ---- collinearity outside lm_robust's own design ----

test_that("iv_robust with a collinear exogenous regressor is the fit without it", {
  for (st in c("HC2", "classical", "CR2")) {
    cluster_ids <- if (st == "CR2") d$cl else NULL
    expect_message(
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
  expect_message(
    duplicated <- lm_lin(y ~ z, covariates = ~ x1 + dup, data = d),
    "dup_c, z:dup_c"
  )
  expect_reduces_to(duplicated, reduced, "duplicated covariate")
  expect_message(
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
    expect_message(
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
    expect_message(
      full <- lm_robust(y ~ x1 + x2 + dup2 + dup3, data = d, weights = w,
                        clusters = cl, se_type = st),
      "dup2, dup3"
    )
    reduced <- lm_robust(y ~ x1 + x2, data = d, weights = w, clusters = cl, se_type = st)
    expect_reduces_to(full, reduced, st)
  }
})

# ---- no estimate exists ----

# An exactly flat first stage: `X` is balanced on `Z` within every level of
# `W`, so the instrument adds nothing and `X`'s fitted values lie in the span
# of the intercept and `W`. This is the DesignLibrary `binary_iv` draw, where a
# binary instrument and a binary regressor give an exactly balanced 2x2 table
# on about 1% of draws at N = 100.
flat_iv_data <- function(n = 80) {
  set.seed(343)
  b <- data.frame(
    W = rep(c(0, 1), each = n / 2),
    Z = rep(rep(c(0, 1), each = n / 4), 2),
    X = rep(rep(c(0, 1), each = n / 8), 4)
  )
  b$Y <- rnorm(n) + 0.5 * b$W
  b
}

test_that("an underidentified regressor is NA and the rest of the fit stands", {
  # Until 2.0.1 each of these was an error. 1.0.6 was worse than either: it
  # warned, let lm_solver's pivot choose, and the pivot took the intercept, so
  # the endogenous regressor came back carrying the intercept's value under its
  # own name. That is the number a reader takes for the LATE, so which column
  # is dropped is an answer here and not a naming question.
  #
  # Too few instruments by count. `hp` keeps a just-identified estimate and
  # `cyl` goes, as stats::lm() drops the later column.
  expect_message(
    full <- iv_robust(mpg ~ hp + cyl | am, data = mtcars),
    "collinear with other regressors"
  )
  expect_reduces_to(full, iv_robust(mpg ~ hp | am, data = mtcars), "too few instruments")

  # Enough instruments by count, but one is a copy of the exogenous regressor,
  # so the endogenous regressor has none. The endogenous `en` sits ahead of the
  # exogenous `x1` in the formula, so column order alone would drop the wrong
  # one; the exogenous columns are moved to the front before the drop set is
  # chosen.
  for (cs in list(list(se_type = "HC2", clusters = NULL),
                  list(se_type = "CR2", clusters = d$cl))) {
    expect_message(
      full <- iv_robust(y ~ en + x1 | x1 + dup, data = d,
                        clusters = cs$clusters, se_type = cs$se_type),
      "returned as NA: en"
    )
    reduced <- lm_robust(y ~ x1, data = d, clusters = cs$clusters,
                         se_type = cs$se_type)
    expect_reduces_to(full, reduced, paste("copied instrument", cs$se_type))
  }

  # A first stage that is flat rather than short of instruments by count.
  b <- flat_iv_data()
  expect_message(full <- iv_robust(Y ~ X | Z, data = b), "returned as NA: X")
  expect_reduces_to(full, lm_robust(Y ~ 1, data = b), "flat first stage")

  expect_message(full <- iv_robust(Y ~ X + W | Z + W, data = b), "returned as NA: X")
  expect_reduces_to(full, lm_robust(Y ~ W, data = b), "flat first stage, exogenous W")
})

test_that("AER::ivreg agrees where it drops the unidentified regressor", {
  skip_if_not_installed("AER")
  expect_message(full <- iv_robust(mpg ~ hp + cyl | am, data = mtcars),
                 "collinear with other regressors")
  aer <- suppressWarnings(AER::ivreg(mpg ~ hp + cyl | am, data = mtcars))
  expect_equal(coef(full), coef(aer), tolerance = DEG_TOL)

  b <- flat_iv_data()
  expect_message(full <- iv_robust(Y ~ X + W | Z + W, data = b), "returned as NA: X")
  expect_equal(coef(full), coef(AER::ivreg(Y ~ X + W | Z + W, data = b)),
               tolerance = DEG_TOL)

  # Where AER drops by column position it keeps the unidentified regressor
  # instead, and this is the one case the two deliberately disagree. `en`'s
  # fitted values lie in the span of the intercept and `x1` to 6e-16, so AER's
  # `en` coefficient is that combination wearing `en`'s name, and its residuals
  # are formed from `en` rather than from what was fitted.
  aer <- AER::ivreg(y ~ en + x1 | x1 + dup, data = d)
  expect_true(is.na(coef(aer)[["x1"]]))
  expect_false(is.na(coef(aer)[["en"]]))
  expect_message(full <- iv_robust(y ~ en + x1 | x1 + dup, data = d),
                 "returned as NA: en")
  expect_true(is.na(coef(full)[["en"]]))
  expect_equal(coef(full)[["x1"]], coef(lm_robust(y ~ x1, data = d))[["x1"]],
               tolerance = DEG_TOL)
})

test_that("the F statistic survives a dropped column that is not the last one", {
  # `coefficients` carries an NA where a column went, while the F statistic's
  # indices and its variance matrix both count positions in the kept set. Where
  # the dropped column was not last the two disagreed, the statistic read the
  # NA, and every robust se_type reported F as NA where lm() and the classical
  # branch returned it. Which of `x1` and `dup2` the pivot takes does not
  # matter: either way the fit spans the reduced model's columns, so it is the
  # reduced model's F statistic.
  for (st in c("HC2", "HC1", "HC3", "stata", "classical")) {
    reduced <- lm_robust(y ~ x1 + x2, data = d, se_type = st)
    for (form in list(y ~ x1 + dup2 + x2, y ~ x1 + x2 + dup2)) {
      expect_message(full <- lm_robust(form, data = d, se_type = st),
                     "collinear with other regressors")
      expect_equal(summary(full)$fstatistic, summary(reduced)$fstatistic,
                   tolerance = DEG_TOL,
                   label = paste(st, deparse(form), "F statistic"))
    }
  }
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

test_that("Horvitz-Thompson refuses a unit observed where its probability was 0", {
  # The inverse-probability weight is infinite, and the estimate came back as
  # NaN or Inf with no message.
  expect_error(horvitz_thompson(y ~ z, data = d, condition_prs = c("0" = 0, "1" = 1)),
               "probability 0 of condition 0 to 40 unit\\(s\\) observed in it")
  expect_error(horvitz_thompson(y ~ z, data = d, condition_prs = c("0" = 1, "1" = 0)),
               "probability 0 of condition 1 to 40 unit\\(s\\) observed in it")
  # Per-unit probabilities: one treated unit given no chance of treatment.
  per_unit <- rep(0.5, n)
  per_unit[which(d$z == 1)[1]] <- 0
  expect_error(horvitz_thompson(y ~ z, data = d, condition_prs = per_unit),
               "probability 0 of condition 1 to 1 unit\\(s\\)")
  # A probability of 0 for a condition no unit is observed in is not a problem
  # for the estimator, and a probability strictly inside (0, 1) is the ordinary
  # case.
  expect_no_error(horvitz_thompson(y ~ z, data = d, condition_prs = c("0" = 0.5, "1" = 0.5)))
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
  messages <- character()
  fit <- withCallingHandlers(
    lm_robust(y ~ x1 + x2, data = tiny, se_type = "classical"),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_equal(sum(is.na(fit$coefficients)), sum(is.na(coef(lm(y ~ x1 + x2, data = tiny)))))
  expect_true(any(grepl("collinear", messages)))
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

# ---- rank deficiency in lm_robust's own design ----

# Built for this section; each test makes its own deficient copy.
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

# The normalization must not cost the rank deficiency the 1e-7 threshold was
# added to catch (#351, #395): an exactly collinear column can otherwise
# survive as a pivot of order 1e-14 and produce a coefficient of order 1e11
# where lm() gives NA.
test_that("exactly collinear columns are still dropped", {
  const <- data.frame(x = rep(1, n), z = rnorm(n), y = rnorm(n))
  expect_message(lm_robust(y ~ x + z, data = const), "collinear")
  fit_const <- suppressWarnings(lm_robust(y ~ x + z, data = const))
  expect_equal(sum(is.na(coef(fit_const))), 1L)
  expect_equal(sum(is.na(coef(fit_const))),
               sum(is.na(coef(lm(y ~ x + z, data = const)))))

  dup <- dat
  dup$x2 <- dup$x1
  expect_message(lm_robust(y ~ x1 + x2 + x3, data = dup), "collinear")
  fit_dup <- suppressWarnings(lm_robust(y ~ x1 + x2 + x3, data = dup))
  expect_equal(sum(is.na(coef(fit_dup))), 1L)
  expect_equal(sum(is.na(coef(fit_dup))),
               sum(is.na(coef(lm(y ~ x1 + x2 + x3, data = dup)))))
})

test_that("a rescaled collinear column is still dropped", {
  dup <- dat
  dup$x2 <- dup$x1 * 1e9
  expect_message(lm_robust(y ~ x1 + x2 + x3, data = dup), "collinear")
  fit <- suppressWarnings(lm_robust(y ~ x1 + x2 + x3, data = dup))
  expect_equal(sum(is.na(coef(fit))), 1L)
  expect_equal(sum(is.na(coef(fit))),
               sum(is.na(coef(lm(y ~ x1 + x2 + x3, data = dup)))))
})

test_that("an all-zero column does not divide by its norm", {
  zeroed <- dat
  zeroed$x2 <- 0
  expect_message(lm_robust(y ~ x1 + x2 + x3, data = zeroed), "collinear")
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
test_that("#351: a constant regressor is detected as collinear, as in lm()", {
  d <- data.frame(y = rnorm(500), x = 1)
  m <- suppressWarnings(lm_robust(y ~ x, data = d))
  expect_true(is.na(m$coefficients[["x"]]))
  expect_equal(unname(is.na(coef(lm(y ~ x, data = d)))), unname(is.na(m$coefficients)))
})

# ---- CR2 with a fixed-effects design short of full rank by two or more ----

test_that("absorbed CR2 matches the dummy expansion when two FE columns drop", {
  # getMeatXtX() compacted the design by removing tossed columns in QR pivot
  # order, which is only correct in descending order: the first left shift
  # renumbers every column after it, so the second removal took the wrong one.
  # One redundant column is exact either way, which is why the nested and
  # disconnected two-factor probes all agreed and this stayed hidden. B3 below
  # is a coarsening of A, so the FE design is short by two.
  skip_if_not_installed("clubSandwich")
  set.seed(3); n <- 400
  d <- data.frame(A = sample(1:12, n, TRUE), C = sample(1:6, n, TRUE))
  d$B3 <- ((d$A - 1) %/% 3) + 1
  d$cl <- sample(1:25, n, TRUE)
  d$x <- rnorm(n)
  d$y <- 0.5 * d$x + d$A * 0.1 + d$C * 0.2 + rnorm(n)

  absorbed <- lm_robust(y ~ x, fixed_effects = ~ A + B3 + C, data = d,
                        clusters = cl, se_type = "CR2")
  dummies <- suppressWarnings(
    lm_robust(y ~ x + factor(A) + factor(B3) + factor(C), data = d,
              clusters = cl, se_type = "CR2")
  )
  # three columns are dropped, so the bug had something to reorder
  expect_gte(sum(is.na(dummies$coefficients)), 2L)

  cs <- clubSandwich::vcovCR(
    lm(y ~ x + factor(A) + factor(B3) + factor(C), data = d),
    cluster = d$cl, type = "CR2"
  )
  # tight, not testthat's default 1.5e-8: these are three routes to one number
  # computed in one session, and the gap the bug left was 3.1%. 1e-9 rather
  # than tighter because this is set on macOS and runs on the CI matrix; the
  # measured gap here is 4.5e-11.
  expect_equal(absorbed$std.error[["x"]], dummies$std.error[["x"]], tolerance = 1e-9)
  expect_equal(absorbed$std.error[["x"]], sqrt(cs["x", "x"]), tolerance = 1e-9)
  expect_equal(absorbed$df[["x"]], dummies$df[["x"]], tolerance = 1e-9)
})

# ---- leverage at or near 1 (estimatr #395) ----

test_that("#395: NaN standard errors from leverage-1 points are explained", {
  set.seed(7)
  N <- 50
  d <- data.frame(x = sample(1:40, N, TRUE), Z = sample(0:1, N, TRUE))
  d$Y <- 0.1 * d$Z + d$x + rnorm(N)
  # This design is rank deficient AND near-saturated, so several warnings fire
  # together: collinearity, leverage, and on some platforms a negative variance
  # diagonal. Collect them instead of nesting expect_message(), which pins the
  # count as well as the content and so breaks whenever another one is added --
  # as it did on all five CI platforms and on none locally.
  ws <- character(0)
  ms <- character(0)
  withCallingHandlers(
    lm_lin(Y ~ Z, covariates = ~ as.factor(x), data = d),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      ms <<- c(ms, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("leverage", ws)))
  expect_true(any(grepl("collinear", ms)))
  expect_false(any(grepl("NaNs produced", ws, fixed = TRUE)))
  # classical SEs do not use leverage, so they stay finite
  m <- suppressWarnings(lm_lin(Y ~ Z, covariates = ~ as.factor(x), data = d,
                               se_type = "classical"))
  expect_false(is.nan(m$std.error[["Z"]]))
})

# The tests below pin what actually goes wrong at a high-leverage point,
# because the #395 test above pins only that a warning fires and the NEWS entry
# for it described a different failure than the one that occurs.
#
# A hat value is a projection diagonal and cannot exceed 1. Leverage of exactly
# 1 costs nothing in the answer: that observation has residual exactly 0, its
# meat contribution is a 0/0 the code resolves to 0, and HC2 returns a finite
# standard error (first test). The trouble is leverage computed marginally ABOVE
# 1, where 1 - h is a small negative number. HC2 then divides by it, half_meat
# takes the square root of the negative result, and every standard error in the
# fit is NaN however small the offending term. HC3 squares the denominator,
# which cancels the sign, so it returned a finite number carrying a spurious
# positive term and said nothing -- the worse of the two failures, because it
# looks like an answer.
#
# lm_variance() now sets the denominator to 0 wherever 1 - h <= 0, which sends
# both through the same isfinite trap that leverage of exactly 1 already used,
# and counts how many observations sit at or near leverage 1 so lm_robust() can
# warn on the condition rather than on a NaN. The count is tolerant, at
# sandwich's `h > 1 - sqrt(eps)`, while the clamp is the strict `1 - h <= 0`,
# because the two answer different questions. The clamp decides the number, and
# on an exactly saturated design the number is the same whichever side of 1 the
# rounding lands on, so a strict count would leave the warning to an ulp. The
# first test below pins the design where that is so.

test_that("#395: leverage of exactly 1 answers finitely, and says that it did", {
  # The single "c" observation is alone in its cell, so it is fitted exactly:
  # leverage is exactly 1 and the residual is exactly 0, and its contribution to
  # the meat is a 0/0 that is correctly resolved to 0 rather than to NaN. The
  # standard error is therefore built from 8 rows where the fit used 9, which is
  # what the warning reports. sandwich warns on this design too and on the same
  # grounds, returning NaN where estimatr drops the row and answers.
  d <- data.frame(
    g = factor(c("a", "a", "a", "a", "b", "b", "b", "b", "c")),
    Z = c(0, 1, 0, 1, 0, 1, 0, 1, 1),
    Y = c(1.0, 2.0, 1.5, 2.5, 3.0, 4.0, 3.5, 4.5, 9.0)
  )
  X <- model.matrix(Y ~ Z + g, d)
  qrx <- qr(X)
  h <- rowSums(qr.Q(qrx)^2)
  e <- as.vector(d$Y - X %*% qr.coef(qrx, d$Y))
  expect_equal(max(h), 1)          # exactly, not approximately
  expect_equal(min(1 - h), 0)
  expect_equal(e[9], 0)

  # The fit's own solver puts that same hat value at 1 + 2.2e-16 rather than at
  # the exact 1 qr() returns above, so this design is the one that decides
  # whether the warning is tolerant or is left to an ulp. Both halves are
  # pinned: the warning fires and names one observation, and the standard error
  # does not depend on which side of 1 the rounding landed.
  expect_warning(
    m <- lm_robust(Y ~ Z + g, data = d, se_type = "HC2"),
    "1 observation has a computed leverage at or near 1"
  )
  expect_false(is.nan(m$std.error[["Z"]]))
  expect_true(is.finite(m$std.error[["Z"]]))

  # Computed here rather than recorded. With observation 9 contributing 0 the
  # meat is the other 8 rows, and the HC2 standard error follows in closed form
  # from a bread and a hat vector that never go near estimatr's solver.
  M <- summary(lm(Y ~ Z + g, data = d))$cov.unscaled
  denom <- 1 - rowSums((X %*% M) * X)
  omega <- ifelse(denom <= 0, 0, e^2 / denom)
  se_hand <- sqrt(diag(M %*% (t(X) %*% (X * omega)) %*% M))
  expect_equal(m$std.error[["Z"]], se_hand[["Z"]])
})

test_that("#395: the leverage guard drops exactly the 1 - h < 0 observations", {
  # Whether a real fit puts a hat value above 1 is a rounding accident and
  # differs by platform, so the guard itself is pinned at the level of the
  # function that implements it, where the leverage is chosen rather than
  # observed. XtX_inv = 0.3 makes h = 0.3 * x^2, so the fourth observation sits
  # at 1.2 and the first three at 0.3.
  X <- matrix(c(1, 1, 1, 2), ncol = 1)
  M <- matrix(0.3, 1, 1)
  ei <- matrix(rep(1, 4), ncol = 1)
  h <- 0.3 * as.vector(X)^2
  expect_equal(h, c(0.3, 0.3, 0.3, 1.2))

  z_variance <- function(se_type) {
    estimatr:::lm_variance(
      X = X, Xunweighted = NULL, XtX_inv = M, ei = ei, weight_mean = 1,
      cluster = NULL, J = 0L, ci = TRUE, se_type = se_type,
      which_covs = TRUE, fe_rank = 0L, fe_leverage = NULL, n_eff = -1L
    )
  }

  # The guarded answers are the three well-behaved rows and nothing else.
  hc2 <- z_variance("HC2")
  hc3 <- z_variance("HC3")
  expect_equal(hc2[["n_leverage_near_one"]], 1L)
  expect_equal(hc3[["n_leverage_near_one"]], 1L)
  expect_equal(sqrt(hc2[["Vcov_hat"]][1, 1]), sqrt(0.09 * 3 / 0.7))
  expect_equal(sqrt(hc3[["Vcov_hat"]][1, 1]), sqrt(0.09 * 3 / 0.49))

  # What the same inputs gave before the guard, computed here rather than
  # recorded: HC2 was NaN throughout, and HC3 was finite, silent, and four times
  # too large, the offending row supplying 94% of the variance.
  expect_true(suppressWarnings(is.nan(sqrt(sum(0.09 * as.vector(X)^2 / (1 - h))))))
  expect_equal(sqrt(sum(0.09 * as.vector(X)^2 / (1 - h)^2)), 3.0904725, tolerance = 1e-7)
  expect_gt(3.0904725 / sqrt(hc3[["Vcov_hat"]][1, 1]), 4)

  # se_types that never read leverage are untouched and report no count.
  for (ty in c("HC0", "HC1", "classical")) {
    expect_equal(z_variance(ty)[["n_leverage_near_one"]], 0L,
                 label = paste(ty, "leverage count"))
  }
})

# These designs are rank deficient, so a collinearity message fires on every fit
# regardless of se_type. Only the leverage warning is of interest here, so fits
# are run through this rather than through expect_silent().
fit_warnings <- function(expr) {
  ws <- character(0)
  ms <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      ms <<- c(ms, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(fit = val, leverage_warning = any(grepl("leverage", ws, fixed = TRUE)))
}

# Leverage above 1 is a floating-point outcome, so it is asserted rather than
# assumed: if a platform's solver returns max(h) <= 1 for this design, the
# condition under test is absent and the guard has nothing to do.
lev_above_one <- function(fml, data) {
  X <- model.matrix(fml, data)
  fit <- estimatr:::lm_solver(X, matrix(data[[all.vars(fml)[1]]], ncol = 1), TRUE)
  keep <- !is.na(fit$beta_hat)
  Xk <- X[, keep, drop = FALSE]
  h <- rowSums((Xk %*% fit$XtX_inv) * Xk)
  list(any = any(1 - h < 0), max = max(h), h = h)
}

test_that("#395: HC2 and HC3 both warn, and neither returns NaN, above leverage 1", {
  set.seed(7)
  N <- 50
  d <- data.frame(x = sample(1:40, N, TRUE), Z = sample(0:1, N, TRUE))
  d$Y <- 0.1 * d$Z + d$x + rnorm(N)
  fml <- Y ~ Z * as.factor(x)

  lev <- lev_above_one(fml, d)
  skip_if_not(lev$any, "this platform's solver does not put any leverage above 1")
  expect_gt(lev$max, 1)

  # Before the guard HC2 was NaN here and HC3 was finite and silent. The
  # asymmetry between them was the defect; both now warn and both answer.
  for (ty in c("HC2", "HC3")) {
    m <- fit_warnings(lm_robust(fml, data = d, se_type = ty))
    expect_true(m$leverage_warning, label = paste(ty, "leverage warning"))
    expect_false(is.nan(m$fit$std.error[["Z"]]), label = paste(ty, "std.error is NaN"))
  }

  # The se_types that never touch leverage neither warn nor change.
  for (ty in c("HC0", "HC1", "classical")) {
    mm <- fit_warnings(lm_robust(fml, data = d, se_type = ty))
    expect_false(mm$leverage_warning, label = paste(ty, "leverage warning"))
  }
})

test_that("#395: the leverage warning counts the observations it dropped", {
  set.seed(7)
  N <- 50
  d <- data.frame(x = sample(1:40, N, TRUE), Z = sample(0:1, N, TRUE))
  d$Y <- 0.1 * d$Z + d$x + rnorm(N)
  fml <- Y ~ Z * as.factor(x)
  skip_if_not(lev_above_one(fml, d)$any, "no leverage above 1 on this platform")

  # Nesting expect_message() would pin the number of warnings as well as their
  # content, which is what broke the first #395 test on all five CI platforms
  # and on none locally. Collect them and read the leverage one out.
  ws <- character(0)
  ms <- character(0)
  withCallingHandlers(
    lm_robust(fml, data = d, se_type = "HC2"),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      ms <<- c(ms, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("collinear", ms)))
  lev <- grep("computed leverage at or near 1", ws, value = TRUE)
  expect_length(lev, 1)

  # The number in the message is the count the C++ made, and the reference is
  # the same criterion applied to a hat vector computed outside it. This design
  # has more than one such observation, so the plural branch of the message is
  # exercised here and the singular branch in the exactly saturated test above.
  expect_match(lev, "^[0-9]+ observations have ")
  expect_equal(as.integer(sub(" .*$", "", lev)),
               sum(lev_above_one(fml, d)$h > 1 - sqrt(.Machine$double.eps)))
})

test_that("#395: CR2 has no analogous hole -- degeneracy goes through its clamp", {
  # CR2 never forms 1 - h_ii. It eigendecomposes
  # (I - H) - H' + Xo MUWTWUM Xo' per cluster and keeps 1/sqrt(eigenvalue) only
  # above 1e-12, so a degenerate cluster contributes 0 rather than a negative or
  # a sign-flipped term. Cluster 11 below is a singleton carrying its own factor
  # level, which is fitted exactly.
  set.seed(11)
  d <- data.frame(g = factor(c(rep("a", 10), rep("b", 10), "c")),
                  Z = c(rep(0:1, 5), rep(0:1, 5), 1))
  d$Y <- rnorm(21) + as.numeric(d$g)
  d$cl <- c(rep(1:5, 2), rep(6:10, 2), 11)

  X <- model.matrix(Y ~ Z + g, d)
  M <- solve(crossprod(X))
  MUWTWUM <- M %*% crossprod(X) %*% M
  block_eigen <- function(cl) {
    ix <- which(d$cl == cl)
    Xb <- X[ix, , drop = FALSE]
    H <- Xb %*% M %*% t(Xb)
    A <- (diag(length(ix)) - H) - t(H) + Xb %*% MUWTWUM %*% t(Xb)
    min(eigen(A, symmetric = TRUE)$values)
  }
  expect_lte(block_eigen(11), 1e-12)     # the clamp fires
  expect_gt(min(vapply(1:10, block_eigen, numeric(1))), 1e-12)

  m <- lm_robust(Y ~ Z + g, data = d, clusters = cl, se_type = "CR2")
  expect_true(is.finite(m$std.error[["Z"]]))
  expect_gt(m$std.error[["Z"]], 0)
})

test_that("absorbed group effects survive a dropped regressor", {
  # Named on its own rather than left to the surface walk, so a failure says
  # which field. `absorbed_group_effects()` formed each group effect as the
  # fitted value minus the row's regressors times the coefficient vector, and
  # that vector carries an NA wherever a column went, so every group came back
  # NA on any rank-deficient fit while the coefficients themselves were right.
  # Present in the released 2.0.0, and unlike the F statistic it does not
  # depend on where the dropped column sits.
  reduced <- lm_robust(y ~ x1 + x2, data = d, fixed_effects = ~ g, se_type = "HC1")
  for (pos in list(y ~ x1 + g_level + x2, y ~ x1 + x2 + g_level)) {
    expect_message(full <- lm_robust(pos, data = d, fixed_effects = ~ g,
                                     se_type = "HC1"),
                   "returned as NA: g_level")
    expect_false(anyNA(full$fixed_effects))
    expect_equal(full$fixed_effects, reduced$fixed_effects, tolerance = DEG_TOL)
  }
})

test_that("#395: a singleton dummy's own standard error is NA, and only its own", {
  # The other half of the #395 clamp. Zeroing a leverage-1 observation's meat
  # term is free for every coefficient that observation carries no information
  # about, and Frisch-Waugh-Lovell says a singleton dummy's neighbours are
  # exactly that. Its own coefficient is the exception: its whole row and
  # column of the meat is zeroed, so 2.0.0 assembled that variance out of rows
  # that say nothing about it and returned a third of the classical standard
  # error, with a p-value and an interval, on a real fit in the reproduction
  # corpus. `sandwich` returns NaN for the entire fit here and 1.0.6 did too,
  # so the number was new in 2.0.0 and wrong; the eighteen it made estimable
  # alongside it were right.
  set.seed(343)
  n <- 60
  d <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    g = c("solo", sample(c("a", "b", "c"), n - 1, replace = TRUE))
  )

  expect_warning(fit <- lm_robust(y ~ x + g, data = d), "NA for gsolo")

  # The singleton alone loses its standard error, and keeps its estimate.
  expect_true(is.na(fit$std.error[["gsolo"]]))
  expect_true(is.na(fit$p.value[["gsolo"]]))
  expect_true(is.na(fit$conf.low[["gsolo"]]))
  expect_false(is.na(fit$coefficients[["gsolo"]]))
  expect_equal(sum(is.na(fit$std.error)), 1)

  # Every other standard error is the one the design without the singleton
  # gives, which is the property that makes dropping the observation correct.
  reduced <- lm_robust(y ~ x + g, data = d[-1, ])
  shared <- c("(Intercept)", "x", "gb", "gc")
  expect_equal(fit$std.error[shared], reduced$std.error[shared], tolerance = DEG_TOL)
  expect_equal(fit$coefficients[shared], reduced$coefficients[shared], tolerance = DEG_TOL)

  # HC3 reaches the clamp by the same route and must answer the same way.
  expect_warning(fit3 <- lm_robust(y ~ x + g, data = d, se_type = "HC3"),
                 "NA for gsolo")
  expect_true(is.na(fit3$std.error[["gsolo"]]))
  expect_equal(sum(is.na(fit3$std.error)), 1)

  # A design with no singleton is untouched: the guard must not fire on the
  # ordinary case it sits in front of.
  expect_silent(clean <- lm_robust(y ~ x + g, data = d[-1, ]))
  expect_false(anyNA(clean$std.error))
})

test_that("#395: the singleton rule holds under weights", {
  # Weights change both the leverage and the meat, and the clamp is applied to
  # the weighted design, so the rule that the singleton alone loses its
  # standard error has to be shown here rather than inferred from the
  # unweighted fit above.
  set.seed(343)
  n <- 60
  d <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    w = runif(n, 0.5, 2),
    g = c("solo", sample(c("a", "b", "c"), n - 1, replace = TRUE))
  )
  d_reduced <- d[-1, ]

  expect_warning(fit <- lm_robust(y ~ x + g, data = d, weights = w),
                 "NA for gsolo")
  expect_true(is.na(fit$std.error[["gsolo"]]))
  expect_false(is.na(fit$coefficients[["gsolo"]]))
  expect_equal(sum(is.na(fit$std.error)), 1)

  # The neighbours keep the weighted design's standard errors without the
  # singleton row, which is the property that makes discarding it correct.
  reduced <- lm_robust(y ~ x + g, data = d_reduced, weights = w)
  shared <- c("(Intercept)", "x", "gb", "gc")
  expect_equal(fit$std.error[shared], reduced$std.error[shared], tolerance = DEG_TOL)
  expect_equal(fit$coefficients[shared], reduced$coefficients[shared], tolerance = DEG_TOL)

  expect_warning(
    fit3 <- lm_robust(y ~ x + g, data = d, weights = w, se_type = "HC3"),
    "NA for gsolo"
  )
  expect_true(is.na(fit3$std.error[["gsolo"]]))
  expect_equal(sum(is.na(fit3$std.error)), 1)
})

test_that("#395: a multivariate fit carries the rule into every outcome block", {
  # The design is shared across outcomes, so lm_variance computes the
  # non-estimable set once and copies it into each of the ny coefficient
  # blocks. The copy is a loop over m * r + j in the C++ that nothing read
  # until this test, and a fit with two outcomes is the only thing that does.
  set.seed(343)
  n <- 60
  d <- data.frame(
    y = rnorm(n),
    y2 = rnorm(n),
    x = rnorm(n),
    g = c("solo", sample(c("a", "b", "c"), n - 1, replace = TRUE))
  )

  expect_warning(fit <- lm_robust(cbind(y, y2) ~ x + g, data = d),
                 "NA for gsolo")
  expect_equal(dim(fit$std.error), c(5L, 2L))
  expect_true(all(is.na(fit$std.error["gsolo", ])))
  expect_equal(sum(is.na(fit$std.error)), 2)
  expect_false(anyNA(fit$coefficients))

  # Each outcome's column is the univariate fit, NA included.
  expect_warning(fy <- lm_robust(y ~ x + g, data = d), "NA for gsolo")
  expect_warning(fy2 <- lm_robust(y2 ~ x + g, data = d), "NA for gsolo")
  expect_equal(unname(fit$std.error[, "y"]), unname(fy$std.error),
               tolerance = DEG_TOL)
  expect_equal(unname(fit$std.error[, "y2"]), unname(fy2$std.error),
               tolerance = DEG_TOL)

  # The warning names the coefficient once rather than once per outcome, which
  # is what the collinear-drop message already does.
  msg <- tryCatch(lm_robust(cbind(y, y2) ~ x + g, data = d),
                  warning = conditionMessage)
  expect_false(grepl("gsolo, gsolo", msg, fixed = TRUE))
})

test_that("#395: lm_lin's centring makes the intercept non-estimable too", {
  # Which coefficients a full-leverage observation alone identifies is a claim
  # about coefficients as written, so it moves with the parametrisation.
  # lm_lin centres its covariates at the grand mean, and the singleton moves
  # that mean, so the intercept joins the singleton's own centred dummy. The
  # treatment effect, which is what the estimator is for, is untouched.
  set.seed(343)
  n <- 60
  d <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    z = rbinom(n, 1, 0.5),
    g = c("solo", sample(c("a", "b", "c"), n - 1, replace = TRUE))
  )

  expect_warning(
    expect_message(fit <- lm_lin(y ~ z, covariates = ~ x + g, data = d),
                   "collinear"),
    "NA for \\(Intercept\\), gsolo_c"
  )
  expect_true(is.na(fit$std.error[["(Intercept)"]]))
  expect_true(is.na(fit$std.error[["gsolo_c"]]))
  expect_false(is.na(fit$std.error[["z"]]))

  # The singleton sits in one treatment arm, so its interaction is collinear
  # and the column is dropped. That third NA is the collinear rule's, not this
  # one's, and it is the only NA coefficient.
  expect_true(is.na(fit$coefficients[["z:gsolo_c"]]))
  expect_equal(sum(is.na(fit$coefficients)), 1)
  expect_equal(sum(is.na(fit$std.error)), 3)

  # Every coefficient the singleton carries no information about keeps its
  # standard error. The criterion is a share of the classical variance, and on
  # this design it separates 8e-3 and 9e-1 from 1e-31 and below.
  finite <- c("z", "x_c", "gb_c", "gc_c", "z:x_c", "z:gb_c", "z:gc_c")
  expect_false(anyNA(fit$std.error[finite]))
})
