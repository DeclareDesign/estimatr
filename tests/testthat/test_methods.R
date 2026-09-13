library(estimatr)

# Every S3 method on every kind of fit, and what each must agree with.
#
# The methods are thin, and thin code is where a mismatch goes unseen: before
# this file confint() on an lh_robust fit errored on every call, nobs() on one
# returned NULL, tidy(conf.level =) on a multivariate fit returned the fit's own
# 95% intervals under a request for 90%, and predict() on an iv_robust fit with
# fixed effects failed with "non-conformable arguments". The first three were
# in 1.0.6. None was covered (review C6).
#
# The fits are ref_surface_fits() from helper-data.R, the sixteen fit types
# test_return_surface.R pins, plus a Horvitz-Thompson fit. Each method is
# checked against the fields it reports rather than only run, and where a
# method refuses a fit type the refusal and its reason are asserted.
#
# Every comparison is within one fit, so the gaps are rounding: the largest
# measured is 6.6e-16.
METH_TOL <- 1e-10

surface <- ref_surface_data()
fits <- ref_surface_fits(surface)
fits$ht <- horvitz_thompson(y ~ zb, data = surface, condition_prs = c("0" = 0.5, "1" = 0.5))

# The lm_robust part of an lh_robust fit carries the per-coefficient fields.
coefficient_part <- function(fit) if (inherits(fit, "lh_robust")) fit$lm_robust else fit
single_outcome <- function(fit) !is.matrix(coefficient_part(fit)$coefficients)

test_that("every fit prints, and so does its summary", {
  for (nm in names(fits)) {
    expect_no_error(capture.output(print(fits[[nm]])), message = nm)
    s <- summary(fits[[nm]])
    expect_no_error(capture.output(print(s)), message = nm)
  }
})

test_that("summary's coefficient table is the fit's estimates and standard errors", {
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    if (inherits(fit, "lh_robust") || !single_outcome(fit)) next
    table <- summary(fit)$coefficients
    expect_equal(unname(table[, 1]), unname(fit$coefficients), tolerance = METH_TOL, label = nm)
    expect_equal(unname(table[, 2]), unname(fit$std.error), tolerance = METH_TOL, label = nm)
  }
})

test_that("tidy reports the fields of the fit", {
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    tidied <- tidy(fit)
    parts <- if (inherits(fit, "lh_robust")) list(fit$lm_robust, fit$lh) else list(fit)
    field <- function(f) unlist(lapply(parts, function(p) as.vector(p[[f]])))
    expect_equal(tidied$estimate, field("coefficients"), tolerance = METH_TOL, label = nm)
    for (f in c("std.error", "p.value", "conf.low", "conf.high", "df")) {
      expect_equal(tidied[[f]], field(f), tolerance = METH_TOL, label = paste(nm, f))
    }
  }
})

test_that("confint is the fit's intervals, and at another level is estimate plus or minus quantile times SE", {
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    parts <- if (inherits(fit, "lh_robust")) list(fit$lm_robust, fit$lh) else list(fit)
    field <- function(f) unlist(lapply(parts, function(p) as.vector(p[[f]])))

    at_fit_level <- confint(fit)
    expect_equal(unname(at_fit_level[, 1]), field("conf.low"), tolerance = METH_TOL, label = nm)
    expect_equal(unname(at_fit_level[, 2]), field("conf.high"), tolerance = METH_TOL, label = nm)

    # Horvitz-Thompson intervals are normal; everything else uses t on the
    # fit's own degrees of freedom.
    quantile <- if (inherits(fit, "horvitz_thompson")) qnorm(0.95) else qt(0.95, field("df"))
    at_90 <- confint(fit, level = 0.9)
    expect_equal(unname(at_90[, 1]), field("coefficients") - quantile * field("std.error"),
                 tolerance = METH_TOL, label = paste(nm, "lower at 0.9"))
    expect_equal(unname(at_90[, 2]), field("coefficients") + quantile * field("std.error"),
                 tolerance = METH_TOL, label = paste(nm, "upper at 0.9"))
    expect_identical(colnames(at_90), c("5 %", "95 %"))

    tidied <- tidy(fit, conf.level = 0.9)
    expect_equal(tidied$conf.low, unname(at_90[, 1]), tolerance = METH_TOL,
                 label = paste(nm, "tidy(conf.level = 0.9)"))
    expect_equal(tidied$conf.high, unname(at_90[, 2]), tolerance = METH_TOL,
                 label = paste(nm, "tidy(conf.level = 0.9)"))

    first <- rownames(at_fit_level)[1]
    expect_equal(confint(fit, parm = first), at_fit_level[first, , drop = FALSE],
                 label = paste(nm, "parm"))
  }
})

test_that("vcov's diagonal is the squared standard errors", {
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    if (inherits(fit, "lh_robust")) next
    expect_equal(unname(sqrt(diag(vcov(fit)))), as.vector(fit$std.error),
                 tolerance = METH_TOL, label = nm)
  }
})

test_that("nobs is the number of observations fitted, and glance is one row", {
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    expect_equal(nobs(fit), nrow(surface), label = paste(nm, "nobs"))
    if (single_outcome(fit)) {
      expect_equal(nrow(glance(fit)), 1L, label = paste(nm, "glance"))
    }
  }
  expect_error(glance(fits$lmr_mv), "multiple responses")
})

# ---- the methods only regression fits have ----

regression_fits <- fits[vapply(fits, function(f) inherits(f, c("lm_robust", "iv_robust")),
                               logical(1))]

test_that("predict() on the fitting data reproduces the fitted values", {
  for (nm in names(regression_fits)) {
    fit <- regression_fits[[nm]]
    if (!single_outcome(fit)) next
    if (nm == "lmr_fe2") {
      # Two-way absorbed effects are identified only in sum, so no new
      # observation can be given its share; the refusal says so.
      expect_error(predict(fit, newdata = surface), "identified in sum")
      next
    }
    expect_equal(unname(predict(fit, newdata = surface)), unname(fit$fitted.values),
                 tolerance = METH_TOL, label = nm)
  }
})

test_that("predict()'s standard errors and intervals are built from the fit's variance", {
  X <- model.matrix(~ x + z, surface)
  for (nm in c("lmr", "lmr_cl", "lmr_w")) {
    fit <- fits[[nm]]
    by_hand <- sqrt(rowSums((X %*% fit$vcov) * X))
    se_fit <- suppressWarnings(predict(fit, newdata = surface, se.fit = TRUE))$se.fit
    expect_equal(unname(se_fit), unname(by_hand), tolerance = METH_TOL, label = nm)

    q <- qt(0.975, fit$df.residual)
    confidence <- predict(fit, newdata = surface, interval = "confidence")$fit
    expect_equal(unname(confidence[, "upr"] - confidence[, "fit"]), unname(q * by_hand),
                 tolerance = METH_TOL, label = paste(nm, "confidence"))
    prediction <- suppressWarnings(predict(fit, newdata = surface, interval = "prediction"))$fit
    expect_equal(unname(prediction[, "upr"] - prediction[, "fit"]),
                 unname(q * sqrt(by_hand^2 + fit$res_var)),
                 tolerance = METH_TOL, label = paste(nm, "prediction"))
  }
  expect_warning(predict(fits$lmr_w, newdata = surface, interval = "prediction"),
                 "constant prediction variance")
  expect_error(predict(fits$lmr_mv, newdata = surface, se.fit = TRUE), "multivariate outcome")
  expect_error(predict(fits$lmr_fe1, newdata = surface, se.fit = TRUE), "no variance estimate")
})

test_that("variable.names, augment, update and model.frame agree with the fit", {
  # update() re-evaluates the stored call, which names the data `d` because
  # that is ref_surface_fits()'s argument.
  d <- surface
  for (nm in names(regression_fits)) {
    fit <- regression_fits[[nm]]
    coef_names <- if (is.matrix(fit$coefficients)) rownames(fit$coefficients) else names(fit$coefficients)
    expect_identical(unname(variable.names(fit)), unname(coef_names), label = nm)

    if (single_outcome(fit)) {
      expect_equal(unname(augment(fit)$.fitted), unname(fit$fitted.values),
                   tolerance = METH_TOL, label = paste(nm, "augment"))
    }

    # update() on an lm_lin fit is not supported: formula() returns the
    # expanded formula, which lm_lin() refuses.
    if (!is.null(fit$scaled_center)) next
    expect_equal(update(fit, . ~ .)$coefficients, fit$coefficients,
                 tolerance = METH_TOL, label = paste(nm, "update"))
  }
  expect_error(augment(fits$lmr_mv), "multiple outcomes")

  for (nm in c("iv", "iv_diag", "iv_fe")) {
    frame <- model.frame(fits[[nm]])
    expect_equal(nrow(frame), fits[[nm]]$nobs, label = nm)
    expect_true(all(c("y", "z", "iv") %in% names(frame)), label = nm)
  }
})
