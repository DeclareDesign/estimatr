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

test_that("tidy, glance and augment return tibbles, as broom's methods do", {
  # 1.x returned plain data frames; every downstream reader indexes by name.
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    expect_s3_class(tidy(fit), "tbl_df")
    if (single_outcome(fit)) expect_s3_class(glance(fit), "tbl_df")
    if (inherits(fit, c("lm_robust", "iv_robust")) && single_outcome(fit)) {
      expect_s3_class(estimatr:::augment.lm_robust(fit), "tbl_df")
      # Two-way absorbed effects refuse predict() with newdata, see below.
      if (nm != "lmr_fe2") {
        expect_s3_class(estimatr:::augment.lm_robust(fit, newdata = surface[1:3, ]), "tbl_df")
      }
    }
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

test_that("vcov is symmetric, update takes a new formula, and tidy names each outcome", {
  expect_equal(vcov(fits$lmr), t(vcov(fits$lmr)), tolerance = 1e-14)

  d <- surface
  grown <- update(fits$lmr, . ~ . + w)
  expect_equal(grown$coefficients, lm_robust(y ~ x + z + w, data = surface)$coefficients,
               tolerance = METH_TOL)

  tidied <- tidy(fits$lmr_mv)
  expect_equal(unique(tidied$outcome), c("y", "y2"))
  expect_equal(tidied$term, rep(c("(Intercept)", "x", "z"), 2))
})

# ---- regressions in particular methods ----

test_that("#123: variable.names() returns the model terms", {
  expect_equal(estimatr:::variable.names.lm_robust(fits$lmr), c("(Intercept)", "x", "z"))
})

test_that("B10: glance() reports the residual df, not the first coefficient's", {
  # x[["df"]] is per-coefficient and under CR2 is Satterthwaite, so the column
  # named df.residual read 8.56 on a 10-cluster fit whose residual df is 98.
  set.seed(1)
  N <- 100
  d <- data.frame(y = rnorm(N), x = rnorm(N), cl = sample(10, N, TRUE))
  m <- lm_robust(y ~ x, clusters = cl, data = d, se_type = "CR2")
  expect_equal(glance(m)$df.residual, m$df.residual)
  expect_false(isTRUE(all.equal(m$df.residual, unname(m$df[[1]]))))
  # and it agrees with the iv_robust method, which always used df.residual
  mi <- iv_robust(y ~ x | x, data = d, clusters = cl, se_type = "CR2")
  expect_equal(glance(mi)$df.residual, mi$df.residual)
})

# augment() is called by its method name below: clubSandwich, loaded by a later
# file, registers methods for these classes too, so bare dispatch would depend
# on file order.

test_that("#377: augment() returns the model frame with .fitted and .resid", {
  a <- estimatr:::augment.lm_robust(fits$lmr)
  expect_true(all(c(".fitted", ".resid") %in% names(a)))
  expect_equal(nrow(a), nrow(surface))
  expect_equal(a$.fitted + a$.resid, surface$y, tolerance = 1e-12)
})

test_that("#377: augment() works for lm_lin, iv_robust and fixed effects", {
  expect_true(".fitted" %in% names(estimatr:::augment.lm_robust(fits$lin)))
  expect_true(".fitted" %in% names(estimatr:::augment.iv_robust(fits$iv)))
  a <- estimatr:::augment.lm_robust(fits$lmr_fe1)
  expect_equal(a$.fitted + a$.resid, surface$y, tolerance = 1e-8)
})

test_that("#377: augment(newdata =) predicts without residuals", {
  a <- estimatr:::augment.lm_robust(fits$lmr, newdata = surface[1:5, ])
  expect_equal(nrow(a), 5L)
  expect_true(".fitted" %in% names(a))
  expect_false(".resid" %in% names(a))
})

# ---- Horvitz-Thompson post-estimation ----

# horvitz_thompson had only print() and tidy(), where its sibling
# difference_in_means had seven methods. confint(), vcov() and glance() errored
# outright; summary() and nobs() fell through to the base defaults and returned
# something that looked like output without being the estimate.

post_dat <- ref_data_post()
ht <- horvitz_thompson(y ~ z, data = post_dat, condition_prs = c("0" = 0.5, "1" = 0.5))

test_that("Horvitz-Thompson vcov and nobs agree with estimatr", {
  e0 <- ref("post_ht")
  expect_equal(as.numeric(vcov(ht)), e0$vcov)
  expect_equal(nobs(ht), e0$nobs)
})

test_that("Horvitz-Thompson confint agrees with estimatr, including a non-default level", {
  e0 <- ref("post_ht")
  expect_equal(unname(confint(ht)), e0$confint)
  # The level argument is the case that exposes the t-versus-z choice: a
  # Horvitz-Thompson fit carries df = NA, so rebuilding the interval off a t
  # quantile returns NA bounds rather than a wrong number.
  expect_equal(unname(confint(ht, level = 0.90)), e0$confint_90)
  expect_false(anyNA(confint(ht, level = 0.90)))
})

test_that("Horvitz-Thompson glance returns estimatr's four columns", {
  e0 <- ref("post_ht")
  g <- generics::glance(ht)
  expect_equal(names(g), c("nobs", "se_type", "condition2", "condition1"))
  expect_equal(names(g), e0$glance_names)
  expect_equal(g$nobs, 40L)
})

test_that("Horvitz-Thompson summary reports a z statistic rather than a t", {
  # The estimator has no degrees of freedom to spend, so the headings must not
  # promise any.
  cols <- colnames(summary(ht)$coefficients)
  expect_true("z value" %in% cols)
  expect_false("t value" %in% cols)
})

test_that("modelsummary gets goodness-of-fit rows for a Horvitz-Thompson fit", {
  skip_if_not_installed("modelsummary")
  # Without glance() this call did not fail. It quietly dropped the GOF rows and
  # printed a coefficient table that looked complete, which is why the missing
  # method was worth finding. modelsummary reads tidy() and glance(), never
  # extract(), so this is the path that matters now.
  out <- modelsummary::modelsummary(ht, output = "data.frame")
  gof <- out$term[out$part == "gof"]
  expect_true("Num.Obs." %in% gof)
  expect_true(all(c("Std.Errors", "condition2", "condition1") %in% gof))
})

# ---- other packages' entry points ----

test_that("commarobust and starprep error with a pointer to the replacement", {
  expect_error(commarobust(lm(y ~ z, data = post_dat)), "removed in estimatr 2\\.0")
  expect_error(commarobust(lm(y ~ z, data = post_dat)), "lm_robust")
  expect_error(starprep(lm(y ~ z, data = post_dat)), "removed in estimatr 2\\.0")
  expect_error(starprep(lm(y ~ z, data = post_dat)), "modelsummary")
})

test_that("extract returns a texreg object for both fit types", {
  skip_if_not_installed("texreg")
  # Exported as plain functions on purpose: texreg looks up extract.<class> by
  # name rather than dispatching on a generic, so S3method() would hide them.
  expect_s4_class(extract.lm_robust(lm_robust(y ~ z + x, data = post_dat)), "texreg")
  expect_s4_class(extract.iv_robust(fits$iv), "texreg")
})

test_that("C3: emmeans works when its namespace is loaded rather than attached", {
  # recover_data.lm_robust called getS3method("recover_data", "lm") with no
  # `envir`, so the lookup searched the caller's path for a generic that lives
  # in emmeans. `emmeans::emmeans(...)` loads the namespace without attaching
  # it, which is the ordinary way to call it, and every such call failed with
  # "no function 'recover_data' could be found" -- surfaced to the user as
  # "Perhaps a 'data' or 'params' argument is needed".
  skip_if_not_installed("emmeans")
  set.seed(1)
  n <- 100
  d <- data.frame(y = rnorm(n), g = factor(sample(3, n, TRUE)))

  em <- as.data.frame(emmeans::emmeans(
    lm_robust(y ~ g, data = d, se_type = "classical"), "g"
  ))
  el <- as.data.frame(emmeans::emmeans(lm(y ~ g, data = d), "g"))
  expect_equal(em$emmean, el$emmean, tolerance = 1e-12)
  # classical standard errors are lm's, so the whole path is checked and not
  # only that it returns something
  expect_equal(em$SE, el$SE, tolerance = 1e-12)

  # and the robust default really does reach emmeans, rather than being
  # silently replaced by lm's
  hc2 <- as.data.frame(emmeans::emmeans(lm_robust(y ~ g, data = d), "g"))
  expect_false(isTRUE(all.equal(hc2$SE, el$SE)))
})

# ---- the rank-deficient surface ----
#
# C3 gave recover_data.lm_robust the `envir` it needed, and the test above
# covers the full-rank path. The branch that made the fix necessary, the
# `pass.it.on` attribute a fit with a dropped column carries, was exercised by
# nothing, and neither was any other method's rank-deficient branch: the NA
# padding in `vcov(complete = TRUE)`, the count `print()` reports, or the
# non-estimable basis emmeans builds from the passed-on design (review C6, C13).
#
# The design has two factors and an empty cell, so `fb:qq` is dropped and the
# a:q cell is not estimable from what was fitted. A collinear column on its own
# would exercise the padding but never the non-estimable basis, because every
# marginal mean would still be estimable.
rank_deficient_data <- function() {
  set.seed(343)
  d <- expand.grid(f = factor(c("a", "b")), q = factor(c("p", "q")), rep = 1:12)
  d <- d[!(d$f == "a" & d$q == "q"), ]
  d$y <- rnorm(nrow(d))
  d$y2 <- rnorm(nrow(d))
  d
}

test_that("vcov pads the dropped coefficient with NA when asked for the whole matrix", {
  d <- rank_deficient_data()
  fit <- suppressMessages(suppressWarnings(lm_robust(y ~ f * q, data = d)))
  expect_lt(fit$rank, fit$k)

  full <- vcov(fit, complete = TRUE)
  expect_equal(dim(full), c(fit$k, fit$k))
  expect_equal(rownames(full), fit$term)
  j <- which(is.na(coef(fit, complete = TRUE)))
  expect_true(all(is.na(full[j, ])))
  expect_true(all(is.na(full[, j])))
  # Every other entry is the matrix the fit holds, in the same order.
  expect_equal(unname(full[-j, -j]), unname(fit$vcov), tolerance = 1e-12)

  # And the incomplete matrix is the estimated coefficients alone.
  expect_equal(dim(vcov(fit, complete = FALSE)), c(fit$rank, fit$rank))
})

test_that("vcov refuses a fit that was asked not to keep its variance", {
  d <- rank_deficient_data()
  fit <- lm_robust(y ~ f, data = d, return_vcov = FALSE)
  expect_error(vcov(fit), "return_vcov = TRUE")
})

test_that("print says how many coefficients are not defined", {
  d <- rank_deficient_data()
  fit <- suppressMessages(suppressWarnings(lm_robust(y ~ f * q, data = d)))
  out <- capture.output(print(summary(fit)))
  expect_true(
    any(grepl("Coefficients: (1 not defined because the design matrix is rank deficient)",
              out, fixed = TRUE))
  )
  # The dropped term keeps its row, all NA, rather than vanishing from the table.
  dropped <- names(which(is.na(coef(fit, complete = TRUE))))
  expect_true(any(grepl(paste0("^", dropped, " +NA"), out)))
})

test_that("an iv fit's summary prints its call and its coefficients", {
  set.seed(343)
  n <- 60
  d <- data.frame(inst = rnorm(n))
  d$en <- d$inst + rnorm(n)
  d$y <- d$en + rnorm(n)
  out <- capture.output(print(summary(iv_robust(y ~ en | inst, data = d))))
  expect_true(any(grepl("iv_robust(formula = y ~ en | inst", out, fixed = TRUE)))
  expect_true(any(grepl("^en ", out)))
})

test_that("emmeans reports a cell the design cannot estimate as non-estimable", {
  skip_if_not_installed("emmeans")
  # The non-estimable basis is built from the design emmeans was handed through
  # `pass.it.on`. Without it the a:q cell would come back as a number built
  # from the dropped column's absence rather than as nonEst.
  d <- rank_deficient_data()
  fit <- suppressMessages(suppressWarnings(lm_robust(y ~ f * q, data = d)))
  em <- as.data.frame(summary(emmeans::emmeans(fit, ~ f * q)))

  missing_cell <- em$f == "a" & em$q == "q"
  expect_equal(sum(missing_cell), 1L)
  expect_true(is.na(em$SE[missing_cell]))
  expect_true(all(!is.na(em$SE[!missing_cell])))

  # The three cells that were observed are the cell means of the data.
  observed <- aggregate(y ~ f + q, data = d, FUN = mean)
  for (i in which(!missing_cell)) {
    target <- observed$y[observed$f == em$f[i] & observed$q == em$q[i]]
    expect_equal(em$emmean[i], target, tolerance = 1e-12,
                 label = paste0("cell ", em$f[i], ":", em$q[i]))
  }
})

test_that("emmeans on a multivariate fit is the per-outcome fits", {
  skip_if_not_installed("emmeans")
  # The n.mult branch builds the basis for every outcome at once with a
  # Kronecker product, and nothing exercised it.
  d <- rank_deficient_data()
  mfit <- lm_robust(cbind(y, y2) ~ f, data = d)
  em <- as.data.frame(summary(emmeans::emmeans(mfit, ~ f | rep.meas)))
  expect_equal(sort(unique(as.character(em$rep.meas))), c("y", "y2"))

  for (outcome in c("y", "y2")) {
    single <- as.data.frame(summary(emmeans::emmeans(
      lm_robust(stats::reformulate("f", response = outcome), data = d), ~ f
    )))
    rows <- em[em$rep.meas == outcome, ]
    rows <- rows[order(rows$f), ]
    single <- single[order(single$f), ]
    expect_equal(rows$emmean, single$emmean, tolerance = 1e-12, label = outcome)
    expect_equal(rows$SE, single$SE, tolerance = 1e-12,
                 label = paste(outcome, "SE"))
  }
})

test_that("C13: update.iv_robust carries extra arguments, both new and already in the call", {
  # Two branches, neither run: an argument already named in the call is
  # replaced in place, and one that is not is appended. The existing tests only
  # gave it a new formula.
  set.seed(343)
  n <- 200
  d <- data.frame(z = rnorm(n), cl = rep(1:20, 10))
  d$en <- d$z + rnorm(n)
  d$y <- d$en + rnorm(n)

  iv <- iv_robust(y ~ en | z, data = d)

  # `data` is already in the call, so it is replaced rather than appended.
  half <- update(iv, data = d[1:100, ])
  expect_equal(half$nobs, 100)
  expect_equal(half$coefficients,
               iv_robust(y ~ en | z, data = d[1:100, ])$coefficients)

  # `se_type` is not in the call, so it is appended to it.
  hc1 <- update(iv, se_type = "HC1")
  expect_equal(hc1$se_type, "HC1")
  expect_equal(hc1$std.error,
               iv_robust(y ~ en | z, data = d, se_type = "HC1")$std.error)

  # Both at once, one of each kind.
  both <- update(iv, data = d[1:100, ], se_type = "HC1")
  expect_equal(both$nobs, 100)
  expect_equal(both$se_type, "HC1")
})

test_that("C13: update.iv_robust(evaluate = FALSE) returns the call rather than the fit", {
  set.seed(343)
  n <- 200
  d <- data.frame(z = rnorm(n))
  d$en <- d$z + rnorm(n)
  d$y <- d$en + rnorm(n)

  iv <- iv_robust(y ~ en | z, data = d)

  cl <- update(iv, se_type = "HC1", evaluate = FALSE)
  expect_true(is.call(cl))
  expect_equal(cl[["se_type"]], "HC1")
  expect_equal(eval(cl)$coefficients, update(iv, se_type = "HC1")$coefficients)
})

test_that("C13: update.iv_robust says what is missing when the fit has no call", {
  set.seed(343)
  n <- 100
  d <- data.frame(z = rnorm(n))
  d$en <- d$z + rnorm(n)
  d$y <- d$en + rnorm(n)

  iv <- iv_robust(y ~ en | z, data = d)
  iv$call <- NULL

  expect_error(update(iv, . ~ .), "need an object with call component")
})

test_that("C13: glance refuses a multivariate iv_robust fit", {
  # glance() is one row, and the branch that says so for iv_robust had never
  # been run: the lm_robust one had.
  set.seed(343)
  n <- 200
  d <- data.frame(z = rnorm(n))
  d$en <- d$z + rnorm(n)
  d$y <- d$en + rnorm(n)
  d$y2 <- d$en + rnorm(n)

  ivm <- iv_robust(cbind(y, y2) ~ en | z, data = d)

  expect_error(glance(ivm), "multiple responses")
  # The single-outcome fit it is built from still glances to one row.
  expect_equal(nrow(glance(iv_robust(y ~ en | z, data = d))), 1L)
})

test_that("C13: extract() reports the F statistic and the cluster count when asked", {
  # Two optional goodness-of-fit rows, neither exercised: `include.fstatistic`
  # is off by default, and `include.nclusts` only fires on a clustered fit.
  set.seed(343)
  n <- 200
  d <- data.frame(x = rnorm(n), z = rnorm(n), cl = rep(1:20, 10))
  d$y <- d$x + d$z + rnorm(n)

  m <- lm_robust(y ~ x + z, data = d, clusters = cl)
  tr <- extract(m, include.fstatistic = TRUE)

  expect_true("F statistic" %in% tr@gof.names)
  expect_true("N Clusters" %in% tr@gof.names)
  expect_equal(tr@gof[tr@gof.names == "F statistic"],
               unname(m$fstatistic[[1]]))
  expect_equal(tr@gof[tr@gof.names == "N Clusters"], m$nclusters)
  # The cluster count is a count, so it prints without decimals.
  expect_false(tr@gof.decimal[tr@gof.names == "N Clusters"])
  expect_true(tr@gof.decimal[tr@gof.names == "F statistic"])

  # An unclustered fit has no cluster row to report.
  expect_false("N Clusters" %in%
                 extract(lm_robust(y ~ x + z, data = d))@gof.names)
})

test_that("C13: augment() refuses data it cannot line the fitted values up against", {
  # `data` is taken at the caller's word rather than checked against the model
  # frame, so the rows have to agree. The case this guards is the ordinary one
  # of handing back the original data after the fit dropped rows for
  # missingness.
  set.seed(343)
  n <- 100
  d <- data.frame(x = rnorm(n))
  d$y <- d$x + rnorm(n)
  d$x[1:10] <- NA

  fit <- lm_robust(y ~ x, data = d)
  expect_equal(fit$nobs, 90)

  expect_error(augment(fit, data = d),
               "fitted values are not available for every row")

  # The model frame it did fit augments without complaint.
  expect_equal(nrow(augment(fit)), 90)
})

test_that("C13: predict() takes supplied weights for the prediction interval", {
  # Two branches: supplying `weights` scales the residual variance by their
  # inverse, and omitting them on a weighted fit warns that the prediction
  # variance is being held constant.
  set.seed(343)
  n <- 200
  d <- data.frame(x = rnorm(n), w = runif(n, 0.5, 2))
  d$y <- d$x + rnorm(n)

  fit <- lm_robust(y ~ x, data = d)
  p <- predict(fit, newdata = d, interval = "prediction", weights = w)
  flat <- predict(fit, newdata = d, interval = "prediction")

  # The point predictions are the same; only the interval width moves.
  expect_equal(p[["fit"]][, "fit"], flat[["fit"]][, "fit"])

  # A unit weighted above 1 gets a narrower interval than the unweighted one,
  # and a unit weighted below 1 a wider one, since pred.var = res_var / w.
  wide <- which.min(d$w)
  narrow <- which.max(d$w)
  width <- function(m, i) m[["fit"]][i, "upr"] - m[["fit"]][i, "lwr"]
  expect_gt(width(p, wide), width(flat, wide))
  expect_lt(width(p, narrow), width(flat, narrow))

  # Written out: the half width is t * sqrt(var_fit + res_var / w).
  tval <- qt(0.025, fit$df.residual, lower.tail = FALSE)
  se_fit <- predict(fit, newdata = d, se.fit = TRUE)[["se.fit"]]
  expect_equal(
    unname(width(p, 1) / 2),
    unname(tval * sqrt(se_fit[1]^2 + fit$res_var / d$w[1]))
  )

  # A weighted fit asked for a prediction interval without weights says so.
  wfit <- lm_robust(y ~ x, data = d, weights = w)
  expect_warning(
    predict(wfit, newdata = d, interval = "prediction"),
    "Assuming constant prediction variance"
  )
})

test_that("C13: predict() refuses an lm_lin design it cannot rebuild by name", {
  # A defensive guard rather than a reachable path: every `newdata` this
  # package can be given rebuilds the design, including a factor treatment, a
  # multi-valued numeric treatment, and a `newdata` missing a level of either.
  # What it protects against is the coefficients and the rebuilt design drifting
  # apart, so the way to fire it is to make them disagree.
  set.seed(343)
  n <- 120
  d <- data.frame(x = rnorm(n), z = rbinom(n, 1, 0.5))
  d$y <- d$x + d$z + rnorm(n)

  fit <- lm_lin(y ~ z, covariates = ~ x, data = d)
  expect_silent(predict(fit, newdata = d))

  # Rename one coefficient, which is what a drift between the two would look
  # like. Without the guard the prediction would be built from the wrong
  # columns and returned as a number.
  drifted <- fit
  names(drifted$coefficients)[names(drifted$coefficients) == "z:x_c"] <- "z:x_centered"

  expect_error(predict(drifted, newdata = d),
               "Cannot rebuild the lm_lin design from `newdata`. Missing: z:x_centered")
})

test_that("C13: emmeans names the outcomes of a multivariate fit that carries no names", {
  # The multivariate basis is built with a Kronecker product and labelled by
  # the coefficient matrix's column names. lm_robust always names them, so the
  # fallback is defensive; what it guarantees is that a fit reaching emmeans
  # without them is labelled by position rather than losing the grouping.
  set.seed(343)
  n <- 200
  d <- data.frame(g = factor(rep(c("p", "q"), n / 2)), x = rnorm(n))
  d$y1 <- rnorm(n)
  d$y2 <- rnorm(n)

  named <- lm_robust(cbind(y1, y2) ~ g, data = d)
  expect_equal(colnames(coef(named)), c("y1", "y2"))

  unnamed <- named
  colnames(unnamed$coefficients) <- NULL

  em_named <- summary(emmeans::emmeans(named, ~ g | rep.meas))
  em_unnamed <- summary(emmeans::emmeans(unnamed, ~ g | rep.meas))

  # The outcomes fall back to 1 and 2, and every estimate is unchanged.
  expect_equal(as.character(unique(em_unnamed$rep.meas)), c("1", "2"))
  expect_equal(em_unnamed$emmean, em_named$emmean)
  expect_equal(em_unnamed$SE, em_named$SE)
})
