library(estimatr)

dat <- ref_data_post()
ht <- horvitz_thompson(y ~ z, data = dat, condition_prs = c("0" = 0.5, "1" = 0.5))

# ---- Horvitz-Thompson post-estimation ----

# horvitz_thompson had only print() and tidy(), where its sibling
# difference_in_means had seven methods. confint(), vcov() and glance() errored
# outright; summary() and nobs() fell through to the base defaults and returned
# something that looked like output without being the estimate.

test_that("vcov and nobs agree with estimatr", {
  e0 <- ref("post_ht")
  expect_equal(as.numeric(vcov(ht)), e0$vcov)
  expect_equal(nobs(ht), e0$nobs)
})

test_that("confint agrees with estimatr, including a non-default level", {
  e0 <- ref("post_ht")
  expect_equal(unname(confint(ht)), e0$confint)
  # The level argument is the case that exposes the t-versus-z choice: a
  # Horvitz-Thompson fit carries df = NA, so rebuilding the interval off a t
  # quantile returns NA bounds rather than a wrong number.
  expect_equal(unname(confint(ht, level = 0.90)), e0$confint_90)
  expect_false(anyNA(confint(ht, level = 0.90)))
})

test_that("glance returns estimatr's four columns", {
  e0 <- ref("post_ht")
  g <- generics::glance(ht)
  expect_equal(names(g), c("nobs", "se_type", "condition2", "condition1"))
  expect_equal(names(g), e0$glance_names)
  expect_equal(g$nobs, 40L)
})

test_that("summary reports a z statistic rather than a t", {
  # The estimator has no degrees of freedom to spend, so the headings must not
  # promise any.
  cols <- colnames(summary(ht)$coefficients)
  expect_true("z value" %in% cols)
  expect_false("t value" %in% cols)
})

test_that("modelsummary gets goodness-of-fit rows for an HT fit", {
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

# ---- removed in 2.0 ----

test_that("commarobust and starprep error with a pointer to the replacement", {
  expect_error(commarobust(lm(y ~ z, data = dat)), "removed in estimatr 2\\.0")
  expect_error(commarobust(lm(y ~ z, data = dat)), "lm_robust")
  expect_error(starprep(lm(y ~ z, data = dat)), "removed in estimatr 2\\.0")
  expect_error(starprep(lm(y ~ z, data = dat)), "modelsummary")
})

# ---- texreg ----

test_that("extract returns a texreg object for both fit types", {
  skip_if_not_installed("texreg")
  # Exported as plain functions on purpose: texreg looks up extract.<class> by
  # name rather than dispatching on a generic, so S3method() would hide them.
  expect_s4_class(extract.lm_robust(lm_robust(y ~ z + x, data = dat)), "texreg")
  expect_true(is.function(extract.iv_robust))
})

# ---- lm_lin prediction ----

# predict() rebuilt the lm_lin design by looking the treatment up by term
# label, so a factor `z` produced X[, "z"] when the design matrix holds `zb`
# and `zc`: subscript out of bounds for every treatment that is not a single
# binary column. The design is now rebuilt the way lm_lin() builds it and lined
# up against the coefficients by name.

lin_cases <- ref_lin_cases()
lin_dat <- ref_data_lin()

for (case in lin_cases) {
  test_that(paste0("predict on an lm_lin fit works: ", case$lbl), {
    fit <- lm_lin(reformulate(case$z, "y"), covariates = case$cov, data = lin_dat)
    p <- predict(fit, newdata = lin_dat)
    expect_length(p, nrow(lin_dat))
    expect_false(anyNA(p))
    # Predicting on the fitting data must reproduce the in-sample fit.
    expect_equal(unname(p), unname(fit$fitted.values))
  })
}

for (case in lin_cases) {
  test_that(paste0("predict agrees with estimatr: ", case$lbl), {
    f <- reformulate(case$z, "y")
    fit_z <- lm_lin(f, covariates = case$cov, data = lin_dat)
    p_z <- predict(fit_z, newdata = lin_dat)
    expect_equal(unname(p_z), ref(paste0("post_predict_", case$lbl)))
  })
}

test_that("predict on a subset of newdata uses the levels seen at fit time", {
  sub <- lin_dat[lin_dat$zf != "c", ]
  fit_z <- lm_lin(y ~ zf, covariates = ~ x, data = lin_dat)
  expect_equal(unname(predict(fit_z, newdata = sub)), ref("post_predict_subset"))
})

test_that("a numeric multi-valued treatment keeps its fit-time levels", {
  fit <- lm_lin(y ~ zn, covariates = ~ x, data = lin_dat)
  # Stored so predict() expands against these rather than against whatever
  # values happen to appear in newdata.
  expect_equal(unname(fit$treatment_vals), c(2, 5))
  expect_null(lm_lin(y ~ zf, covariates = ~ x, data = lin_dat)$treatment_vals)
})

test_that("newdata carrying only one treatment level still predicts", {
  # The stored terms keep the fit's factor levels, so the design rebuilds with
  # the unused indicator columns at zero rather than losing them. Worth pinning:
  # this is the case where a design rebuilt from `newdata` alone would silently
  # come out narrower than the coefficient vector.
  nd <- lin_dat
  nd$zf <- factor(rep("a", nrow(nd)), levels = "a")
  fit_z <- lm_lin(y ~ zf, covariates = ~ x, data = lin_dat)
  expect_equal(unname(predict(fit_z, newdata = nd)), ref("post_predict_one_level"))
})


test_that("C3: emmeans works when its namespace is loaded rather than attached", {
  # recover_data.lm_robust called getS3method("recover_data", "lm") with no
  # `envir`, so the lookup searched the caller's path for a generic that lives
  # in emmeans. `emmeans::emmeans(...)` loads the namespace without attaching
  # it, which is the ordinary way to call it, and every such call failed with
  # "no function 'recover_data' could be found" -- surfaced to the user as
  # "Perhaps a 'data' or 'params' argument is needed".
  skip_if_not_installed("emmeans")
  set.seed(1); n <- 100
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
