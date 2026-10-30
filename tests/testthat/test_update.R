library(estimatr)

# update() on the fits whose stored formula is not what the estimator reads:
# lm_lin() splits the formula into a treatment and covariates, and lh_robust()
# passes it through its dots. Before 2.1.0 update() on an lm_lin() fit
# expanded `. ~ . + w` against formula(fit), `y ~ z + x`, and lm_lin() refused
# the result. Also here: the lh_robust() forwarders, and the model.frame()
# refusal on difference_in_means() and horvitz_thompson() fits, which
# model.frame.default() answered from the caller's workspace.
UPD_TOL <- 1e-10

d <- ref_surface_data()

test_that("update() on an lm_lin() fit adds and removes covariates as a fresh fit would", {
  fit <- lm_lin(y ~ z, covariates = ~ x, data = d)
  same <- update(fit, . ~ .)
  expect_equal(same$coefficients, fit$coefficients, tolerance = UPD_TOL)

  grown <- update(fit, . ~ . + w)
  fresh <- lm_lin(y ~ z, covariates = ~ x + w, data = d)
  expect_equal(grown$coefficients, fresh$coefficients, tolerance = UPD_TOL)
  expect_equal(grown$std.error, fresh$std.error, tolerance = UPD_TOL)
  expect_equal(formula(grown), formula(fresh))

  swapped <- update(grown, . ~ . - x)
  expect_equal(swapped$coefficients,
               lm_lin(y ~ z, covariates = ~ w, data = d)$coefficients,
               tolerance = UPD_TOL)
})

test_that("update() on an lm_lin() fit refuses to drop the treatment or every covariate", {
  fit <- lm_lin(y ~ z, covariates = ~ x, data = d)
  expect_error(update(fit, . ~ . - z), "cannot remove the treatment `z`")
  expect_error(update(fit, . ~ . - x),
               "cannot remove every covariate.*lm_robust\\(y ~ z\\)")
})

test_that("update() on an lm_lin() fit takes new data and new arguments", {
  fit <- lm_lin(y ~ z, covariates = ~ x, data = d)
  half <- d[seq_len(nrow(d) %/% 2), ]
  expect_equal(update(fit, data = half)$coefficients,
               lm_lin(y ~ z, covariates = ~ x, data = half)$coefficients,
               tolerance = UPD_TOL)
  expect_equal(update(fit, se_type = "HC3")$std.error,
               lm_lin(y ~ z, covariates = ~ x, data = d, se_type = "HC3")$std.error,
               tolerance = UPD_TOL)
})

test_that("update() on an lh_robust() fit refits with the hypothesis", {
  skip_if_not_installed("car")
  fit <- lh_robust(y ~ x + z, data = d, linear_hypothesis = "z = 0")
  grown <- update(fit, . ~ . + w)
  fresh <- lh_robust(y ~ x + z + w, data = d, linear_hypothesis = "z = 0")
  expect_equal(coef(grown), coef(fresh), tolerance = UPD_TOL)
  expect_equal(grown$lh$coefficients, fresh$lh$coefficients, tolerance = UPD_TOL)
  expect_identical(getCall(grown)$linear_hypothesis, "z = 0")
})

test_that("an lh_robust() fit answers the fit-level generics with its lm_robust part", {
  skip_if_not_installed("car")
  skip_if_not_installed("sandwich")
  fit <- lh_robust(y ~ x + z, data = d, linear_hypothesis = "z = 0")
  part <- fit$lm_robust
  expect_identical(coef(fit), coef(part))
  expect_identical(vcov(fit), vcov(part))
  expect_identical(df.residual(fit), df.residual(part))
  expect_identical(predict(fit, newdata = d[1:3, ]), predict(part, newdata = d[1:3, ]))
  expect_identical(model.matrix(fit), model.matrix(part))
  expect_identical(formula(fit), formula(part))
  expect_equal(sandwich::vcovHC(fit, type = "HC2"), vcov(part), tolerance = UPD_TOL)
  expect_error(cooks.distance(fit), "case-deletion diagnostics")

  # The hypothesis table on its own
  expect_named(coef(fit$lh), "z = 0")
  expect_equal(unname(coef(fit$lh)), fit$lh$coefficients, tolerance = UPD_TOL)
  expect_equal(unname(sqrt(diag(vcov(fit$lh)))), fit$lh$std.error, tolerance = UPD_TOL)

  # summary()'s table is the fit's with the hypothesis rows appended
  tbl <- coef(summary(fit))
  expect_identical(rownames(tbl), c(names(coef(part)), "z = 0"))
  expect_equal(unname(tbl[, "Estimate"]), c(unname(coef(part)), fit$lh$coefficients),
               tolerance = UPD_TOL)
})

test_that("model.frame(), terms(), and model.matrix() refuse dim and ht fits with their variables in scope", {
  # model.frame.default() evaluated the formula in the caller's environment:
  # with a `y` and a `zb` in scope it returned a frame built from them.
  y <- d$y
  zb <- d$zb
  dim <- difference_in_means(y ~ zb, data = d)
  ht <- horvitz_thompson(y ~ zb, data = d, condition_prs = c("0" = 0.5, "1" = 0.5))
  for (fit in list(dim, ht)) {
    expect_error(model.frame(fit), "rows used are not stored on the fit")
    expect_error(terms(fit), "not built from a design matrix")
    expect_error(model.matrix(fit), "not built from a design matrix")
  }
})
