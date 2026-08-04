library(estimatrZero)

set.seed(1)
dat <- data.frame(y = rnorm(40), z = rep(0:1, 20), x = rnorm(40))
ht <- horvitz_thompson(y ~ z, data = dat, condition_prs = c("0" = 0.5, "1" = 0.5))

# ---- Horvitz-Thompson post-estimation ----

# horvitz_thompson had only print() and tidy(), where its sibling
# difference_in_means had seven methods. confint(), vcov() and glance() errored
# outright; summary() and nobs() fell through to the base defaults and returned
# something that looked like output without being the estimate.

test_that("vcov and nobs agree with estimatr", {
  skip_if_not_installed("estimatr")
  ref <- estimatr::horvitz_thompson(y ~ z, data = dat, condition_prs = 0.5)
  expect_equal(as.numeric(vcov(ht)), as.numeric(vcov(ref)))
  expect_equal(nobs(ht), nobs(ref))
})

test_that("confint agrees with estimatr, including a non-default level", {
  skip_if_not_installed("estimatr")
  ref <- estimatr::horvitz_thompson(y ~ z, data = dat, condition_prs = 0.5)
  expect_equal(unname(confint(ht)), unname(confint(ref)))
  # The level argument is the case that exposes the t-versus-z choice: a
  # Horvitz-Thompson fit carries df = NA, so rebuilding the interval off a t
  # quantile returns NA bounds rather than a wrong number.
  expect_equal(unname(confint(ht, level = 0.90)), unname(confint(ref, level = 0.90)))
  expect_false(anyNA(confint(ht, level = 0.90)))
})

test_that("glance returns estimatr's four columns", {
  skip_if_not_installed("estimatr")
  ref <- estimatr::horvitz_thompson(y ~ z, data = dat, condition_prs = 0.5)
  g <- generics::glance(ht)
  expect_equal(names(g), c("nobs", "se_type", "condition2", "condition1"))
  expect_equal(names(g), names(generics::glance(ref)))
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
