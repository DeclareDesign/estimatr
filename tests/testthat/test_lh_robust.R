library(estimatr)

# lh_robust's own behaviour and the regressions filed against it.
#
# The estimate and standard error of a hypothesis against its own lm_robust
# fit are in test_equivalence.R; its methods in test_methods.R; the refusal of
# a multivariate outcome in test_errors.R.

lh_data <- function() {
  set.seed(42)
  n <- 200
  data.frame(y = rnorm(n), x = rnorm(n), z = rbinom(n, 1, 0.5), cl = rep(1:20, 10))
}

dat <- lh_data()

test_that("#405: lh_robust CIs match lm_robust CIs with clusters", {
  m_cl <- lm_robust(y ~ x + z, data = dat, clusters = cl)
  lh_x <- lh_robust(y ~ x + z, data = dat, clusters = cl, linear_hypothesis = "x=0")
  expect_equal(lh_x$lh$conf.low, unname(confint(m_cl)["x", "2.5 %"]), tolerance = 1e-10)
  expect_equal(lh_x$lh$conf.high, unname(confint(m_cl)["x", "97.5 %"]), tolerance = 1e-10)
})

test_that("#320: lh_robust returns the joint hypothesis", {
  lh2 <- lh_robust(y ~ x + z, data = dat, linear_hypothesis = c("x=0", "z=0"))
  expect_true(all(c("value", "p.value") %in% names(lh2$joint_hypothesis)))
  expect_equal(unname(lh2$joint_hypothesis["numdf"]), 2)
  expect_gte(lh2$joint_hypothesis["p.value"], 0)
  expect_lte(lh2$joint_hypothesis["p.value"], 1)
})

test_that("lh_robust handles an intercept-only model", {
  # estimatr 1.0.6 errors here with "missing value where TRUE/FALSE needed".
  # Its df warning guards with `length(fit$df) > 0 && var(fit$df > 0)`, and
  # var() of a length-one vector is NA, so the `if` has nothing to branch on.
  # The unreleased origin/lh-fixes branch changes that 0 to a 1; this rewrite
  # resolves df per hypothesis by name instead and never calls var(), so the
  # case works rather than being patched. Pinned because a future change to the
  # df logic could reintroduce it silently.
  # Live case: the italian_village_continued design in the ResearchDesigns
  # library, which fits age ~ 1 and cannot run under CRAN estimatr.
  intercept_only <- data.frame(age = dat$y)
  m <- lh_robust(age ~ 1, data = intercept_only,
                 linear_hypothesis = "(Intercept) = 20")
  expect_equal(nrow(m$lh), 1L)
  expect_equal(m$lh$df, unname(m$lm_robust$df["(Intercept)"]))
  expect_false(is.na(m$lh$p.value))
})

test_that("lh_robust survives a single-coefficient fit", {
  # estimatr 1.0.6 errored here. Its degrees-of-freedom check read
  #   if (length(lm_robust_fit$df) > 0 && var(lm_robust_fit$df > 0))
  # with the "> 0" inside var() rather than outside, so it took the variance of
  # a logical vector; and var() of a length-one vector is NA either way, so &&
  # was handed NA and stopped. Every one-coefficient fit hit it. That is the
  # shape Declaration 9.2 of Blair, Coppock and Humphreys (2023) uses, and it
  # kept book.declaredesign.org from rebuilding between 2025-02 and 2026-08.
  # The rewrite resolves df per hypothesis and has no such check; this pins
  # that the one-coefficient path stays reachable.
  set.seed(20260828)
  d <- data.frame(age = rnorm(3, mean = 30, sd = 10))
  m <- lh_robust(age ~ 1, data = d, linear_hypothesis = "(Intercept) = 20")

  td <- tidy(m$lh)
  expect_equal(nrow(td), 1L)
  # the hypothesis value is the estimated intercept minus 20
  expect_equal(td$estimate, mean(d$age) - 20)
  # and the df must come from the fit rather than arriving as NA
  expect_equal(td$df, nrow(d) - 1L)
  expect_true(all(is.finite(c(td$std.error, td$statistic, td$p.value))))
})
