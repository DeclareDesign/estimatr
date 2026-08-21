library(estimatr)

# The shape of a fitted object, pinned against estimatr 1.0.6.
#
# Six fields went missing from 2.0 without anything reporting it, and one that
# survived came back wrong. A NAMESPACE diff could not see any of it, since no
# export changed; the estimator comparisons could not either, since they name
# the fields they check; and a grep of all 35 reverse dependencies found
# nothing, since the calls are unchanged. It surfaced because a revdep read
# `felevels` on a fixed-effects fit and got NULL.
#
# A missing field is the bad case rather than the loud one: `fit$felevels`
# returns NULL, so the reader gets a wrong answer instead of an error. Hence a
# test over the whole surface rather than over a chosen list of fields.

surface_d <- ref_surface_data()
surface_fits <- ref_surface_fits(surface_d)

# Compare every recorded field, descending into nested fits, and return how
# many were actually compared so a test cannot pass by checking nothing.
compare_surface <- function(got, recorded, path) {
  n <- 0L
  for (f in setdiff(names(recorded), SURFACE_SKIP)) {
    a <- recorded[[f]]
    b <- got[[f]]
    if (is.null(b)) next
    if (is.list(a) && !is.data.frame(a) && is.list(b)) {
      n <- n + compare_surface(b, a, paste0(path, "$", f))
    } else if (is.numeric(a) || is.character(a)) {
      expect_equal(unname(b), unname(a), tolerance = REF_TOL,
                   label = paste0(path, "$", f))
      n <- n + 1L
    }
  }
  n
}

# Not answers: `call` differs by construction, `terms` and `formula` carry
# environments, `weights` and `outcome` are inputs, and `residuals`,
# `joint_hypothesis` and `treatment_vals` are additions 1.0.6 does not have.
SURFACE_SKIP <- c("call", "terms", "formula", "terms_regressors", "outcome",
                  "weights", "model")

for (nm in names(surface_fits)) {
  test_that(paste0("return fields survive from 1.0.6: ", nm), {
    recorded <- ref(paste0("surface_", nm))
    got <- surface_fits[[nm]]
    expect_equal(setdiff(names(recorded), c(names(got), SURFACE_SKIP)), character(0))
  })

  test_that(paste0("return values agree with 1.0.6: ", nm), {
    recorded <- ref(paste0("surface_", nm))
    got <- surface_fits[[nm]]
    n_checked <- compare_surface(got, recorded, nm)
    # lh_robust's components are themselves fits, so a comparison that only
    # looked one level down would check nothing at all and report a green
    # empty test. Assert that something was actually compared.
    expect_gt(n_checked, 0)
  })
}

# The specific fields the revdep run found missing, named so a failure says
# which one rather than only that the surface shrank.
test_that("fixed-effects fits carry felevels, and it holds level names", {
  # Names rather than the per-factor level counts, which is what the one place
  # that did set this field was putting in it. Order follows the character sort
  # the FE matrix imposes, which is what 1.0.6 returns too.
  expect_equal(surface_fits$lmr_fe1$felevels, ref("surface_lmr_fe1")$felevels)
  expect_type(surface_fits$lmr_fe1$felevels$id, "character")
  expect_length(surface_fits$lmr_fe1$felevels$id, nlevels(surface_d$id))
  expect_named(surface_fits$lmr_fe2$felevels, c("id", "tt"))
  # iv_robust absorbs the same way and had lost the field entirely
  expect_equal(surface_fits$iv_fe$felevels, ref("surface_iv_fe")$felevels)
})

test_that("fixed-effects fits carry the unprojected tss alongside proj_tss", {
  fe1 <- surface_fits$lmr_fe1
  expect_equal(fe1$tss, sum((surface_d$y - mean(surface_d$y))^2))
  expect_false(isTRUE(all.equal(fe1$tss, fe1$proj_tss)))
  expect_equal(surface_fits$iv_fe$tss, sum((surface_d$y - mean(surface_d$y))^2))
})

test_that("proj_fstatistic counts every regressor, since FE absorbed the intercept", {
  # Two regressors under test, so numdf is 2. It was 1: the intercept was
  # subtracted although absorption had already removed it from the design.
  expect_equal(unname(surface_fits$lmr_fe1$proj_fstatistic["numdf"]), 2)
  expect_equal(unname(surface_fits$lmr_fe2$proj_fstatistic["numdf"]), 2)
  # With one regressor the subtraction took numdf to zero and the statistic
  # was dropped rather than returned.
  expect_equal(unname(surface_fits$iv_fe$proj_fstatistic["numdf"]), 1)
  expect_false(is.null(surface_fits$iv_fe$proj_fstatistic))
})

test_that("iv_robust with fixed effects returns the absorbed group effects", {
  fe <- surface_fits$iv_fe$fixed_effects
  expect_length(fe, nlevels(surface_d$id))
  expect_equal(fe, ref("surface_iv_fe")$fixed_effects, tolerance = REF_TOL)
})

test_that("lm_lin carries treatment_levels, including the baseline", {
  expect_equal(surface_fits$lin$treatment_levels, c(0, 1))
  # treatment_vals is the other thing: the non-baseline levels predict()
  # expands against. The two are not interchangeable.
  expect_equal(surface_fits$lin_multi$treatment_levels, c(0, 2, 5))
  expect_equal(unname(surface_fits$lin_multi$treatment_vals), c(2, 5))
})

test_that("an unclustered blocked difference_in_means reports nclusters as 0", {
  # Dropping the field instead makes it indistinguishable from NULL to a reader.
  expect_equal(surface_fits$dim_bl$nclusters, 0)
})

test_that("felevels lists the levels in the fit, not the levels in the data", {
  # Downstream code sizes a degrees-of-freedom correction off
  # `length(fit$felevels[[term]])`, so this is a number rather than a label. A
  # level whose every row was dropped for missingness is not in the fit and
  # must not be counted. eventstudyr does exactly this, and counting the
  # absent levels moved its standard errors in the fourth significant figure.
  set.seed(2)
  n <- 200
  d <- data.frame(y = rnorm(n), x = rnorm(n), grp = rep(1:20, each = 10))
  idvar <- "grp"

  complete <- lm_robust(y ~ x, data = d, fixed_effects = ~ get(idvar), se_type = "HC1")
  expect_named(complete$felevels, "get(idvar)")
  expect_equal(complete$felevels[["get(idvar)"]], as.character(1:20))

  # Drop every row of groups 1 and 2, and part of group 3. Eighteen groups
  # remain in the estimation sample.
  d$y[d$grp %in% 1:2] <- NA
  d$y[d$grp == 3][1:3] <- NA
  dropped <- lm_robust(y ~ x, data = d, fixed_effects = ~ get(idvar), se_type = "HC1")
  expect_equal(dropped$felevels[["get(idvar)"]], as.character(3:20))
  expect_length(dropped$felevels[["get(idvar)"]], 18)
  expect_equal(length(dropped$felevels[["get(idvar)"]]),
               length(unique(d$grp[!is.na(d$y)])))

  # The name survives the row removal. 1.0.6 fell back to `V1` here, which is
  # the one difference from it that this field still carries, and it moves no
  # estimate.
  expect_named(dropped$felevels, "get(idvar)")

  # iv_robust absorbs through the same path
  d$z <- d$x + rnorm(n)
  iv <- iv_robust(y ~ x | z, data = d, fixed_effects = ~ get(idvar), se_type = "HC1")
  expect_named(iv$felevels, "get(idvar)")
  expect_equal(iv$felevels[["get(idvar)"]], as.character(3:20))
})

# ---- names on the n-length fields ----

# `fitted.values` carries the model frame's row names and `residuals` does not,
# which is what 1.0.6 returned in every configuration. The distinction is worth
# pinning rather than leaving to chance: the names arrive only because
# `X %*% beta` inherits the design matrix's rownames, so they appear or vanish
# depending on which branch built the vector. 2.0 briefly named residuals on
# fits without fixed effects and not on fits with them, which made a fit with
# absorbed fixed effects SMALLER than the same fit without them -- 8.0 MB
# against 14.4 MB at n = 100,000, all of it row-name strings.
#
# The comparisons above call unname() on both sides, by design, so they cannot
# see this. Hence a test of its own.

test_that("fitted.values is named and residuals is not, in every configuration", {
  set.seed(1)
  n <- 300
  d <- data.frame(
    y = rnorm(n), y2 = rnorm(n), x = rnorm(n),
    g = factor(sample(6, n, TRUE)), cl = factor(sample(8, n, TRUE)),
    w = runif(n, 1, 2), iv = rnorm(n)
  )
  d$zz <- d$iv + rnorm(n)
  d$yy <- d$zz + rnorm(n)
  d$y[4] <- NA  # a dropped row, so the names are not simply seq_len(n)

  fits <- list(
    plain     = lm_robust(y ~ x, data = d),
    weighted  = lm_robust(y ~ x, weights = w, data = d),
    clustered = lm_robust(y ~ x, clusters = cl, data = d),
    fe        = lm_robust(y ~ x, fixed_effects = ~ g, data = d),
    fe_wt     = lm_robust(y ~ x, fixed_effects = ~ g, weights = w, data = d),
    fe_cl     = suppressWarnings(
      lm_robust(y ~ x, fixed_effects = ~ g, clusters = cl, data = d)),
    iv        = iv_robust(yy ~ zz | iv, data = d),
    iv_fe     = iv_robust(yy ~ zz | iv, fixed_effects = ~ g, data = d)
  )

  for (nm in names(fits)) {
    expect_false(is.null(names(fits[[nm]]$fitted.values)), info = nm)
    expect_null(names(fits[[nm]]$residuals), info = nm)
  }

  # and the same for a multivariate outcome, where the fields are matrices
  mv <- lm_robust(cbind(y, y2) ~ x, data = d)
  mv_fe <- lm_robust(cbind(y, y2) ~ x, fixed_effects = ~ g, data = d)
  for (f in list(mv, mv_fe)) {
    expect_false(is.null(rownames(f$fitted.values)))
    expect_null(rownames(f$residuals))
  }
})
