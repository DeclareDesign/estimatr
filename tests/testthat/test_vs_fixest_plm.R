library(estimatr)

# Absorbed fixed effects against fixest and plm.
#
# Both packages absorb fixed effects by machinery that shares no code with this
# one, which makes them the right reference for the absorption itself rather
# than for the variance formula. What they corroborate is that alternating
# projections lands on the same coefficients and the same small-sample
# corrections as two independent implementations.
#
# These values are recorded rather than computed live, unlike sandwich and
# clubSandwich. fixest and plm both change their small-sample defaults between
# releases, and a live comparison would then fail on someone else's release
# note rather than on anything about this package. The recording pins the
# answer those versions gave; data-raw/make_external_reference.R regenerates
# it and records the versions used.

expect_ext_equal <- function(fit, key, tol = EXT_TOL) {
  target <- ext_ref(key)
  expect_equal(unname(fit$coefficients), unname(target$coefficients),
               tolerance = tol, label = paste0(key, ": coefficients"))
  expect_equal(unname(fit$std.error), unname(target$std.error),
               tolerance = tol, label = paste0(key, ": std.error"))
}

d <- ext_data_fe()

test_that("one-way absorption matches fixest", {
  expect_ext_equal(
    lm_robust(y ~ x + z, fixed_effects = ~ g, data = d, se_type = "classical"),
    "fixest_fe1_iid"
  )
  expect_ext_equal(
    lm_robust(y ~ x + z, fixed_effects = ~ g, data = d, se_type = "HC1"),
    "fixest_fe1_hetero"
  )
  expect_ext_equal(
    lm_robust(y ~ x + z, fixed_effects = ~ g, clusters = cl, data = d, se_type = "stata"),
    "fixest_fe1_cluster"
  )
  expect_ext_equal(
    lm_robust(y ~ x + z, fixed_effects = ~ g, weights = w, data = d, se_type = "HC1"),
    "fixest_fe1_w_hetero"
  )
})

# Two-way absorption is iterative in both packages, so they agree only to the
# convergence tolerance rather than to machine precision. EXT_TOL_ITER is set
# from the worst case measured across seeds; see helper-external.R.
test_that("two-way absorption matches fixest", {
  expect_ext_equal(
    lm_robust(y ~ x + z, fixed_effects = ~ g + g2, data = d, se_type = "HC1"),
    "fixest_fe2_hetero", tol = EXT_TOL_ITER
  )
  expect_ext_equal(
    lm_robust(y ~ x + z, fixed_effects = ~ g + g2, clusters = cl, data = d,
              se_type = "stata"),
    "fixest_fe2_cluster", tol = EXT_TOL_ITER
  )
})

# ---- how absorbed parameters are counted ----
#
# When every absorbed level sits inside one cluster, the two conventions for
# the small-sample correction diverge: fixest's default `fixef.K = "nested"`
# drops the nested absorbed parameters from K, and `fixef.K = "full"` keeps
# them. This package keeps them, which is what Stata's areg does.
#
# On data where the fixed effect is not nested in the cluster the two
# conventions agree, so the tests above cannot tell them apart and the choice
# would be untested. That is the whole reason this data set exists.

test_that("absorbed parameters are counted in the cluster correction", {
  dn <- ext_data_fe_nested()
  fit <- lm_robust(y ~ x + z, fixed_effects = ~ g, clusters = cl, data = dn,
                   se_type = "stata")
  expect_ext_equal(fit, "fixest_nested_cluster")

  # And is not the other convention. Without this the test above would pass if
  # both packages silently switched to dropping nested parameters.
  other <- ext_ref("fixest_nested_cluster_K_nested")
  expect_false(isTRUE(all.equal(unname(fit$std.error),
                                unname(other$std.error), tolerance = 1e-6)))
})

test_that("the within estimator matches plm with Arellano's variance", {
  # plm's `method = "arellano"` with `type = "HC0"` and no cluster adjustment is
  # CR0 on the absorbed design, reached through panel machinery rather than
  # through this package's absorption.
  dn <- ext_data_fe_nested()
  expect_ext_equal(
    lm_robust(y ~ x + z, fixed_effects = ~ g, clusters = cl, data = dn, se_type = "CR0"),
    "plm_within_arellano_hc0"
  )
})

test_that("the recorded reference names the versions it came from", {
  v <- external_reference_versions()
  expect_true(all(c("fixest", "plm", "R") %in% names(v)))
  expect_true(all(nzchar(v)))
})
