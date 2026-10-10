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

# ---- the 2SLS path ----
#
# Both stages have the fixed effects absorbed out of them. The other external
# references reach weighted 2SLS (AER and sandwich) and clustered 2SLS
# (clubSandwich), but none of them absorbs, so until these cells existed the
# 2SLS-with-absorption path had no outside implementation to answer to at all
# (review C4). fixest names the fitted endogenous coefficient `fit_en` and
# estimatr names it `en`, so both vectors are read by position.

test_that("2SLS with absorbed fixed effects matches fixest", {
  d_iv <- ext_data_fe_iv()
  expect_ext_equal(
    iv_robust(y ~ en + x | inst + x, data = d_iv, fixed_effects = ~ g,
              se_type = "classical"),
    "fixest_iv_fe1_iid"
  )
  expect_ext_equal(
    iv_robust(y ~ en + x | inst + x, data = d_iv, fixed_effects = ~ g,
              se_type = "HC1"),
    "fixest_iv_fe1_hetero"
  )
  expect_ext_equal(
    iv_robust(y ~ en + x | inst + x, data = d_iv, fixed_effects = ~ g,
              clusters = cl, se_type = "stata"),
    "fixest_iv_fe1_cluster"
  )
  expect_ext_equal(
    iv_robust(y ~ en + x | inst + x, data = d_iv, fixed_effects = ~ g,
              weights = w, se_type = "HC1"),
    "fixest_iv_fe1_w_hetero"
  )
})

test_that("two-way absorption on the 2SLS path matches fixest", {
  d_iv <- ext_data_fe_iv()
  expect_ext_equal(
    iv_robust(y ~ en + x | inst + x, data = d_iv, fixed_effects = ~ g + cl,
              se_type = "HC1"),
    "fixest_iv_fe2_hetero", tol = EXT_TOL_ITER
  )
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

test_that("the recorded reference is the one these versions produced", {
  v <- external_reference_versions()
  # Pinned rather than checked for existence. Every value in the fixture is
  # frozen, so a regeneration under another release of any of these packages
  # has to arrive as a diff in this file rather than as a reference that has
  # quietly become a different one (review C10).
  expect_equal(v[["fixest"]], "0.14.2")
  expect_equal(v[["plm"]], "2.6.7")
  expect_equal(v[["blkvar"]], "0.0.1.6")
  expect_equal(v[["randomizr"]], "2.0.1")
  # R's version is provenance rather than a convention any of these answers
  # depends on, so it is recorded and not pinned.
  expect_true(nzchar(v[["R"]]))

  # The 1.0.6 recording's own accessor, which was defined and called by
  # nothing. The whole of test_vs_estimatr.R is a comparison against this
  # version, and nothing said which version that was.
  expect_equal(reference_estimatr_version(), "1.0.6")
})
