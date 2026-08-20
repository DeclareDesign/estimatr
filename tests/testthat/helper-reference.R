# Recorded answers from estimatr 1.0.6.
#
# The comparison tests used to call `estimatr::` live. That worked while this
# package was still named estimatrZero and estimatr 1.0.6 was a separate
# install; under its own name, `Suggests: estimatr` is a dependency on itself,
# so the reference has to be a recording rather than a call.
#
# `data-raw/make_estimatr_reference.R` produces the fixture. It is run by hand
# against an installed estimatr 1.0.6 and is not part of the check.

.reference_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- readRDS(test_path("fixtures", "estimatr_reference.rds"))
    }
    cached
  }
})

# A missing key is an error, never a skip. A skip here would be invisible on a
# green run, and the whole point of the recording is that it cannot go missing
# quietly.
ref <- function(key) {
  values <- .reference_fixture()$values
  if (!key %in% names(values)) {
    stop(
      "no recorded estimatr reference for key '", key, "'.\n",
      "Regenerate the fixture with data-raw/make_estimatr_reference.R ",
      "against an installed estimatr ", .reference_fixture()$estimatr_version, ".",
      call. = FALSE
    )
  }
  values[[key]]
}

# Tolerance for any comparison against the recorded fixture.
#
# The recording is a fixed set of numbers produced on one machine. Comparing
# them to values computed on another machine has a floor set by BLAS and LAPACK
# rather than by this package, so an exact or near-exact comparison is a test
# of the CI runner's linear algebra. That is not hypothetical: the first CI run
# after the freeze failed 11 assertions on Ubuntu and Windows and none on
# macOS, which is where the fixture was recorded, and every failure was too
# small for waldo to render.
#
# 1e-9 is set from the worst case rather than by taste. The tightest quantity
# compared here is the weighted `adj.r.squared`, whose recorded value is
# -7.1e-4, so one ulp of it is a relative difference of 3.1e-13. 1e-9 leaves
# about 3,000 ulps of headroom there and far more everywhere else, while
# remaining several orders of magnitude tighter than any difference that would
# mean the estimators disagree.
#
# Comparisons between two things computed in the same session keep their tight
# tolerances: those are exact by construction and have no platform floor.
REF_TOL <- 1e-9

reference_estimatr_version <- function() {
  .reference_fixture()$estimatr_version
}
