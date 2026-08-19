# Recorded answers from estimatr 1.0.6.
#
# The comparison tests used to call `estimatr::` live, which was right while
# this package was named estimatrZero and estimatr was a different package on
# CRAN. Once it is renamed, `Suggests: estimatr` is a dependency on itself, so
# the reference has to be a recording rather than a call.
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

reference_estimatr_version <- function() {
  .reference_fixture()$estimatr_version
}
