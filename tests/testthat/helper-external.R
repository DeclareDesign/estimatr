# Data and tolerances for the external validation tests.
#
# These tests answer a different question from test_vs_estimatr.R. That file
# pins 2.0 against estimatr 1.0.6, which establishes that the rewrite did not
# change any answer; it cannot establish that 1.0.6 was right. The files here
# compare against implementations with no shared lineage: sandwich and
# clubSandwich live, fixest and plm from a recording, and Stata from output
# frozen in 2019.
#
# The surface is organised by variance estimator rather than by function,
# because the variance estimator is what an external reference can actually
# corroborate and is what a user would be harmed by getting wrong. Each cell
# below exists because some external package computes the same quantity; cells
# where no reference computes the same quantity are covered by the internal
# tests instead (HC2/HC3 with absorbed effects, for instance, has no external
# implementation and is pinned against dummy expansion in test_fe_leverage.R).
#
# Builders live here rather than in helper-data.R because
# data-raw/make_external_reference.R sources this file too. If the recorder and
# the tests built their data separately the recording would drift and nothing
# would report it.

ext_data_ols <- function() {
  set.seed(42)
  n <- 200
  d <- data.frame(
    y  = rnorm(n),
    x  = rnorm(n),
    z  = rbinom(n, 1, 0.5),
    w  = runif(n, 0.5, 2),
    cl = rep(1:20, 10)
  )
  # Unbalanced clusters as well as balanced: CR2 and the Satterthwaite degrees
  # of freedom are both sensitive to cluster size variation, and a balanced
  # design hides errors in the cluster-size bookkeeping.
  d$clu <- sample(13, n, replace = TRUE)
  d
}

ext_data_fe <- function() {
  set.seed(42)
  n <- 300
  d <- data.frame(
    y  = rnorm(n),
    x  = rnorm(n),
    z  = rnorm(n),
    w  = runif(n, 0.5, 2),
    g  = factor(rep(1:20, length.out = n)),
    cl = rep(1:15, each = 20)
  )
  # Drawn rather than cycled, so that g2 is not a deterministic function of g.
  d$g2 <- factor(sample(6, n, replace = TRUE))
  d
}

# The cluster variable is the fixed effect, so every absorbed level is nested
# inside a cluster. This is the only configuration that distinguishes the two
# conventions for counting absorbed parameters in the small-sample correction,
# and estimatr counts them (fixest's `fixef.K = "full"`). On non-nested data
# the two conventions agree and the choice is untested.
ext_data_fe_nested <- function() {
  set.seed(11)
  n <- 300
  data.frame(
    y  = rnorm(n),
    x  = rnorm(n),
    z  = rnorm(n),
    g  = factor(rep(1:15, each = 20)),
    cl = rep(1:15, each = 20),
    tm = rep(1:20, 15)
  )
}

# mtcars, with the weight column the Stata do-files construct. Base R data, so
# the frozen Stata output can never go stale against a changing dataset.
ext_data_stata <- function() {
  d <- mtcars
  d$w <- d$drat / 5
  d
}

# Tolerances.
#
# LIVE_TOL applies where both sides are computed in this session. Those agree
# to around 1e-16 relative in practice; 1e-10 leaves six orders of headroom for
# the two packages' different operation orders on a different BLAS, while
# remaining far tighter than any difference that would mean they disagree.
LIVE_TOL <- 1e-10

# EXT_TOL applies to the recorded fixest and plm values, and has the same
# platform floor as REF_TOL in helper-reference.R: a recording compared on
# another machine is limited by that machine's linear algebra.
EXT_TOL <- 1e-9

# Multi-way absorption is iterative in both packages, so the two answers agree
# only to the convergence tolerance rather than to machine precision. Measured
# across five seeds the relative difference ranged from 1.5e-11 to 6.4e-10;
# 1e-7 is set from that worst case with two orders of headroom. One-way
# absorption is direct and agrees exactly, so it keeps EXT_TOL.
EXT_TOL_ITER <- 1e-7

# Stata's `file write` truncated its output to a handful of significant digits,
# and for most cells that print precision, not either implementation, sets the
# floor on how tightly the comparison can be made. The precision varies by
# several orders of magnitude across the tables: `45.459797` pins eight
# significant digits, while `.00004096` pins four. A single tolerance would
# therefore be far too loose for the first cell or too tight for the second, so
# each comparison derives its own from the number actually written down.
#
# Half an ulp of the last printed digit is the quantisation bound. On top of it
# sits a floor for genuine differences in operation order between this package
# and Stata, which are visible in the well-printed cells: the classical F
# statistic agrees to 1.4e-7 relative where quantisation alone would allow only
# 1.1e-8. 1e-6 covers every such case observed across the three tables with an
# order of headroom.
STATA_FLOOR <- 1e-6

# `printed` is the number as Stata wrote it, kept as a string precisely so that
# the trailing digits are still countable here.
stata_rel_tol <- function(printed) {
  decimals <- ifelse(
    grepl(".", printed, fixed = TRUE),
    nchar(sub("^[^.]*\\.", "", printed)),
    0
  )
  value <- as.numeric(printed)
  half_ulp <- 0.5 * 10^(-decimals)
  half_ulp / abs(value) + STATA_FLOOR
}

# Compare one computed value against one frozen Stata string.
expect_equal_stata <- function(actual, printed, label) {
  expect_equal(
    unname(actual), as.numeric(printed),
    tolerance = stata_rel_tol(printed), label = label
  )
}

.external_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- readRDS(test_path("fixtures", "external_reference.rds"))
    }
    cached
  }
})

# A missing key is an error, never a skip, for the reason given in
# helper-reference.R: a skip is invisible on a green run.
ext_ref <- function(key) {
  values <- .external_fixture()$values
  if (!key %in% names(values)) {
    stop(
      "no recorded external reference for key '", key, "'.\n",
      "Regenerate the fixture with data-raw/make_external_reference.R.",
      call. = FALSE
    )
  }
  values[[key]]
}

external_reference_versions <- function() .external_fixture()$versions

# Read one of the frozen Stata tables. Stata wrote coefficients in its own
# order, which is covariates first and `_cons` last; R puts the intercept
# first. Every caller below names the columns, so the reordering is done at the
# comparison rather than here.
read_stata_fixture <- function(file, col.names) {
  read.table(
    test_path("fixtures", "stata", file),
    col.names = col.names,
    colClasses = "character",
    stringsAsFactors = FALSE
  )
}
