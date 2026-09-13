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

# A design that is hostile in every way the variance estimators can be while
# staying full rank, for the external comparisons that otherwise run only on
# standard normals. Each feature is one a well-conditioned fixture hides:
#
# - Column norms about eight orders of magnitude apart: income in dollars
#   beside a share below 0.02. Rank detection that compared every pivot with
#   the largest in the matrix, as estimatr did until 8b78aef, dropped the share
#   from this design as collinear.
# - A quadratic in age in raw years, so the regressors are correlated the way
#   they are in applied work.
# - Cluster sizes from 1 to 74, eight of them singletons, in shuffled order.
# - Weights spanning more than three orders of magnitude.
# - One income far in the tail, with leverage 0.9993 in a cluster of three.
# - Heteroskedastic errors, and an endogenous regressor for the 2SLS cells.
#
# The quadratic is in age rather than calendar year on purpose. In raw years the
# scaled condition number is 5.6e4 against 42 here, and at that conditioning
# centring the year, which cannot move the true standard errors, moves
# clubSandwich's weighted CR2 by 2.8e-2, estimatr's by 2.4e-5, and a CR2 written
# from its definition by 4.5e-8. A fixture that hard tests clubSandwich rather
# than this package.
ext_data_hard <- function() {
  set.seed(20260913)
  n <- 400
  sizes <- c(rep(1, 8), rep(2, 6), 3, 3, 4, 5, 6, 8, 10, 13, 17, 22, 28, 36, 45, 58, 74)
  sizes <- c(sizes, n - sum(sizes))
  cl <- rep(seq_along(sizes), sizes)
  cl <- sample(seq_along(sizes))[cl]
  d <- data.frame(
    cl = cl,
    age = runif(n, 18, 90),
    income = exp(rnorm(n, 10.5, 1)),
    share = runif(n, 0, 0.02),
    z = rbinom(n, 1, 0.4),
    w = exp(rnorm(n, 0, 1.2)),
    inst = rnorm(n)
  )
  d <- d[sample(n), ]
  d$age2 <- d$age^2
  d$income[which.max(d$income)] <- 5e7
  d$en <- 0.4 * d$inst + 50 * d$share + rnorm(n)
  noise_sd <- 2000 * (1 + d$z + abs(d$age - 50) / 20)
  d$y <- 5000 + 0.02 * d$income + 800 * d$z + 1e5 * d$share + 30 * (d$age - 50) -
    0.5 * (d$age - 50)^2 + 900 * d$en + rnorm(n) * noise_sd
  rownames(d) <- NULL
  d
}

# CR2 written out from its definition, with the working model an identity
# (Pustejovsky and Tipton 2018, equations 4 and 5): for each cluster g,
# A_g = [(I - H)(I - H)']_gg^(-1/2) with H = X (X'WX)^-1 X'W. The hat matrix
# comes from a QR of the column-scaled, weighted design, so the reference loses
# no more precision than that decomposition must. Full-rank designs only: the
# QR's pivot has to be the identity for R to be read off in column order.
cr2_by_definition <- function(X, y, cluster, w = rep(1, length(y))) {
  s <- sqrt(colSums(X^2))
  Xs <- sweep(X, 2, s, "/")
  q <- qr(Xs * sqrt(w))
  stopifnot(identical(q$pivot, seq_len(ncol(X))))
  R_inv <- backsolve(qr.R(q), diag(ncol(X)))
  bread <- tcrossprod(R_inv)
  e <- as.vector(y - Xs %*% (bread %*% crossprod(Xs, w * y)))
  resid_maker <- diag(length(y)) - Xs %*% bread %*% t(Xs * w)
  meat <- matrix(0, ncol(X), ncol(X))
  for (g in unique(cluster)) {
    i <- which(cluster == g)
    ev <- eigen(tcrossprod(resid_maker[i, , drop = FALSE]), symmetric = TRUE)
    adjustment <- ev$vectors %*% (t(ev$vectors) / sqrt(ev$values))
    u <- crossprod(Xs[i, , drop = FALSE] * w[i], adjustment %*% e[i])
    meat <- meat + tcrossprod(u)
  }
  bread %*% meat %*% bread / tcrossprod(s)
}

# clubSandwich cannot be held to LIVE_TOL on ext_data_hard(), and the shortfall
# is on its side: against cr2_by_definition() estimatr agrees to 9e-12 and
# clubSandwich to 1.4e-6, and unweighted the gap closes to 1e-12 when the one
# high-leverage observation is removed. The 2SLS gap closes the same way, from
# 1.4e-6 to 1.6e-13. The worst measured difference between
# estimatr and clubSandwich is 1.3e-4, on the unweighted Satterthwaite degrees
# of freedom, and 2e-3 keeps ~15x over it. That is loose, and it is the
# definition that holds these cells tight; clubSandwich is here to catch an
# error the definition shares with this package.
HARD_CLUB_TOL <- 2e-3

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
