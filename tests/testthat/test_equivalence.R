library(estimatr)

# Two routes to one number.
#
# Several arguments and estimators reach the same quantity through different
# code: a multivariate outcome and one fit per outcome, iv_robust with every
# regressor its own instrument and lm_robust, weights that are all 1 and no
# weights, a cluster per observation and the heteroskedasticity-consistent
# estimator it reduces to, subset= and the rows already filtered. Each pair
# below must agree, and where a pair is not supposed to agree the difference
# is asserted instead, with the reason.
#
# Agreement between two routes is not evidence either is right; that is what
# the test_vs_*.R files are for. What a pair catches is a route that has
# drifted from its twin, which is how try_cholesky came to skip rank detection
# while the default path did not.
#
# Both sides are computed in one session. Every pair below holds at 1e-12; the
# largest single element measured is 9.9e-13, a fitted value on the Cholesky
# path, and EQ_TOL keeps 100x over that.
EQ_TOL <- 1e-10

eq_data <- function() {
  set.seed(1)
  n <- 240
  d <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    z = rbinom(n, 1, 0.5),
    w = runif(n, 0.5, 2),
    cl = sample(25, n, replace = TRUE),
    g = factor(sample(letters[1:8], n, replace = TRUE)),
    h = factor(sample(5, n, replace = TRUE))
  )
  d$y <- 1 + d$x1 - d$x2 + 0.5 * d$z + rnorm(n) * (1 + abs(d$x1))
  d$y2 <- -2 + 0.3 * d$x1 + d$z + rnorm(n)
  d$zc <- as.integer(d$cl %% 2 == 0)
  # Twelve blocks of twenty with ten treated in each, and matched pairs.
  d$blk <- rep(1:12, each = 20)
  d$zb <- as.integer(ave(runif(n), d$blk, FUN = function(u) rank(u) <= 10))
  d$pair <- rep(seq_len(n / 2), each = 2)
  d$zp <- rep(c(0L, 1L), n / 2)
  d$en <- d$x1 + rnorm(n)
  d$inst <- d$x1 + d$en + rnorm(n)
  d$one <- 1
  d$id <- seq_len(n)
  d
}

d <- eq_data()
n <- nrow(d)

expect_same_inference <- function(a, b, label) {
  expect_equal(unname(a$coefficients), unname(b$coefficients), tolerance = EQ_TOL,
               label = paste(label, "coefficients"))
  expect_equal(unname(a$std.error), unname(b$std.error), tolerance = EQ_TOL,
               label = paste(label, "standard errors"))
  expect_equal(unname(a$df), unname(b$df), tolerance = EQ_TOL, label = paste(label, "df"))
}

# ---- a multivariate outcome is one fit per outcome ----
#
# cbind(y, y2) shares the design and the decomposition across outcomes and
# stacks the variance as a Kronecker product. None of that may change any
# single outcome's answer. Multivariate outcomes with clusters had no coverage
# at all before this file (review C4).

test_that("a multivariate fit equals one fit per outcome, for every variance path", {
  fits <- list(
    HC0 = function(f) lm_robust(f, data = d, se_type = "HC0"),
    HC1 = function(f) lm_robust(f, data = d, se_type = "HC1"),
    HC2 = function(f) lm_robust(f, data = d, se_type = "HC2"),
    HC3 = function(f) lm_robust(f, data = d, se_type = "HC3"),
    classical = function(f) lm_robust(f, data = d, se_type = "classical"),
    CR0 = function(f) lm_robust(f, data = d, clusters = cl, se_type = "CR0"),
    CR2 = function(f) lm_robust(f, data = d, clusters = cl, se_type = "CR2"),
    stata = function(f) lm_robust(f, data = d, clusters = cl, se_type = "stata"),
    weighted_HC2 = function(f) lm_robust(f, data = d, weights = w, se_type = "HC2"),
    weighted_CR2 = function(f) lm_robust(f, data = d, weights = w, clusters = cl, se_type = "CR2"),
    fe_two_way_HC3 = function(f) lm_robust(f, data = d, fixed_effects = ~ g + h, se_type = "HC3"),
    fe_CR2 = function(f) lm_robust(f, data = d, fixed_effects = ~ g, clusters = cl, se_type = "CR2"),
    fe_weighted_CR0 = function(f) lm_robust(f, data = d, weights = w, fixed_effects = ~ g,
                                            clusters = cl, se_type = "CR0")
  )
  for (nm in names(fits)) {
    both <- fits[[nm]](cbind(y, y2) ~ x1 + x2 + z)
    one <- fits[[nm]](y ~ x1 + x2 + z)
    two <- fits[[nm]](y2 ~ x1 + x2 + z)
    for (field in c("coefficients", "std.error", "df", "fitted.values")) {
      expect_equal(unname(both[[field]]), unname(cbind(one[[field]], two[[field]])),
                   tolerance = EQ_TOL, label = paste(nm, field))
    }
    expect_equal(unname(both$residuals), unname(cbind(one$residuals, two$residuals)),
                 tolerance = EQ_TOL, label = paste(nm, "residuals"))
    expect_equal(unname(both$r.squared), c(one$r.squared, two$r.squared),
                 tolerance = EQ_TOL, label = paste(nm, "R-squared"))
  }
})

test_that("multivariate lm_lin and iv_robust equal one fit per outcome", {
  both <- lm_lin(cbind(y, y2) ~ z, covariates = ~ x1, data = d)
  one <- lm_lin(y ~ z, covariates = ~ x1, data = d)
  two <- lm_lin(y2 ~ z, covariates = ~ x1, data = d)
  expect_equal(unname(both$coefficients), unname(cbind(one$coefficients, two$coefficients)),
               tolerance = EQ_TOL)
  expect_equal(unname(both$std.error), unname(cbind(one$std.error, two$std.error)),
               tolerance = EQ_TOL)

  both <- iv_robust(cbind(y, y2) ~ en + x2 | inst + x2, data = d)
  one <- iv_robust(y ~ en + x2 | inst + x2, data = d)
  two <- iv_robust(y2 ~ en + x2 | inst + x2, data = d)
  expect_equal(unname(both$coefficients), unname(cbind(one$coefficients, two$coefficients)),
               tolerance = EQ_TOL)
  expect_equal(unname(both$std.error), unname(cbind(one$std.error, two$std.error)),
               tolerance = EQ_TOL)
})

# ---- arguments that should only withhold output ----

test_that("return_vcov, ci and se_type = 'none' withhold output without changing the rest", {
  for (st in c("HC2", "CR2", "classical")) {
    cluster_ids <- if (st == "CR2") d$cl else NULL
    full <- lm_robust(y ~ x1 + x2 + z, data = d, se_type = st, clusters = cluster_ids)

    no_vcov <- lm_robust(y ~ x1 + x2 + z, data = d, se_type = st, clusters = cluster_ids,
                         return_vcov = FALSE)
    expect_same_inference(no_vcov, full, paste(st, "return_vcov = FALSE"))
    expect_null(no_vcov$vcov)

    no_ci <- lm_robust(y ~ x1 + x2 + z, data = d, se_type = st, clusters = cluster_ids,
                       ci = FALSE)
    expect_equal(no_ci$std.error, full$std.error, tolerance = EQ_TOL)
    expect_true(all(is.na(no_ci$p.value)))
    expect_true(all(is.na(no_ci$conf.low)) && all(is.na(no_ci$conf.high)))

    none <- lm_robust(y ~ x1 + x2 + z, data = d, se_type = "none", clusters = cluster_ids)
    expect_equal(none$coefficients, full$coefficients, tolerance = EQ_TOL)
    expect_true(all(is.na(none$std.error)))
  }
})

test_that("lm_robust_fit on a model matrix equals lm_robust on the formula", {
  X <- model.matrix(~ x1 + x2 + z, d)
  cases <- list(
    HC2 = list(weights = NULL, cluster = NULL, se_type = "HC2"),
    weighted_HC3 = list(weights = d$w, cluster = NULL, se_type = "HC3"),
    classical = list(weights = NULL, cluster = NULL, se_type = "classical"),
    CR2 = list(weights = NULL, cluster = d$cl, se_type = "CR2"),
    weighted_CR2 = list(weights = d$w, cluster = d$cl, se_type = "CR2")
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    direct <- lm_robust_fit(y = d$y, X = X, weights = cs$weights, cluster = cs$cluster,
                            se_type = cs$se_type, has_int = TRUE)
    formula_fit <- lm_robust(y ~ x1 + x2 + z, data = d, weights = cs$weights,
                             clusters = cs$cluster, se_type = cs$se_type)
    expect_same_inference(direct, formula_fit, nm)
    expect_equal(direct$vcov, formula_fit$vcov, tolerance = EQ_TOL, label = paste(nm, "vcov"))
  }
})

# ---- estimators that reduce to one another ----

test_that("iv_robust with every regressor its own instrument is lm_robust", {
  for (st in c("HC0", "HC1", "HC2", "HC3", "classical", "CR0", "CR2", "stata")) {
    cluster_ids <- if (st %in% c("CR0", "CR2", "stata")) d$cl else NULL
    for (weighted in c(FALSE, TRUE)) {
      weights <- if (weighted) d$w else NULL
      label <- paste0(st, if (weighted) ", weighted")
      iv <- iv_robust(y ~ x1 + z | x1 + z, data = d, se_type = st,
                      clusters = cluster_ids, weights = weights)
      ols <- lm_robust(y ~ x1 + z, data = d, se_type = st,
                       clusters = cluster_ids, weights = weights)
      expect_same_inference(iv, ols, label)
      expect_equal(iv$r.squared, ols$r.squared, tolerance = EQ_TOL, label = paste(label, "R-squared"))
    }
  }
  for (st in c("HC1", "HC2", "CR2")) {
    cluster_ids <- if (st == "CR2") d$cl else NULL
    iv <- iv_robust(y ~ x1 + z | x1 + z, data = d, se_type = st, clusters = cluster_ids,
                    fixed_effects = ~ g + h)
    ols <- lm_robust(y ~ x1 + z, data = d, se_type = st, clusters = cluster_ids,
                     fixed_effects = ~ g + h)
    expect_same_inference(iv, ols, paste(st, "two-way fixed effects"))
  }
})

test_that("difference_in_means is the regression it describes, design by design", {
  # Welch: the HC2 standard error, with Welch-Satterthwaite rather than n - 2
  # degrees of freedom, which is what t.test() reports.
  dim <- difference_in_means(y ~ z, data = d)
  ols <- lm_robust(y ~ z, data = d, se_type = "HC2")
  expect_equal(dim$coefficients[["z"]], ols$coefficients[["z"]], tolerance = EQ_TOL)
  expect_equal(dim$std.error[["z"]], ols$std.error[["z"]], tolerance = EQ_TOL)
  welch <- t.test(d$y[d$z == 1], d$y[d$z == 0])
  expect_equal(unname(dim$df), unname(welch$parameter), tolerance = EQ_TOL)
  expect_equal(unname(dim$p.value), welch$p.value, tolerance = EQ_TOL)

  # Clustered and weighted designs are delegated to lm_robust outright, df
  # included.
  pairs <- list(
    clustered = list(difference_in_means(y ~ zc, data = d, clusters = cl),
                     lm_robust(y ~ zc, data = d, clusters = cl, se_type = "CR2")),
    weighted = list(difference_in_means(y ~ z, data = d, weights = w),
                    lm_robust(y ~ z, data = d, weights = w, se_type = "HC2")),
    weighted_clustered = list(difference_in_means(y ~ zc, data = d, weights = w, clusters = cl),
                              lm_robust(y ~ zc, data = d, weights = w, clusters = cl, se_type = "CR2"))
  )
  for (nm in names(pairs)) {
    a <- pairs[[nm]][[1]]
    b <- pairs[[nm]][[2]]
    expect_equal(unname(a$coefficients), unname(b$coefficients[2]), tolerance = EQ_TOL, label = nm)
    expect_equal(unname(a$std.error), unname(b$std.error[2]), tolerance = EQ_TOL, label = nm)
    expect_equal(unname(a$df), unname(b$df[2]), tolerance = EQ_TOL, label = nm)
  }

  # Blocked, with every block big: Lin's estimator with the block indicators as
  # the covariates, estimate and standard error both.
  blocked <- difference_in_means(y ~ zb, data = d, blocks = blk)
  lin <- lm_lin(y ~ zb, covariates = ~ factor(blk), data = d)
  expect_equal(unname(blocked$coefficients), unname(lin$coefficients["zb"]), tolerance = EQ_TOL)
  expect_equal(unname(blocked$std.error), unname(lin$std.error["zb"]), tolerance = EQ_TOL)

  # Matched pairs: the paired t-test on the within-pair differences.
  pairs_fit <- difference_in_means(y ~ zp, data = d, blocks = pair)
  diffs <- d$y[d$zp == 1] - d$y[d$zp == 0]
  expect_equal(unname(pairs_fit$coefficients), mean(diffs), tolerance = EQ_TOL)
  expect_equal(unname(pairs_fit$std.error), sd(diffs) / sqrt(length(diffs)), tolerance = EQ_TOL)
  expect_equal(unname(pairs_fit$df), length(diffs) - 1, tolerance = EQ_TOL)
})

test_that("lh_robust is the linear combination of its own lm_robust fit", {
  combination <- c(0, 2, 0, 1)
  for (st in c("HC2", "CR2", "classical")) {
    cluster_ids <- if (st == "CR2") d$cl else NULL
    lh <- lh_robust(y ~ x1 + x2 + z, data = d, se_type = st, clusters = cluster_ids,
                    linear_hypothesis = "z + 2*x1 = 0.5")
    ols <- lm_robust(y ~ x1 + x2 + z, data = d, se_type = st, clusters = cluster_ids)
    estimate <- sum(combination * ols$coefficients) - 0.5
    se <- sqrt(drop(t(combination) %*% ols$vcov %*% combination))

    expect_same_inference(lh$lm_robust, ols, paste(st, "lm_robust component"))
    expect_equal(unname(lh$lh$coefficients), estimate, tolerance = EQ_TOL, label = st)
    expect_equal(unname(lh$lh$std.error), se, tolerance = EQ_TOL, label = st)
    expect_equal(unname(lh$lh$p.value), 2 * pt(-abs(estimate / se), lh$lh$df),
                 tolerance = EQ_TOL, label = st)
  }
})

# ---- degenerate settings of an argument ----

test_that("weights that are all 1 are no weights", {
  for (st in c("HC0", "HC1", "HC2", "HC3", "classical", "CR0", "CR2", "stata")) {
    cluster_ids <- if (st %in% c("CR0", "CR2", "stata")) d$cl else NULL
    weighted <- lm_robust(y ~ x1 + x2 + z, data = d, se_type = st, clusters = cluster_ids,
                          weights = one)
    unweighted <- lm_robust(y ~ x1 + x2 + z, data = d, se_type = st, clusters = cluster_ids)
    expect_same_inference(weighted, unweighted, st)
    expect_equal(weighted$r.squared, unweighted$r.squared, tolerance = EQ_TOL)
  }
  expect_same_inference(
    lm_robust(y ~ x1 + x2 + z, data = d, fixed_effects = ~ g + h, se_type = "HC3", weights = one),
    lm_robust(y ~ x1 + x2 + z, data = d, fixed_effects = ~ g + h, se_type = "HC3"),
    "two-way fixed effects, HC3"
  )
  expect_same_inference(
    difference_in_means(y ~ zb, data = d, blocks = blk, weights = one),
    difference_in_means(y ~ zb, data = d, blocks = blk),
    "blocked difference in means"
  )

  # The unblocked difference in means is the exception, and says so: any
  # weights route it through lm_robust, so the degrees of freedom are n - 2
  # rather than Welch's, while the estimate and standard error are unchanged.
  expect_message(
    weighted <- difference_in_means(y ~ z, data = d, weights = one),
    "Welch-Satterthwaite approximation will not be used"
  )
  unweighted <- difference_in_means(y ~ z, data = d)
  expect_equal(weighted$coefficients, unweighted$coefficients, tolerance = EQ_TOL)
  expect_equal(weighted$std.error, unweighted$std.error, tolerance = EQ_TOL)
  expect_equal(unname(weighted$df), n - 2)
})

test_that("a cluster per observation is the heteroskedasticity-consistent estimator", {
  for (p in list(c("CR0", "HC0"), c("CR2", "HC2"), c("stata", "HC1"))) {
    clustered <- lm_robust(y ~ x1 + x2 + z, data = d, se_type = p[1], clusters = id)
    robust <- lm_robust(y ~ x1 + x2 + z, data = d, se_type = p[2])
    expect_equal(clustered$vcov, robust$vcov, tolerance = EQ_TOL,
                 label = paste(p[1], "with singleton clusters"))
  }
  for (p in list(c("CR0", "HC0"), c("stata", "HC1"))) {
    clustered <- lm_robust(y ~ x1 + x2 + z, data = d, se_type = p[1], clusters = id, weights = w)
    robust <- lm_robust(y ~ x1 + x2 + z, data = d, se_type = p[2], weights = w)
    expect_equal(clustered$vcov, robust$vcov, tolerance = EQ_TOL,
                 label = paste(p[1], "with singleton clusters, weighted"))
  }

  # Weighted CR2 is built against an identity working model and weighted HC2
  # against precision weights (?lm_robust; pinned against clubSandwich's
  # inverse_var in test_vs_clubsandwich.R), so with a cluster per observation
  # they are two estimators, not one. Measured gap 2.4e-2.
  clustered <- lm_robust(y ~ x1 + x2 + z, data = d, se_type = "CR2", clusters = id, weights = w)
  robust <- lm_robust(y ~ x1 + x2 + z, data = d, se_type = "HC2", weights = w)
  expect_gt(max(abs(clustered$std.error / robust$std.error - 1)), 1e-3)
})

test_that("subset= is the same as filtering the rows first", {
  keep <- d$x2 > -1
  pairs <- list(
    CR2 = list(lm_robust(y ~ x1 + x2 + z, data = d, clusters = cl, subset = x2 > -1),
               lm_robust(y ~ x1 + x2 + z, data = d[keep, ], clusters = cl)),
    weighted_HC2 = list(lm_robust(y ~ x1 + x2 + z, data = d, weights = w, subset = x2 > -1),
                        lm_robust(y ~ x1 + x2 + z, data = d[keep, ], weights = w)),
    fe_two_way = list(lm_robust(y ~ x1 + z, data = d, fixed_effects = ~ g + h, subset = x2 > -1),
                      lm_robust(y ~ x1 + z, data = d[keep, ], fixed_effects = ~ g + h)),
    iv_CR2 = list(iv_robust(y ~ en + x2 | inst + x2, data = d, clusters = cl, subset = x2 > -1),
                  iv_robust(y ~ en + x2 | inst + x2, data = d[keep, ], clusters = cl)),
    lin_weighted = list(lm_lin(y ~ z, covariates = ~ x1 + x2, data = d, weights = w, subset = x2 > -1),
                        lm_lin(y ~ z, covariates = ~ x1 + x2, data = d[keep, ], weights = w)),
    dim_clustered = list(difference_in_means(y ~ zc, data = d, clusters = cl, subset = x2 > -1),
                         difference_in_means(y ~ zc, data = d[keep, ], clusters = cl)),
    dim_blocked = list(difference_in_means(y ~ zb, data = d, blocks = blk, subset = x2 > -1),
                       difference_in_means(y ~ zb, data = d[keep, ], blocks = blk))
  )
  for (nm in names(pairs)) {
    expect_same_inference(pairs[[nm]][[1]], pairs[[nm]][[2]], nm)
  }
})

test_that("a missing value anywhere the fit reads drops exactly that row", {
  dm <- d
  dm$y[c(3, 50)] <- NA
  dm$x1[7] <- NA
  dm$w[9] <- NA
  dm$cl[11] <- NA
  dm$g[13] <- NA
  dm$blk[15] <- NA
  # Each fit, the columns it reads, and the variable whose missingness it warns
  # about because that variable is not in the formula.
  cases <- list(
    CR2 = list(function(dd) lm_robust(y ~ x1 + x2 + z, data = dd, clusters = cl),
               c("y", "x1", "x2", "z", "cl"), "cluster"),
    weighted_HC2 = list(function(dd) lm_robust(y ~ x1 + x2 + z, data = dd, weights = w),
                        c("y", "x1", "x2", "z", "w"), "weights"),
    fe_two_way = list(function(dd) lm_robust(y ~ x1 + z, data = dd, fixed_effects = ~ g + h),
                      c("y", "x1", "z", "g", "h"), "fixed_effects"),
    iv_CR2 = list(function(dd) iv_robust(y ~ en + x2 | inst + x2, data = dd, clusters = cl),
                  c("y", "en", "x2", "inst", "cl"), "cluster"),
    lin_weighted = list(function(dd) lm_lin(y ~ z, covariates = ~ x1 + x2, data = dd, weights = w),
                        c("y", "z", "x1", "x2", "w"), "weights"),
    dim_clustered = list(function(dd) difference_in_means(y ~ zc, data = dd, clusters = cl),
                         c("y", "zc", "cl"), "cluster"),
    dim_blocked = list(function(dd) difference_in_means(y ~ zb, data = dd, blocks = blk),
                       c("y", "zb", "blk"), "block")
  )
  for (nm in names(cases)) {
    fit <- cases[[nm]][[1]]
    complete <- stats::complete.cases(dm[, cases[[nm]][[2]]])
    expect_warning(
      with_missing <- fit(dm),
      paste0("missingness in the ", cases[[nm]][[3]])
    )
    filtered <- fit(d[complete, ])
    expect_same_inference(with_missing, filtered, nm)
    expect_equal(with_missing$nobs, filtered$nobs, label = paste(nm, "nobs"))
  }
})

# ---- the two solver paths ----

# Built for this section; the rank-deficient copies are made inside the test.
set.seed(42)
n <- 200
dat <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  x3 = rnorm(n),
  z = rnorm(n),
  cl = rep(1:20, each = 10)
)
dat$y <- 1 + dat$x1 + dat$x2 + dat$x3 + rnorm(n)

# The fallback is shared code, but each estimator reaches it through its own
# fit call, and until now only unclustered, unweighted lm_robust() exercised
# it. A wiring mistake at one entry point would be invisible everywhere else.
test_that("the Cholesky path is wired correctly at every entry point", {
  dup <- dat
  dup$x2 <- dup$x1
  dup$w <- runif(n, 0.5, 2)
  dup$Z <- rep(c(0, 1), length.out = n)

  same_both_ways <- function(fitter, label) {
    q <- suppressWarnings(fitter(FALSE))
    ch <- suppressWarnings(fitter(TRUE))
    expect_equal(coef(ch), coef(q), tolerance = 1e-10,
                 info = paste("coefficients,", label))
    expect_equal(ch$std.error, q$std.error, tolerance = 1e-10,
                 info = paste("standard errors,", label))
    expect_equal(sum(is.na(coef(ch))), 1L,
                 info = paste("dropped count,", label))
  }

  same_both_ways(function(tc) lm_robust(y ~ x1 + x2 + x3, data = dup,
                                        clusters = cl, se_type = "CR2",
                                        try_cholesky = tc), "clustered CR2")
  same_both_ways(function(tc) lm_robust(y ~ x1 + x2 + x3, data = dup,
                                        weights = w, try_cholesky = tc),
                 "weighted")
  same_both_ways(function(tc) iv_robust(y ~ x1 + x2 + x3 | x1 + x2 + z,
                                        data = dup, try_cholesky = tc),
                 "iv_robust")

  # lm_lin builds its own centered interactions, so the design it hands the
  # solver is not the one written in the formula.
  lin_q <- lm_lin(y ~ Z, ~ x1 + x3, data = dup, try_cholesky = FALSE)
  lin_ch <- lm_lin(y ~ Z, ~ x1 + x3, data = dup, try_cholesky = TRUE)
  expect_equal(coef(lin_ch), coef(lin_q), tolerance = 1e-10)
  expect_equal(lin_ch$std.error, lin_q$std.error, tolerance = 1e-10)

  # A full-rank fixed-effects fit never reaches the fallback, so it checks
  # that the fast path is right where it actually runs.
  fe_q <- lm_robust(y ~ x1 + x3, data = dat, fixed_effects = ~ cl,
                    try_cholesky = FALSE)
  fe_ch <- lm_robust(y ~ x1 + x3, data = dat, fixed_effects = ~ cl,
                     try_cholesky = TRUE)
  expect_equal(coef(fe_ch), coef(fe_q), tolerance = 1e-10)
  expect_equal(fe_ch$std.error, fe_q$std.error, tolerance = 1e-10)
})

test_that("on a full-rank design try_cholesky gives the default path's answer everywhere", {
  # The test above holds a rank-deficient design, where the Cholesky path falls
  # back to the QR. This one holds every variance path at full rank, which is
  # where the fast path actually runs. Worst measured gap 9.9e-13.
  paths <- list(
    HC0 = function(tc) lm_robust(y ~ x1 + x2 + z, data = d, se_type = "HC0", try_cholesky = tc),
    HC1 = function(tc) lm_robust(y ~ x1 + x2 + z, data = d, se_type = "HC1", try_cholesky = tc),
    HC2 = function(tc) lm_robust(y ~ x1 + x2 + z, data = d, se_type = "HC2", try_cholesky = tc),
    HC3 = function(tc) lm_robust(y ~ x1 + x2 + z, data = d, se_type = "HC3", try_cholesky = tc),
    classical = function(tc) lm_robust(y ~ x1 + x2 + z, data = d, se_type = "classical",
                                       try_cholesky = tc),
    CR0 = function(tc) lm_robust(y ~ x1 + x2 + z, data = d, clusters = cl, se_type = "CR0",
                                 try_cholesky = tc),
    CR2 = function(tc) lm_robust(y ~ x1 + x2 + z, data = d, clusters = cl, se_type = "CR2",
                                 try_cholesky = tc),
    stata = function(tc) lm_robust(y ~ x1 + x2 + z, data = d, clusters = cl, se_type = "stata",
                                   try_cholesky = tc),
    weighted_HC2 = function(tc) lm_robust(y ~ x1 + x2 + z, data = d, weights = w,
                                          try_cholesky = tc),
    weighted_CR2 = function(tc) lm_robust(y ~ x1 + x2 + z, data = d, weights = w, clusters = cl,
                                          try_cholesky = tc),
    fe_two_way_HC3 = function(tc) lm_robust(y ~ x1 + x2 + z, data = d, fixed_effects = ~ g + h,
                                            se_type = "HC3", try_cholesky = tc),
    fe_CR2 = function(tc) lm_robust(y ~ x1 + x2 + z, data = d, fixed_effects = ~ g,
                                    clusters = cl, se_type = "CR2", try_cholesky = tc),
    multivariate_CR2 = function(tc) lm_robust(cbind(y, y2) ~ x1 + x2 + z, data = d,
                                              clusters = cl, try_cholesky = tc),
    iv_CR2 = function(tc) iv_robust(y ~ en + x2 | inst + x2, data = d, clusters = cl,
                                    try_cholesky = tc),
    iv_weighted = function(tc) iv_robust(y ~ en + x2 | inst + x2, data = d, weights = w,
                                         try_cholesky = tc),
    lin_weighted_clustered = function(tc) lm_lin(y ~ z, covariates = ~ x1 + x2, data = d,
                                                 weights = w, clusters = cl, try_cholesky = tc)
  )
  for (nm in names(paths)) {
    qr_fit <- paths[[nm]](FALSE)
    cholesky_fit <- paths[[nm]](TRUE)
    expect_same_inference(cholesky_fit, qr_fit, nm)
    expect_equal(unname(cholesky_fit$fitted.values), unname(qr_fit$fitted.values),
                 tolerance = EQ_TOL, label = paste(nm, "fitted values"))
  }

  diagnosed <- lapply(c(FALSE, TRUE), function(tc) {
    iv_robust(y ~ en + x2 | inst + x2, data = d, diagnostics = TRUE, try_cholesky = tc)
  })
  for (field in c("diagnostic_first_stage_fstatistic", "diagnostic_endogeneity_test")) {
    expect_equal(diagnosed[[2]][[field]], diagnosed[[1]][[field]], tolerance = EQ_TOL,
                 label = field)
  }
})

# ---- absorbed fixed effects are the dummy regression ----

dat <- ref_data_fe()
n <- nrow(dat)
test_that("FE demeaning gives same coefs as dummy regression", {
  m_fe    <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  m_dummy <- lm_robust(y ~ z + x + factor(bl), data = dat)
  expect_equal(unname(coef(m_fe)["z"]), unname(coef(m_dummy)["z"]), tolerance = 1e-9)
  expect_equal(unname(coef(m_fe)["x"]), unname(coef(m_dummy)["x"]), tolerance = 1e-9)
})

test_that("FE residuals equal dummy regression residuals", {
  m_fe    <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  m_dummy <- lm_robust(y ~ z + x + factor(bl), data = dat)
  resid_fe    <- dat$y - m_fe$fitted.values
  resid_dummy <- dat$y - m_dummy$fitted.values
  expect_equal(resid_fe, resid_dummy, tolerance = 1e-9)
})

test_that("HC2 and HC3 now work with one-way FE and match the dummy regression", {
  # These two used to assert an error. The restriction was wrong for any
  # number of FE factors; see test_fe_leverage.R for the identity.
  for (se in c("HC2", "HC3")) {
    fe  <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl, se_type = se)
    dum <- lm_robust(y ~ z + factor(bl), data = dat, se_type = se)
    expect_equal(unname(fe$std.error), unname(dum$std.error["z"]), tolerance = 1e-9)
  }
})

test_that("HC2 and HC3 with two-way FE match the dummy regression", {
  dat2 <- dat
  dat2$bl2 <- factor(rep(1:4, length.out = nrow(dat2)))
  for (se in c("HC2", "HC3")) {
    fe  <- lm_robust(y ~ z, data = dat2, fixed_effects = ~ bl + bl2, se_type = se)
    dum <- lm_robust(y ~ z + factor(bl) + bl2, data = dat2, se_type = se)
    expect_equal(unname(fe$std.error), unname(dum$std.error["z"]), tolerance = 1e-9)
  }
})

test_that("CR2 with FE matches the dummy regression", {
  fe  <- lm_robust(y ~ z, data = dat, fixed_effects = ~ bl, clusters = cl,
                   se_type = "CR2")
  dum <- lm_robust(y ~ z + factor(bl), data = dat, clusters = cl, se_type = "CR2")
  expect_equal(unname(fe$std.error), unname(dum$std.error["z"]), tolerance = 1e-9)
  expect_equal(unname(fe$df), unname(dum$df["z"]), tolerance = 1e-9)
})
test_that("two-way FE converges and gives sensible results", {
  # Two-way FE: block + cluster
  m_2way <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl + cl, se_type = "HC1")
  m_dum  <- lm_robust(y ~ z + factor(bl) + factor(cl), data = dat, se_type = "HC1")
  expect_equal(unname(coef(m_2way)["z"]), unname(coef(m_dum)["z"]), tolerance = 1e-7)
})

test_that("iv_robust with FE coefs match FWL manually", {
  # Demean y, z (endogenous), iv manually then run 2SLS
  y_dm  <- dat$y  - ave(dat$y,  dat$bl, FUN = mean)
  z_dm  <- dat$z  - ave(dat$z,  dat$bl, FUN = mean)
  iv_dm <- dat$iv - ave(dat$iv, dat$bl, FUN = mean)
  dat_dm <- data.frame(y=y_dm, z=z_dm, iv=iv_dm)

  m_fe  <- iv_robust(y ~ z | iv, data = dat, fixed_effects = ~bl, se_type = "HC1")
  m_man <- iv_robust(y ~ z | iv, data = dat_dm, se_type = "HC1")
  expect_equal(unname(coef(m_fe)["z"]), unname(coef(m_man)["z"]), tolerance = 1e-9)
})


test_that("B4: the fixed-effects R-squared is weighted and per outcome", {
  # The FE branch used mean(yoriginal) and raw residuals, so a weighted fit
  # reported an unweighted R-squared, and a multivariate outcome was pooled
  # into one number where the same model with dummies gives one per column.
  set.seed(2); n <- 300
  d <- data.frame(x = rnorm(n), bl = sample(8, n, TRUE), w = runif(n, 0.2, 3))
  d$y <- d$x + d$bl * 0.3 + rnorm(n)
  d$y2 <- d$y + rnorm(n)

  wf <- lm_robust(y ~ x, fixed_effects = ~ bl, data = d, weights = w)
  lw <- summary(lm(y ~ x + factor(bl), data = d, weights = w))
  expect_equal(wf$r.squared, lw$r.squared, tolerance = 1e-12)
  expect_equal(wf$adj.r.squared, lw$adj.r.squared, tolerance = 1e-12)

  uf <- lm_robust(y ~ x, fixed_effects = ~ bl, data = d)
  lu <- summary(lm(y ~ x + factor(bl), data = d))
  expect_equal(uf$r.squared, lu$r.squared, tolerance = 1e-12)

  mv <- lm_robust(cbind(y, y2) ~ x, fixed_effects = ~ bl, data = d)
  dm <- lm_robust(cbind(y, y2) ~ x + factor(bl), data = d)
  expect_equal(length(mv$r.squared), 2L)
  expect_equal(unname(mv$r.squared), unname(dm$r.squared), tolerance = 1e-12)

  ivw <- iv_robust(y ~ x | x, fixed_effects = ~ bl, data = d, weights = w)
  expect_equal(ivw$r.squared, lw$r.squared, tolerance = 1e-12)
})

test_that("absorbed fixed effects with a multivariate outcome are the dummy regression", {
  set.seed(43)
  N <- 40
  d <- data.frame(Y = rnorm(N), Y2 = rnorm(N), Z = rbinom(N, 1, 0.5), X = rnorm(N),
                  B = factor(rep(1:4, each = 10)))
  dummies <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B), data = d)
  absorbed <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B, data = d)
  expect_equal(unname(absorbed$coefficients), unname(dummies$coefficients[c("Z", "X"), ]),
               tolerance = 1e-10)
  expect_equal(unname(absorbed$fitted.values), unname(dummies$fitted.values), tolerance = 1e-8)
})

# ---- Horvitz-Thompson's ways of stating one design ----
#
# condition_prs takes a named probability vector, a per-unit vector, a
# two-column matrix, or a randomizr declaration, and each is turned into
# per-unit probabilities by its own branch. A simple design stated each way is
# one design and must give one estimate, whatever the conditions are called.

test_that("one simple design stated four ways is one Horvitz-Thompson estimate", {
  set.seed(9)
  n_ht <- 60
  p <- 0.4
  ht <- data.frame(y = rnorm(n_ht), z = rbinom(n_ht, 1, p))
  named <- horvitz_thompson(y ~ z, data = ht, condition_prs = c("0" = 1 - p, "1" = p))

  per_unit <- horvitz_thompson(y ~ z, data = ht, condition_prs = rep(p, n_ht))
  two_column <- horvitz_thompson(y ~ z, data = ht,
                                 condition_prs = cbind("0" = rep(1 - p, n_ht), "1" = rep(p, n_ht)))
  ht_labelled <- ht
  ht_labelled$z <- ifelse(ht$z == 1, "treat", "control")
  labelled <- horvitz_thompson(y ~ z, data = ht_labelled,
                               condition_prs = c(control = 1 - p, treat = p),
                               condition1 = "control", condition2 = "treat")
  for (form in list(per_unit = per_unit, two_column = two_column, labelled = labelled)) {
    expect_equal(unname(form$coefficients), unname(named$coefficients), tolerance = EQ_TOL)
    expect_equal(unname(form$std.error), unname(named$std.error), tolerance = EQ_TOL)
  }

  skip_if_not_installed("randomizr")
  declared <- horvitz_thompson(y ~ z, data = ht,
                               condition_prs = randomizr::declare_ra(N = n_ht, prob = p, simple = TRUE))
  expect_equal(unname(declared$coefficients), unname(named$coefficients), tolerance = EQ_TOL)
  expect_equal(unname(declared$std.error), unname(named$std.error), tolerance = EQ_TOL)
})

