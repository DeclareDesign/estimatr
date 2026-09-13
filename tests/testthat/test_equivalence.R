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
# Both sides are computed in one session. The worst gap measured across every
# pair below is 1.5e-13, and EQ_TOL keeps ~700x over it.
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
