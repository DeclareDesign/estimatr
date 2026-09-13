library(estimatr)

# Change the input in a way whose effect on the answer is known exactly, and
# check the effect.
#
# Every other layer of this suite needs a reference. The 1.0.6 recording cannot
# see an error 1.0.6 already had, and the outside implementations are only as
# hostile as the data they are run on, which elsewhere is well conditioned
# throughout. Neither saw the two defects fixed on 2026-09-11 and 2026-09-12:
# rank detection that depended on a column's units, and a Cholesky path that
# did not detect rank at all. A transformation with a known effect needs no
# reference and does not care where the code came from.
#
# Column scaling, the property the first of those defects broke, is the last
# section. The rank decisions themselves, and the sweep toward singularity, are
# in test_degenerate.R; the two solver paths' agreement is in
# test_equivalence.R.
#
# Its first run found two more: multi-way demeaning that stopped early on an
# outcome in small units, and IV diagnostics that decided endogeneity by name.
#
# Both sides of every comparison are computed in one session, so the floor is
# operation order rather than BLAS. Every assertion below holds at 1e-12, and
# INV_TOL keeps 100x over that for the platforms where the linear algebra
# differs. Multi-way fixed-effect order has its own constant, set where it is
# used.
INV_TOL <- 1e-10

inv_data <- function() {
  set.seed(1)
  n <- 240
  d <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    z = rbinom(n, 1, 0.5),
    w = runif(n, 0.5, 2),
    cl = sample(25, n, replace = TRUE),
    g = factor(sample(letters[1:8], n, replace = TRUE)),
    h = factor(sample(5, n, replace = TRUE)),
    inst = rnorm(n),
    inst2 = rnorm(n)
  )
  d$en <- d$inst + 0.5 * d$inst2 + rnorm(n)
  d$y <- 1 + d$x1 - d$x2 + 0.5 * d$z + d$en + rnorm(n) * (1 + abs(d$x1))
  # Assigned by cluster for the clustered difference in means, and within pairs
  # for the matched-pair one.
  d$zc <- as.integer(d$cl %% 2 == 0)
  d$pair <- rep(seq_len(n / 2), each = 2)
  d$zp <- rep(c(0L, 1L), n / 2)
  d
}

d <- inv_data()
n <- nrow(d)

# One fit per estimator and variance path. Each takes the data and nothing
# else, so a test can hand it a transformed copy.
fitters <- list(
  HC0 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, se_type = "HC0"),
  HC2 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, se_type = "HC2"),
  HC3 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, se_type = "HC3"),
  classical = function(d) lm_robust(y ~ x1 + x2 + z, data = d, se_type = "classical"),
  CR2 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, clusters = cl),
  stata_clustered = function(d) lm_robust(y ~ x1 + x2 + z, data = d, clusters = cl, se_type = "stata"),
  weighted_HC2 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, weights = w),
  weighted_CR2 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, weights = w, clusters = cl),
  fe_HC2 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, fixed_effects = ~ g),
  fe_two_way_HC2 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, fixed_effects = ~ g + h),
  fe_CR2 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, fixed_effects = ~ g,
                                 clusters = cl, se_type = "CR2"),
  fe_weighted_CR0 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, fixed_effects = ~ g,
                                          clusters = cl, weights = w, se_type = "CR0"),
  iv_HC2 = function(d) iv_robust(y ~ en + x1 | inst + inst2 + x1, data = d),
  iv_CR2 = function(d) iv_robust(y ~ en + x1 | inst + inst2 + x1, data = d, clusters = cl),
  iv_fe = function(d) iv_robust(y ~ en + x1 | inst + inst2 + x1, data = d, fixed_effects = ~ g),
  lin = function(d) lm_lin(y ~ z, covariates = ~ x1 + x2, data = d),
  lin_weighted_clustered = function(d) lm_lin(y ~ z, covariates = ~ x1 + x2, data = d,
                                              weights = w, clusters = cl),
  dim = function(d) difference_in_means(y ~ z, data = d),
  dim_blocked = function(d) difference_in_means(y ~ z, data = d, blocks = g),
  dim_clustered = function(d) difference_in_means(y ~ zc, data = d, clusters = cl),
  dim_pairs = function(d) difference_in_means(y ~ zp, data = d, blocks = pair),
  dim_weighted = function(d) difference_in_means(y ~ z, data = d, weights = w)
)

expect_same_inference <- function(fit, base, label, tolerance = INV_TOL) {
  expect_equal(fit$coefficients, base$coefficients, tolerance = tolerance,
               label = paste(label, "coefficients"))
  expect_equal(fit$std.error, base$std.error, tolerance = tolerance,
               label = paste(label, "standard errors"))
  expect_equal(fit$df, base$df, tolerance = tolerance, label = paste(label, "df"))
}

# ---- the outcome's units ----
#
# y -> a + b y moves the intercept to a + b * intercept, every other
# coefficient to b times itself and every standard error to |b| times itself,
# and leaves the degrees of freedom and p-values alone.
#
# The shifts are kept within a few multiples of the outcome's spread. A shift
# far beyond it cancels digits in the residuals, which is a floor set by double
# precision rather than by this package: at y -> 1e8 + y the coefficients agree
# only to 4e-7, and at y -> 1e6 + 1e-6 y to 3e-3.

outcome_maps <- list(
  c(a = 0, b = -1),
  c(a = 0, b = 1e-6),
  c(a = -5, b = 1e6),
  c(a = 7, b = 0.5)
)

test_that("an affine change to the outcome carries through every estimator exactly", {
  for (nm in names(fitters)) {
    base <- fitters[[nm]](d)
    for (ab in outcome_maps) {
      a <- ab[["a"]]
      b <- ab[["b"]]
      moved <- d
      moved$y <- a + b * d$y
      fit <- fitters[[nm]](moved)

      back <- fit$coefficients / b
      intercept <- names(back) == "(Intercept)"
      back[intercept] <- (fit$coefficients[intercept] - a) / b

      label <- sprintf("%s, y -> %g + %g y,", nm, a, b)
      expect_equal(back, base$coefficients, tolerance = INV_TOL,
                   label = paste(label, "coefficients"))
      expect_equal(fit$std.error / abs(b), base$std.error, tolerance = INV_TOL,
                   label = paste(label, "standard errors"))
      expect_equal(fit$df, base$df, tolerance = INV_TOL, label = paste(label, "df"))
      # The intercept is the one coefficient a shift moves by more than a
      # factor, so its test statistic is the one that changes.
      expect_equal(fit$p.value[!intercept], base$p.value[!intercept], tolerance = INV_TOL,
                   label = paste(label, "p-values"))
      if (!is.null(base$fitted.values)) {
        expect_equal((fit$fitted.values - a) / b, base$fitted.values, tolerance = INV_TOL,
                     label = paste(label, "fitted values"))
      }
      if (!is.null(base$residuals)) {
        expect_equal(fit$residuals / b, base$residuals, tolerance = INV_TOL,
                     label = paste(label, "residuals"))
      }
    }
  }
})

test_that("IV diagnostics and linear hypotheses follow the outcome's units", {
  diagnostic_fields <- c("diagnostic_first_stage_fstatistic",
                         "diagnostic_endogeneity_test",
                         "diagnostic_overid_test")
  for (st in c("classical", "HC2", "CR2")) {
    cluster_ids <- if (st == "CR2") d$cl else NULL
    base <- iv_robust(y ~ en + x1 | inst + inst2 + x1, data = d, se_type = st,
                      clusters = cluster_ids, diagnostics = TRUE)
    for (ab in outcome_maps) {
      moved <- d
      moved$y <- ab[["a"]] + ab[["b"]] * d$y
      fit <- iv_robust(y ~ en + x1 | inst + inst2 + x1, data = moved, se_type = st,
                       clusters = cluster_ids, diagnostics = TRUE)
      for (field in diagnostic_fields) {
        expect_equal(fit[[field]], base[[field]], tolerance = INV_TOL,
                     label = sprintf("%s, %s, y -> %g + %g y", field, st, ab[["a"]], ab[["b"]]))
      }
    }
  }

  base <- lh_robust(y ~ x1 + x2 + z, data = d, linear_hypothesis = "z + 2*x1 = 0")
  for (ab in outcome_maps) {
    moved <- d
    moved$y <- ab[["a"]] + ab[["b"]] * d$y
    fit <- lh_robust(y ~ x1 + x2 + z, data = moved, linear_hypothesis = "z + 2*x1 = 0")
    expect_equal(fit$lh$coefficients / ab[["b"]], base$lh$coefficients, tolerance = INV_TOL)
    expect_equal(fit$lh$std.error / abs(ab[["b"]]), base$lh$std.error, tolerance = INV_TOL)
  }
})

# Horvitz-Thompson divides each arm's weighted total by N rather than by the
# arm's own weight, so it follows the outcome's scale and not its location: a
# shift of a moves the estimate by a times the difference of the two arms'
# summed inverse probabilities over N. That is the estimator, and the second
# assertion is here so that a Hajek normalization cannot arrive unannounced.

test_that("Horvitz-Thompson follows the outcome's scale, and its location as the estimator says", {
  p <- 0.4
  prs <- c("0" = 1 - p, "1" = p)
  base <- horvitz_thompson(y ~ z, data = d, condition_prs = prs)
  for (b in c(-1, 1e-6, 1e6)) {
    moved <- d
    moved$y <- b * d$y
    fit <- horvitz_thompson(y ~ z, data = moved, condition_prs = prs)
    expect_equal(fit$coefficients / b, base$coefficients, tolerance = INV_TOL)
    expect_equal(fit$std.error / abs(b), base$std.error, tolerance = INV_TOL)
  }

  a <- 7
  moved <- d
  moved$y <- a + d$y
  fit <- horvitz_thompson(y ~ z, data = moved, condition_prs = prs)
  expect_equal(
    fit$coefficients,
    base$coefficients + a * (sum(d$z / p) - sum((1 - d$z) / (1 - p))) / n,
    tolerance = INV_TOL
  )

  set.seed(4)
  shuffled <- horvitz_thompson(y ~ z, data = d[sample(n), ], condition_prs = prs)
  expect_equal(shuffled$coefficients, base$coefficients, tolerance = INV_TOL)
  expect_equal(shuffled$std.error, base$std.error, tolerance = INV_TOL)
})

# ---- rows in any order ----
#
# Clustering and fixed effects both sort the data internally, and the fitted
# values and residuals have to be put back in the order the rows arrived in.

test_that("the order of the rows changes no estimate, and fits come back in row order", {
  set.seed(2)
  perm <- sample(n)
  shuffled <- d[perm, ]
  for (nm in names(fitters)) {
    base <- fitters[[nm]](d)
    fit <- fitters[[nm]](shuffled)
    expect_same_inference(fit, base, nm)
    if (!is.null(base$fitted.values)) {
      expect_equal(fit$fitted.values, base$fitted.values[perm], tolerance = INV_TOL,
                   label = paste(nm, "fitted values"))
    }
    if (!is.null(base$residuals)) {
      expect_equal(fit$residuals, base$residuals[perm], tolerance = INV_TOL,
                   label = paste(nm, "residuals"))
    }
  }
})

# ---- terms in any order ----

test_that("the order of regressors and instruments changes no estimate", {
  pairs <- list(
    HC2 = list(
      function(d) lm_robust(y ~ x1 + x2 + z, data = d),
      function(d) lm_robust(y ~ z + x2 + x1, data = d)
    ),
    CR2 = list(
      function(d) lm_robust(y ~ x1 + x2 + z, data = d, clusters = cl),
      function(d) lm_robust(y ~ z + x2 + x1, data = d, clusters = cl)
    ),
    fe_CR2 = list(
      function(d) lm_robust(y ~ x1 + x2 + z, data = d, fixed_effects = ~ g,
                            clusters = cl, se_type = "CR2"),
      function(d) lm_robust(y ~ z + x2 + x1, data = d, fixed_effects = ~ g,
                            clusters = cl, se_type = "CR2")
    ),
    iv_CR2 = list(
      function(d) iv_robust(y ~ en + x1 | inst + inst2 + x1, data = d, clusters = cl),
      function(d) iv_robust(y ~ x1 + en | x1 + inst2 + inst, data = d, clusters = cl)
    )
  )
  for (nm in names(pairs)) {
    a <- pairs[[nm]][[1]](d)
    b <- pairs[[nm]][[2]](d)
    terms <- names(a$coefficients)
    expect_equal(b$coefficients[terms], a$coefficients, tolerance = INV_TOL,
                 label = paste(nm, "coefficients"))
    expect_equal(b$std.error[terms], a$std.error, tolerance = INV_TOL,
                 label = paste(nm, "standard errors"))
    expect_equal(b$df[terms], a$df, tolerance = INV_TOL, label = paste(nm, "df"))
  }
})

# Reordering multi-way fixed effects leaves the coefficients within 5e-13, but
# the standard errors move by up to 9.3e-10, measured across twenty seeds on
# two- and three-way designs with HC1, HC2 and HC3, weighted and not. That is
# three orders of magnitude above what any one-way or unabsorbed fit shows
# here, so it gets its own constant rather than loosening INV_TOL: 1e-8 keeps
# ~11x over the worst case.
INV_TOL_MULTIWAY <- 1e-8

test_that("the order of multi-way fixed effects changes no estimate", {
  pairs <- list(
    two_way_HC2 = list(~ g + h, ~ h + g, "HC2", FALSE),
    two_way_weighted_HC3 = list(~ g + h, ~ h + g, "HC3", TRUE),
    two_way_HC1 = list(~ g + h, ~ h + g, "HC1", FALSE)
  )
  for (nm in names(pairs)) {
    p <- pairs[[nm]]
    weights <- if (p[[4]]) d$w else NULL
    a <- lm_robust(y ~ x1 + x2 + z, data = d, fixed_effects = p[[1]],
                   se_type = p[[3]], weights = weights)
    b <- lm_robust(y ~ z + x2 + x1, data = d, fixed_effects = p[[2]],
                   se_type = p[[3]], weights = weights)
    terms <- names(a$coefficients)
    expect_equal(b$coefficients[terms], a$coefficients, tolerance = INV_TOL,
                 label = paste(nm, "coefficients"))
    expect_equal(b$std.error[terms], a$std.error, tolerance = INV_TOL_MULTIWAY,
                 label = paste(nm, "standard errors"))
    expect_equal(b$r.squared, a$r.squared, tolerance = INV_TOL,
                 label = paste(nm, "R-squared"))
  }
})

# ---- the weights' units ----
#
# Multiplying every weight by a constant changes nothing: not the estimates,
# not the variance, and not the R-squared. test_lm_lin_equivalence.R checks
# this for lm_lin alone.

test_that("multiplying every weight by a constant changes nothing", {
  weighted <- list(
    HC0 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, weights = w, se_type = "HC0"),
    HC1 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, weights = w, se_type = "HC1"),
    HC2 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, weights = w, se_type = "HC2"),
    HC3 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, weights = w, se_type = "HC3"),
    classical = function(d) lm_robust(y ~ x1 + x2 + z, data = d, weights = w, se_type = "classical"),
    CR0 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, weights = w, clusters = cl, se_type = "CR0"),
    CR2 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, weights = w, clusters = cl),
    stata = function(d) lm_robust(y ~ x1 + x2 + z, data = d, weights = w, clusters = cl, se_type = "stata"),
    fe_HC2 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, weights = w, fixed_effects = ~ g),
    iv_HC2 = function(d) iv_robust(y ~ en + x1 | inst + inst2 + x1, data = d, weights = w),
    iv_CR2 = function(d) iv_robust(y ~ en + x1 | inst + inst2 + x1, data = d, weights = w, clusters = cl),
    dim = function(d) difference_in_means(y ~ z, data = d, weights = w),
    dim_blocked = function(d) difference_in_means(y ~ z, data = d, weights = w, blocks = g),
    dim_clustered = function(d) difference_in_means(y ~ zc, data = d, weights = w, clusters = cl)
  )
  for (nm in names(weighted)) {
    base <- weighted[[nm]](d)
    for (k in c(1e-6, 7, 1e6)) {
      moved <- d
      moved$w <- k * d$w
      fit <- weighted[[nm]](moved)
      label <- sprintf("%s, weights times %g,", nm, k)
      expect_same_inference(fit, base, label)
      if (!is.null(base$r.squared)) {
        expect_equal(fit$r.squared, base$r.squared, tolerance = INV_TOL,
                     label = paste(label, "R-squared"))
      }
    }
  }
})

# ---- how groups are labelled ----
#
# A cluster, block or fixed-effect identifier is a label. Integer codes,
# character strings, a factor whose levels run in the opposite order, large
# non-integer doubles and a random relabelling must all give the same fit.

relabellings <- list(
  character = function(v) paste0("id_", v),
  large_double = function(v) as.numeric(v) * 1e9 + 0.5,
  reversed_factor = function(v) factor(v, levels = rev(sort(unique(v)))),
  shuffled_codes = function(v) {
    set.seed(3)
    u <- sort(unique(v))
    unname(setNames(sample(length(u)) + 1000L, u)[as.character(v)])
  }
)

test_that("relabelling the clusters changes nothing", {
  clustered <- list(
    CR0 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, clusters = cl, se_type = "CR0"),
    CR2 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, clusters = cl),
    stata = function(d) lm_robust(y ~ x1 + x2 + z, data = d, clusters = cl, se_type = "stata"),
    iv_CR2 = function(d) iv_robust(y ~ en + x1 | inst + inst2 + x1, data = d, clusters = cl),
    lin_CR2 = function(d) lm_lin(y ~ z, covariates = ~ x1, data = d, clusters = cl),
    dim = function(d) difference_in_means(y ~ zc, data = d, clusters = cl)
  )
  for (nm in names(clustered)) {
    base <- clustered[[nm]](d)
    for (rl in names(relabellings)) {
      moved <- d
      moved$cl <- relabellings[[rl]](d$cl)
      expect_same_inference(clustered[[nm]](moved), base, paste(nm, rl))
    }
  }
})

test_that("relabelling fixed effects and blocks changes nothing", {
  grouped <- list(
    fe_HC2 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, fixed_effects = ~ g),
    fe_CR2 = function(d) lm_robust(y ~ x1 + x2 + z, data = d, fixed_effects = ~ g,
                                   clusters = cl, se_type = "CR2"),
    fe_two_way_HC3 = function(d) lm_robust(y ~ x1 + x2 + z, data = d,
                                           fixed_effects = ~ g + h, se_type = "HC3"),
    dim_blocked = function(d) difference_in_means(y ~ z, data = d, blocks = g)
  )
  for (nm in names(grouped)) {
    base <- grouped[[nm]](d)
    for (rl in names(relabellings)) {
      moved <- d
      moved$g <- relabellings[[rl]](as.integer(d$g))
      expect_same_inference(grouped[[nm]](moved), base, paste(nm, rl))
    }
  }
})

# ---- reparameterising the design ----
#
# Replacing the regressors X by X A for an invertible A spans the same column
# space, so fitted values, residuals and hat values are unchanged, the
# coefficients become A^-1 b, and every variance estimate becomes
# A^-1 V A^-T. The A below shifts, rescales and mixes columns at once, which
# covers the translation and cross-column cases that column scaling alone does
# not.
#
# The CR2 degrees of freedom are not compared: Satterthwaite's approximation is
# computed per linear combination, so a mixed coefficient has its own.

reparam <- function(d) {
  d$u1 <- 3 + 2 * d$x1 - d$x2
  d$u2 <- -1 + d$x1 + 0.5 * d$x2
  d$u3 <- 0.25 * d$x1 + d$z
  d
}
# Columns of A are u1, u2 and u3 written in (Intercept), x1, x2, z.
A <- cbind(c(1, 0, 0, 0), c(3, 2, -1, 0), c(-1, 1, 0.5, 0), c(0, 0.25, 0, 1))
A_inv <- solve(A)

test_that("mixing the regressors transforms coefficients and variance by A", {
  dr <- reparam(d)
  for (st in c("HC0", "HC1", "HC2", "HC3", "classical", "CR0", "CR2", "stata")) {
    cluster_ids <- if (st %in% c("CR0", "CR2", "stata")) dr$cl else NULL
    for (weighted in c(FALSE, TRUE)) {
      weights <- if (weighted) dr$w else NULL
      label <- paste0(st, if (weighted) ", weighted")
      fx <- lm_robust(y ~ x1 + x2 + z, data = dr, se_type = st,
                      clusters = cluster_ids, weights = weights)
      fu <- lm_robust(y ~ u1 + u2 + u3, data = dr, se_type = st,
                      clusters = cluster_ids, weights = weights)
      expect_equal(unname(fu$coefficients), drop(A_inv %*% fx$coefficients),
                   tolerance = INV_TOL, label = paste(label, "coefficients"))
      expect_equal(unname(fu$vcov), unname(A_inv %*% fx$vcov %*% t(A_inv)),
                   tolerance = INV_TOL, label = paste(label, "vcov"))
      expect_equal(fu$fitted.values, fx$fitted.values, tolerance = INV_TOL,
                   label = paste(label, "fitted values"))
    }
  }
})

test_that("under fixed effects the shifts are absorbed and only the mixing remains", {
  dr <- reparam(d)
  M_inv <- solve(A[2:4, 2:4])
  cases <- list(
    two_way_HC1 = list(fe = ~ g + h, se = "HC1", cl = NULL),
    two_way_HC2 = list(fe = ~ g + h, se = "HC2", cl = NULL),
    two_way_HC3 = list(fe = ~ g + h, se = "HC3", cl = NULL),
    one_way_CR2 = list(fe = ~ g, se = "CR2", cl = dr$cl)
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    fx <- lm_robust(y ~ x1 + x2 + z, data = dr, fixed_effects = cs$fe,
                    se_type = cs$se, clusters = cs$cl)
    fu <- lm_robust(y ~ u1 + u2 + u3, data = dr, fixed_effects = cs$fe,
                    se_type = cs$se, clusters = cs$cl)
    expect_equal(unname(fu$coefficients), drop(M_inv %*% fx$coefficients),
                 tolerance = INV_TOL, label = paste(nm, "coefficients"))
    expect_equal(unname(fu$vcov), unname(M_inv %*% fx$vcov %*% t(M_inv)),
                 tolerance = INV_TOL, label = paste(nm, "vcov"))
  }
})

# For 2SLS, X -> X A with the instruments unchanged gives the same
# transformation, and mixing the instruments within their own span changes
# nothing at all, the diagnostics included. The Wu-Hausman test is the one that
# leaned on rank detection before 2026-09-11, which is why it is here.
#
# The last spelling spans the exogenous x1 without naming it. Endogeneity used
# to be decided by name, which counted x1 as endogenous there and returned a
# Wu-Hausman statistic of 0.144 on 2 degrees of freedom against 0.109 on 1.

test_that("mixing the regressors of a 2SLS fit transforms it by A", {
  di <- d
  di$ue <- 3 + d$en + 2 * d$x1
  di$ux <- -1 + 0.5 * d$x1
  A_iv_inv <- solve(cbind(c(1, 0, 0), c(3, 1, 2), c(-1, 0, 0.5)))
  for (st in c("HC0", "HC2", "HC3", "CR2")) {
    cluster_ids <- if (st == "CR2") di$cl else NULL
    fx <- iv_robust(y ~ en + x1 | inst + inst2 + x1, data = di, se_type = st,
                    clusters = cluster_ids)
    fu <- iv_robust(y ~ ue + ux | inst + inst2 + x1, data = di, se_type = st,
                    clusters = cluster_ids)
    expect_equal(unname(fu$coefficients), drop(A_iv_inv %*% fx$coefficients),
                 tolerance = INV_TOL, label = paste(st, "coefficients"))
    expect_equal(unname(fu$vcov), unname(A_iv_inv %*% fx$vcov %*% t(A_iv_inv)),
                 tolerance = INV_TOL, label = paste(st, "vcov"))
  }
})

test_that("mixing the instruments within their span changes nothing, diagnostics included", {
  di <- d
  di$v1 <- 1 + d$inst + 2 * d$inst2 + 3 * d$x1
  di$v2 <- d$inst2 - d$x1
  di$v3 <- d$x1 + d$inst
  for (st in c("classical", "HC2", "CR2")) {
    cluster_ids <- if (st == "CR2") di$cl else NULL
    fz <- iv_robust(y ~ en + x1 | inst + inst2 + x1, data = di, se_type = st,
                    clusters = cluster_ids, diagnostics = TRUE)
    respelled <- list(
      named = iv_robust(y ~ en + x1 | v1 + v2 + x1, data = di, se_type = st,
                        clusters = cluster_ids, diagnostics = TRUE),
      unnamed = iv_robust(y ~ en + x1 | v1 + v2 + v3, data = di, se_type = st,
                          clusters = cluster_ids, diagnostics = TRUE)
    )
    for (nm in names(respelled)) {
      fv <- respelled[[nm]]
      expect_same_inference(fv, fz, paste(st, nm))
      for (field in c("diagnostic_first_stage_fstatistic",
                      "diagnostic_endogeneity_test",
                      "diagnostic_overid_test")) {
        expect_equal(fv[[field]], fz[[field]], tolerance = INV_TOL,
                     label = paste(st, nm, field))
      }
    }
  }
})

# ---- shifting a regressor ----
#
# Column normalization makes rank detection indifferent to a column's scale,
# and it can do nothing about its location: x + c for large c is genuinely
# close to the intercept. So the invariance holds up to where the design stops
# being full rank in lm()'s judgement, and past that point the requirement is
# that estimatr and lm() agree about it.

test_that("shifting a regressor moves only the intercept", {
  base <- lm_robust(y ~ x1 + x2 + z, data = d)
  slopes <- c("x1", "x2", "z")
  for (shift in c(10, 100)) {
    moved <- d
    moved$x1 <- d$x1 + shift
    fit <- lm_robust(y ~ x1 + x2 + z, data = moved)
    expect_equal(fit$coefficients[slopes], base$coefficients[slopes], tolerance = INV_TOL)
    expect_equal(fit$std.error[slopes], base$std.error[slopes], tolerance = INV_TOL)
    expect_equal(
      fit$coefficients[["(Intercept)"]],
      base$coefficients[["(Intercept)"]] - shift * base$coefficients[["x1"]],
      tolerance = INV_TOL
    )
  }
})

test_that("far enough out, both paths drop what lm() drops", {
  for (shift in c(1e5, 1e7, 1e9)) {
    moved <- d
    moved$x1 <- d$x1 + shift
    dropped_by_lm <- is.na(coef(lm(y ~ x1 + x2 + z, data = moved)))
    for (tc in c(FALSE, TRUE)) {
      fit <- suppressWarnings(lm_robust(y ~ x1 + x2 + z, data = moved, try_cholesky = tc))
      expect_equal(is.na(fit$coefficients), dropped_by_lm,
                   label = sprintf("dropped columns at shift %g, try_cholesky = %s", shift, tc))
    }
  }
})

# ---- a column's units ----


# Rank detection must not depend on the units a column is measured in.
# stats::lm() gets this from LINPACK dqrdc2, which compares each column's
# remaining norm against its own original norm. Eigen's setThreshold()
# compares every pivot against the largest pivot in the matrix, so before
# lm_solver() and getMeatXtX() normalized their columns, one regressor in
# large units pushed the others under the threshold and they were dropped as
# collinear on a full-rank design: a wrong answer with a warning, not an
# error. The property below is stronger than any fixed design.

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

powers <- 0:12

# Refit with x2 multiplied by 10^k, then undo the scaling on the x2 row. Every
# coefficient and standard error must come back to what the unscaled fit gave.
expect_scale_invariant <- function(fitter) {
  base <- fitter(dat)
  for (k in powers) {
    scaled <- dat
    scaled$x2 <- scaled$x2 * 10^k
    fit <- fitter(scaled)

    coefs <- coef(fit)
    ses <- fit$std.error
    coefs["x2"] <- coefs["x2"] * 10^k
    ses["x2"] <- ses["x2"] * 10^k

    expect_false(anyNA(coefs), label = paste("no coefficient dropped at k =", k))
    expect_equal(coefs, coef(base), tolerance = 1e-8,
                 info = paste("coefficients at k =", k))
    expect_equal(ses, base$std.error, tolerance = 1e-8,
                 info = paste("standard errors at k =", k))
  }
}

test_that("lm_robust is invariant to column scaling", {
  expect_scale_invariant(function(d) lm_robust(y ~ x1 + x2 + x3, data = d))
})

test_that("clustered lm_robust is invariant to column scaling", {
  expect_scale_invariant(
    function(d) lm_robust(y ~ x1 + x2 + x3, clusters = cl, data = d)
  )
})

# getMeatXtX() is the path HC2, HC3, and CR2 read the hat values off, and it
# carries the same threshold as lm_solver(). A fix applied to only one of the
# two leaves the variance read off a rank the coefficients were not fitted at.
test_that("every hat-value se_type is invariant to column scaling", {
  for (se_type in c("HC0", "HC1", "HC2", "HC3", "classical")) {
    local({
      this_type <- se_type
      expect_scale_invariant(
        function(d) lm_robust(y ~ x1 + x2 + x3, se_type = this_type, data = d)
      )
    })
  }
  expect_scale_invariant(
    function(d) lm_robust(y ~ x1 + x2 + x3, clusters = cl,
                          se_type = "CR0", data = d)
  )
})

test_that("iv_robust is invariant to column scaling", {
  expect_scale_invariant(
    function(d) iv_robust(y ~ x1 + x2 + x3 | x1 + z + x3, data = d)
  )
  expect_scale_invariant(
    function(d) iv_robust(y ~ x1 + x2 + x3 | x1 + z + x3,
                          clusters = cl, data = d)
  )
})

# The Cholesky path is a second rank determination, and Eigen's LLT reports
# success on a numerically singular Gram matrix, so info() alone never caught
# a rank-deficient design. Normalizing the columns makes each L_ii the
# column's own residual norm, which is dqrdc2's test, and the path falls back
# to the QR below the same 1e-7.
test_that("try_cholesky reaches the same answer as the QR path", {
  expect_scale_invariant(
    function(d) lm_robust(y ~ x1 + x2 + x3, data = d, try_cholesky = TRUE)
  )

  qr_fit <- lm_robust(y ~ x1 + x2 + x3, data = dat, try_cholesky = FALSE)
  ch_fit <- lm_robust(y ~ x1 + x2 + x3, data = dat, try_cholesky = TRUE)
  expect_equal(coef(ch_fit), coef(qr_fit), tolerance = 1e-10)
  expect_equal(ch_fit$std.error, qr_fit$std.error, tolerance = 1e-10)
})

test_that("a full-rank design in large units agrees with lm", {
  d <- dat
  d$x2 <- d$x2 * 1e9
  fit <- lm_robust(y ~ x1 + x2 + x3, data = d, se_type = "classical")
  ref <- lm(y ~ x1 + x2 + x3, data = d)
  expect_false(anyNA(coef(fit)))
  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-10)
  expect_equal(unname(fit$std.error),
               unname(summary(ref)$coefficients[, 2]), tolerance = 1e-10)
})

