library(estimatr)

# Every comparison in this file is between two routes to one number, computed
# in one session on one BLAS: the absorbed fit and the same model with the
# dummies written out. testthat's default 1.5e-8 is far looser than that can be
# held, and a drift below it would pass in silence, so the comparisons here are
# `expect_same()`. Use plain `expect_equal()` for anything that is not an
# absorbed-versus-expanded pair; comparisons against the 1.0.6 recording stay at
# REF_TOL, which is a different question.
#
# 1e-9 is set from the measured worst case across all 14,125 values compared
# here, which is 5.5e-11 on a disconnected two-way design where the
# pseudo-inverse is doing the work. That leaves ~18x of headroom, and this
# number is set on macOS and has to survive Ubuntu, Windows and R-devel, where
# the linear algebra is not the same. An earlier draft used 1e-10, which left
# 1.8x and would have been a coin flip on the CI matrix for no gain: the point
# of this constant is to catch a drift three orders of magnitude below
# testthat's default, and 1e-9 still does that with 15x to spare.
FE_TOL <- 1e-9
expect_same <- function(object, expected, ...) {
  expect_equal(object, expected, ..., tolerance = FE_TOL)
}

# HC2 and HC3 with absorbed fixed effects.
#
# These used to be refused, on the grounds that they need hat values from the
# full [dummies | X] design and absorption only leaves the demeaned ones. The
# hat values in fact split exactly, for ANY number of FE factors:
#
#   P_[X | D] = P_D + P_{M_D X}
#
# so h_ii = h_ii(demeaned X) + diag(P_D)_i, and the second term costs a
# factorisation of order min(levels) rather than a dummy hat matrix. With one
# factor P_D is diagonal and the term is just w_i / sum(w over i's group).
#
# The tests below pin the identity against the two references that matter:
# writing the dummies out by hand, and estimatr itself.

configs <- ref_data_leverage()

# ---- against explicit dummies ----

for (nm in names(configs)) {
  for (se in c("HC2", "HC3")) {
    test_that(paste0("absorbed FE equals explicit dummies: ", nm, ", ", se), {
      d <- configs[[nm]]
      fe  <- lm_robust(y ~ x + z, fixed_effects = ~ g, data = d, se_type = se)
      dum <- lm_robust(y ~ x + z + factor(g), data = d, se_type = se)
      expect_same(unname(fe$std.error),
                   unname(dum$std.error[c("x", "z")]))
      expect_same(unname(fe$coefficients),
                   unname(dum$coefficients[c("x", "z")]))
    })
  }
}

test_that("the identity survives weights", {
  for (nm in names(configs)) {
    d <- configs[[nm]]
    for (se in c("HC2", "HC3")) {
      fe  <- lm_robust(y ~ x + z, fixed_effects = ~ g, data = d,
                       se_type = se, weights = wts)
      dum <- lm_robust(y ~ x + z + factor(g), data = d,
                       se_type = se, weights = wts)
      expect_same(unname(fe$std.error), unname(dum$std.error[c("x", "z")]),
                   info = paste(nm, se))
    }
  }
})

# ---- against estimatr ----

test_that("absorbed FE agrees with estimatr for HC2 and HC3", {
  for (nm in names(configs)) {
    d <- configs[[nm]]
    for (se in c("HC2", "HC3")) {
      a <- lm_robust(y ~ x + z, fixed_effects = ~ g, data = d, se_type = se)
      b <- ref(paste0("lev_", nm, "_", se))
      # `b` is the 1.0.6 recording, not a same-session computation, so these
      # stay at REF_TOL: comparing a fixture more tightly would test the
      # runner's BLAS rather than this package.
      expect_equal(unname(a$std.error), unname(b$std.error),
                   tolerance = REF_TOL, info = paste(nm, se))
      # df, and so every p-value and interval, has to match too. Comparing
      # only the standard error is what let the absorbed rank go missing from
      # the residual degrees of freedom without a test noticing.
      expect_equal(unname(a$df), unname(b$df), tolerance = REF_TOL, info = paste(nm, se))
      expect_equal(unname(a$p.value), unname(b$p.value),
                   tolerance = REF_TOL, info = paste(nm, se))
      expect_equal(unname(a$conf.low), unname(b$conf.low),
                   tolerance = REF_TOL, info = paste(nm, se))
    }
  }
})

test_that("absorbed FE agrees with estimatr under weights", {
  d <- configs$unbalanced
  for (se in c("HC2", "HC3")) {
    a <- lm_robust(y ~ x + z, fixed_effects = ~ g, data = d,
                   se_type = se, weights = wts)
    b <- ref(paste0("levw_unbalanced_", se))
    expect_equal(unname(a$std.error), unname(b$std.error), tolerance = REF_TOL)
    expect_equal(unname(a$df), unname(b$df), tolerance = REF_TOL)
    expect_equal(unname(a$p.value), unname(b$p.value), tolerance = REF_TOL)
  }
})

# ---- defaults ----

test_that("one-way FE keeps the package default of HC2", {
  # It used to fall back to "stata" because HC2 was unaffordable. It is not
  # unaffordable any more, and HC2 is also what estimatr returns for this call,
  # so falling back would be a silent disagreement with the released package.
  d <- configs$balanced
  expect_same(lm_robust(y ~ x + z, fixed_effects = ~ g, data = d)$se_type, "HC2")
})

test_that("multi-way FE keeps the HC2 default too, and stays silent", {
  # This used to fall back to HC1 and warn, because HC2 meant expanding the
  # dummies. The identity covers any number of factors, so there is nothing
  # left to trade away: HC2 is the default here as it is everywhere else, and
  # as it was in 1.0.6.
  set.seed(3)
  d <- data.frame(y = rnorm(300), x = rnorm(300),
                  a = factor(rep(1:15, each = 20)),
                  b = factor(rep(1:20, times = 15)))
  expect_no_warning(fit <- lm_robust(y ~ x, fixed_effects = ~ a + b, data = d))
  expect_same(fit$se_type, "HC2")
})

# ---- what the dummy expansion supplies, and what it costs ----

# The one-way identity above is a shortcut, not the only route. Where it does
# not apply, the fixed effects are expanded into dummies and the hat values
# come off the full design, which is what 1.0.6 did for every one of these
# cases. The answers below must therefore equal the explicit-dummy fit exactly.

test_that("two-way FE gets HC2 and HC3 from the leverage identity", {
  set.seed(3)
  d <- data.frame(y = rnorm(300), x = rnorm(300),
                  a = factor(rep(1:15, each = 20)),
                  b = factor(rep(1:20, times = 15)))
  for (se in c("HC2", "HC3")) {
    fe  <- lm_robust(y ~ x, fixed_effects = ~ a + b, data = d, se_type = se)
    dum <- lm_robust(y ~ x + a + b, data = d, se_type = se)
    expect_same(unname(fe$std.error), unname(dum$std.error["x"]))
    expect_same(unname(fe$df), unname(dum$df["x"]))
  }
})

test_that("CR2 with one FE factor matches the dummy expansion", {
  # CR2's adjustment comes from cluster-level blocks of the hat matrix, not
  # from h_ii, so the one-way identity says nothing about it and the dummies
  # have to be built even here.
  d <- configs$balanced
  d$cl <- factor(rep(1:10, length.out = nrow(d)))
  fe  <- lm_robust(y ~ x + z, fixed_effects = ~ g, clusters = cl, data = d,
                   se_type = "CR2")
  dum <- lm_robust(y ~ x + z + g, clusters = cl, data = d, se_type = "CR2")
  expect_same(unname(fe$std.error), unname(dum$std.error[c("x", "z")]))
  expect_same(unname(fe$df), unname(dum$df[c("x", "z")]))
})

test_that("weighted CR2 with FE is still refused, as in 1.0.6", {
  d <- configs$balanced
  d$cl <- factor(rep(1:10, length.out = nrow(d)))
  d$w  <- runif(nrow(d), 1, 2)
  expect_error(
    lm_robust(y ~ x + z, fixed_effects = ~ g, clusters = cl, weights = w,
              data = d, se_type = "CR2"),
    "cannot be combined with `weights`"
  )
})

test_that("the clustered FE default warns rather than switching silently", {
  d <- configs$balanced
  d$cl <- factor(rep(1:10, length.out = nrow(d)))
  # 1.0.6 defaulted to CR2 here.
  expect_warning(
    fit <- lm_robust(y ~ x + z, fixed_effects = ~ g, clusters = cl, data = d),
    "`se_type` defaults to"
  )
  expect_same(fit$se_type, "CR0")
  expect_no_warning(
    lm_robust(y ~ x + z, fixed_effects = ~ g, clusters = cl, data = d, se_type = "CR0")
  )
})

test_that("the default warnings fire once per session, not once per call", {
  # helper-warning-frequency.R forces every frequency warning to signal so the
  # rest of the suite does not depend on file order. This one test needs the
  # real behaviour, so it puts the option back and clears the spent budget.
  old <- options(rlib_warning_verbosity = "default")
  on.exit(options(old), add = TRUE)
  rlang::reset_warning_verbosity("estimatr_fe_clustered_default")

  d <- configs$balanced
  d$cl <- factor(rep(1:10, length.out = nrow(d)))
  expect_warning(lm_robust(y ~ x + z, fixed_effects = ~ g, clusters = cl, data = d),
                 "`se_type` defaults to")
  expect_no_warning(lm_robust(y ~ x + z, fixed_effects = ~ g, clusters = cl, data = d))
  expect_no_warning(lm_robust(y ~ x + z, fixed_effects = ~ g, clusters = cl, data = d))
})

test_that("one-way unclustered FE keeps the HC2 default and stays silent", {
  d <- configs$balanced
  expect_no_warning(fit <- lm_robust(y ~ x + z, fixed_effects = ~ g, data = d))
  expect_same(fit$se_type, "HC2")
})

# ---- the leverage identity itself ----

test_that("the FE leverage vector is the group weight share", {
  d <- configs$unbalanced
  # clean_model_data() carries the fixed effects as an integer code matrix
  # with the level names alongside; build that shape here.
  g_fac <- as.factor(d$g)
  md <- list(
    fixed_effects = matrix(as.integer(g_fac), ncol = 1L,
                           dimnames = list(NULL, "g")),
    fe_level_names = list(g = levels(g_fac)),
    outcome = as.matrix(d$y),
    design_matrix = cbind("(Intercept)" = 1, x = d$x),
    weights = d$wts,
    terms = terms(y ~ x, data = d)
  )
  out <- estimatr:::demean_fes(md)
  # demean_fes() carries the group codes; the leverage is computed from them
  # when the requested se_type wants it.
  expect_same(estimatr:::fe_leverage(out$fe_codes, d$wts)$leverage,
               d$wts / ave(d$wts, d$g, FUN = sum))
})

test_that("no leverage vector is produced for two-way FE", {
  set.seed(3)
  d <- data.frame(y = rnorm(300), x = rnorm(300),
                  a = factor(rep(1:15, each = 20)),
                  b = factor(rep(1:20, times = 15)))
  md <- list(
    fixed_effects = d[c("a", "b")],
    outcome = as.matrix(d$y),
    design_matrix = cbind("(Intercept)" = 1, x = d$x),
    weights = NULL,
    terms = terms(y ~ x, data = d)
  )
  expect_null(estimatr:::demean_fes(md)$fe_leverage)
})


# ---- the multi-way identity itself ----

# diag(P_D) is what fe_leverage() returns. Checked against the projection
# written out in full, across factor counts, orderings, weights, and a
# DISCONNECTED design, where D is rank deficient and the pseudo-inverse is
# what makes the projection still come out right.

dense_fe_design <- function(codes) {
  do.call(cbind, lapply(seq_along(codes), function(k) {
    f <- factor(codes[[k]], levels = seq_len(max(codes[[k]])))
    m <- model.matrix(~ 0 + f)
    if (k == 1) m else m[, -1, drop = FALSE]
  }))
}

dense_fe_leverage <- function(codes, w = NULL) {
  n <- length(codes[[1]])
  if (is.null(w)) w <- rep(1, n)
  D <- dense_fe_design(codes)
  Dw <- sqrt(w) * D
  qr_D <- qr(Dw)
  Q <- qr.Q(qr_D)[, seq_len(qr_D$rank), drop = FALSE]
  rowSums(Q^2)
}

fe_cfgs <- list(
  "one factor"        = c(20L),
  "two factors"       = c(30L, 7L),
  "wide then narrow"  = c(60L, 3L),
  "narrow then wide"  = c(3L, 60L),
  "equal widths"      = c(15L, 15L),
  "three factors"     = c(20L, 8L, 4L),
  "four factors"      = c(12L, 7L, 5L, 3L)
)

for (nm in names(fe_cfgs)) {
  for (wtd in c(FALSE, TRUE)) {
    test_that(paste0("fe_leverage equals diag(P_D): ", nm,
                     if (wtd) ", weighted" else ""), {
      set.seed(17)
      n <- 900
      codes <- lapply(fe_cfgs[[nm]], function(g) {
        z <- sample(g, n, TRUE); z[seq_len(g)] <- seq_len(g); as.integer(z)
      })
      w <- if (wtd) runif(n, 0.2, 5) else NULL
      expect_same(estimatr:::fe_leverage(codes, w)$leverage,
                   dense_fe_leverage(codes, w))
    })
  }
}

test_that("fe_leverage handles a disconnected design", {
  # No observation links the first block of levels to the second, so D is rank
  # deficient and a plain solve would fail.
  set.seed(23)
  codes <- list(as.integer(c(sample(1:10, 300, TRUE), sample(11:20, 300, TRUE))),
                as.integer(c(sample(1:4, 300, TRUE),  sample(5:8, 300, TRUE))))
  expect_same(estimatr:::fe_leverage(codes)$leverage, dense_fe_leverage(codes))
  w <- runif(600, 0.3, 3)
  expect_same(estimatr:::fe_leverage(codes, w)$leverage,
               dense_fe_leverage(codes, w))
  # D loses one column per extra connected component, and the rank comes back
  # from the same eigendecomposition that produced the leverage.
  expect_same(estimatr:::fe_leverage(codes)$rank, qr(dense_fe_design(codes))$rank)
  # rank-only must agree with the full call, and skip the vector
  cheap <- estimatr:::fe_leverage(codes, leverage = FALSE)
  expect_null(cheap$leverage)
  expect_same(cheap$rank, estimatr:::fe_leverage(codes)$rank)
})

test_that("a singleton fixed-effect group has leverage exactly one", {
  codes <- list(as.integer(c(1L, sample(2:20, 499, TRUE))),
                as.integer(sample(1:5, 500, TRUE)))
  expect_same(estimatr:::fe_leverage(codes)$leverage[1], 1)
})

test_that("multi-way HC2 and HC3 equal the explicit-dummy fit", {
  set.seed(29)
  n <- 800
  d <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n),
                  wt = runif(n, 0.3, 3),
                  a = factor(sample(25, n, TRUE)),
                  b = factor(sample(6, n, TRUE)),
                  cc = factor(sample(4, n, TRUE)))
  for (se in c("HC2", "HC3")) {
    fe  <- lm_robust(y ~ x + z, fixed_effects = ~ a + b, data = d, se_type = se)
    dum <- lm_robust(y ~ x + z + a + b, data = d, se_type = se)
    expect_same(unname(fe$std.error), unname(dum$std.error[c("x", "z")]))

    fe3  <- lm_robust(y ~ x + z, fixed_effects = ~ a + b + cc, data = d, se_type = se)
    dum3 <- lm_robust(y ~ x + z + a + b + cc, data = d, se_type = se)
    expect_same(unname(fe3$std.error), unname(dum3$std.error[c("x", "z")]))

    few  <- lm_robust(y ~ x + z, fixed_effects = ~ a + b, data = d,
                      weights = wt, se_type = se)
    dumw <- lm_robust(y ~ x + z + a + b, data = d, weights = wt, se_type = se)
    expect_same(unname(few$std.error), unname(dumw$std.error[c("x", "z")]))
  }
})

test_that("iv_robust gets multi-way HC2 and HC3 from the identity too", {
  set.seed(31)
  n <- 900
  d <- data.frame(x = rnorm(n), iv = rnorm(n),
                  a = factor(sample(20, n, TRUE)), b = factor(sample(5, n, TRUE)))
  d$z <- d$iv + rnorm(n)
  d$y <- d$z + rnorm(n)
  for (se in c("HC2", "HC3")) {
    fe  <- iv_robust(y ~ z | iv, fixed_effects = ~ a + b, data = d, se_type = se)
    dum <- iv_robust(y ~ z + a + b | iv + a + b, data = d, se_type = se)
    expect_same(unname(fe$std.error), unname(dum$std.error["z"]))
    expect_same(unname(fe$coefficients), unname(dum$coefficients["z"]))
  }
})


# ---- multi-way against estimatr 1.0.6 ----

# 2.0 answers these from the leverage identity where 1.0.6 expanded the
# dummies. The numbers must be the same either way, degrees of freedom
# included.

test_that("multi-way FE agrees with estimatr, including df", {
  fe <- ref_data_fe()
  for (se in c("HC2", "HC3")) {
    a <- lm_robust(y ~ z + x, data = fe, fixed_effects = ~ bl + cl, se_type = se)
    b <- ref(paste0("fe_2way_", se))
    expect_equal(unname(a$std.error), unname(b$std.error), tolerance = REF_TOL)
    expect_equal(unname(a$df), unname(b$df), tolerance = REF_TOL)
    expect_equal(unname(a$p.value), unname(b$p.value), tolerance = REF_TOL)
  }
})

test_that("multi-way weighted and three-way agree with estimatr", {
  fe <- ref_data_fe_multiway()
  for (se in c("HC2", "HC3")) {
    a <- lm_robust(y ~ z + x, data = fe, fixed_effects = ~ bl + cl,
                   weights = w, se_type = se)
    b <- ref(paste0("fe_2way_w_", se))
    expect_equal(unname(a$std.error), unname(b$std.error), tolerance = REF_TOL)
    expect_equal(unname(a$df), unname(b$df), tolerance = REF_TOL)

    a3 <- lm_robust(y ~ z + x, data = fe, fixed_effects = ~ bl + cl + c3,
                    se_type = se)
    b3 <- ref(paste0("fe_3way_", se))
    expect_equal(unname(a3$std.error), unname(b3$std.error), tolerance = REF_TOL)
    expect_same(unname(a3$df), unname(b3$df))
  }
})

test_that("the multi-way default returns what 1.0.6 returned", {
  fe <- ref_data_fe_multiway()
  a <- lm_robust(y ~ z + x, data = fe, fixed_effects = ~ bl + cl)
  b <- ref("fe_2way_default")
  expect_same(a$se_type, "HC2")
  expect_equal(unname(a$std.error), unname(b$std.error), tolerance = REF_TOL)
  expect_same(unname(a$df), unname(b$df))
})


# ---- a fixed-effect design that is rank deficient ----

test_that("a nested fixed-effect factor matches the explicit-dummy fit", {
  # `c3` is a deterministic function of `cl`, so it adds no columns the other
  # factors do not already span and the FE design is rank deficient by 4.
  #
  # 1.0.6 does not get this right. It expands the dummies and lets a pivoted
  # QR drop the redundant columns, but reads the hat values off the padded
  # design, so its absorbed answer disagrees with its OWN explicit-dummy fit
  # (0.152538 against 0.149296 here). Taking diag(P_D) through a pseudo-inverse
  # is exact whatever the rank, so 2.0 returns the dummy-regression answer.
  set.seed(42)
  n <- 200
  d <- data.frame(y = rnorm(n), x = rnorm(n), z = rbinom(n, 1, 0.5),
                  bl = rep(1:20, each = 10), cl = rep(1:10, 20),
                  c3 = rep(1:5, length.out = 200))
  for (se in c("HC2", "HC3")) {
    fe  <- lm_robust(y ~ z + x, data = d, fixed_effects = ~ bl + cl + c3,
                     se_type = se)
    # The explicit fit drops the redundant columns and says so; the absorbed
    # fit has no columns to drop, and handles the deficiency in the projection.
    expect_warning(
      dum <- lm_robust(y ~ z + x + factor(bl) + factor(cl) + factor(c3),
                       data = d, se_type = se),
      "collinear"
    )
    expect_same(unname(fe$std.error), unname(dum$std.error[c("z", "x")]))
    expect_same(unname(fe$df), unname(dum$df[c("z", "x")]))
  }
})

test_that("the exact FE rank is used, not the nominal level count", {
  set.seed(42)
  n <- 200
  d <- data.frame(y = rnorm(n), x = rnorm(n),
                  bl = rep(1:20, each = 10), cl = rep(1:10, 20),
                  c3 = rep(1:5, length.out = 200))
  # nominal would be (20 + 10 + 5) - 3 + 1 = 33; the true rank is 29.
  fit <- lm_robust(y ~ x, data = d, fixed_effects = ~ bl + cl + c3)
  expect_same(fit$df.residual, n - 1L - 29L)
})


# ---- absorbing must never change the answer, however degenerate the design ----

# Absorbing fixed effects is a computational choice, so an absorbed fit and the
# same fit with the dummies written out have to agree for EVERY se_type, not
# just the ones that read hat values. That fails if the rank correction is
# taken from the nominal level count, because a factor spanned by the others
# contributes fewer than `levels - 1` columns.
#
# Canvassed against the ecosystem on these designs: base `lm()` and `plm`
# (twoways within) both report the exact rank; `fixest::feols` reports the
# nominal one and issues no note, so it differs from its own dummy regression
# here. estimatr follows lm().

degenerate_designs <- list(
  "nested factor" = function() {
    set.seed(42)
    n <- 200
    data.frame(y = rnorm(n), x = rnorm(n),
               bl = rep(1:20, each = 10), cl = rep(1:10, 20),
               # a deterministic function of `cl`, so it spans nothing new
               c3 = rep(1:5, length.out = n))
  },
  "disconnected two-way" = function() {
    set.seed(7)
    m <- 300
    data.frame(y = rnorm(2 * m), x = rnorm(2 * m),
               bl = c(sample(1:10, m, TRUE), sample(11:20, m, TRUE)),
               cl = c(sample(1:4, m, TRUE), sample(5:8, m, TRUE)),
               c3 = 1L)
  }
)

for (nm in names(degenerate_designs)) {
  test_that(paste0("absorbed equals explicit dummies on a ", nm), {
    d <- degenerate_designs[[nm]]()
    fe_form <- if (nm == "nested factor") ~ bl + cl + c3 else ~ bl + cl
    dum_form <- if (nm == "nested factor") {
      y ~ x + factor(bl) + factor(cl) + factor(c3)
    } else {
      y ~ x + factor(bl) + factor(cl)
    }
    for (se in c("HC0", "HC1", "HC2", "HC3", "classical")) {
      fe <- lm_robust(y ~ x, fixed_effects = fe_form, data = d, se_type = se)
      dum <- suppressWarnings(lm_robust(dum_form, data = d, se_type = se))
      expect_same(unname(fe$std.error), unname(dum$std.error["x"]),
                   info = paste(nm, se))
      expect_same(unname(fe$df), unname(dum$df["x"]), info = paste(nm, se))
      expect_same(unname(fe$df.residual), unname(dum$df.residual),
                   info = paste(nm, se))
    }
  })
}

test_that("the FE rank matches a pivoted QR of the dummy design", {
  d <- degenerate_designs[["nested factor"]]()
  D <- model.matrix(~ 0 + factor(bl) + factor(cl) + factor(c3), data = d)
  # rank of the FE design alone, less nothing: qr() is the reference
  fit <- lm_robust(y ~ x, fixed_effects = ~ bl + cl + c3, data = d)
  expect_same(fit$df.residual, nrow(d) - 1L - qr(D)$rank)
  # and the nominal count would have been wrong
  nominal <- (20 + 10 + 5) - 3 + 1
  expect_lt(qr(D)$rank, nominal)
})

test_that("a singleton FE group agrees absorbed and expanded, for every se_type", {
  # One group holds a single observation, so its leverage is 1 and the fit is
  # exact for it. The absorbed fit has always returned a finite standard error
  # here, because that observation contributes nothing to the meat. The
  # explicit-dummy fit used to return NaN instead: the wider design puts its
  # computed hat value marginally above 1, and HC2 divided by the resulting
  # negative number. The leverage guard resolves both to the same zero
  # contribution, so the two routes now agree rather than disagreeing by an
  # accident of which design matrix was formed (estimatr #395).
  set.seed(11)
  k <- 200
  d <- data.frame(y = rnorm(k), x = rnorm(k), g = c(1L, sample(2:20, k - 1, TRUE)))
  expect_same(sum(d$g == 1L), 1L)

  fe <- lm_robust(y ~ x, fixed_effects = ~ g, data = d, se_type = "HC2")
  expect_true(is.finite(fe$std.error[["x"]]))
  expect_same(fe$df.residual, lm(y ~ x + factor(g), data = d)$df.residual)

  for (se in c("HC1", "HC2", "HC3")) {
    a <- lm_robust(y ~ x, fixed_effects = ~ g, data = d, se_type = se)
    # the expanded design is where the hat value lands above 1, so it warns
    b <- suppressWarnings(
      lm_robust(y ~ x + factor(g), data = d, se_type = se)
    )
    expect_true(is.finite(b$std.error[["x"]]), info = se)
    expect_same(unname(a$std.error), unname(b$std.error["x"]), info = se)
  }
})


test_that("a single-level fixed effect alongside others is absorbed, not an error", {
  # `const` has one level, so it contributes an intercept and `g`'s dummies
  # already span it. The fit is therefore the one-way fit, for every se_type.
  # 1.0.6 refused this combination outright.
  set.seed(4)
  n <- 200
  d <- data.frame(y = rnorm(n), x = rnorm(n), g = sample(10, n, TRUE), const = 1L)
  for (se in c("HC0", "HC1", "HC2", "HC3", "classical")) {
    both <- lm_robust(y ~ x, fixed_effects = ~ g + const, data = d, se_type = se)
    one <- lm_robust(y ~ x, fixed_effects = ~ g, data = d, se_type = se)
    expect_same(unname(both$std.error), unname(one$std.error), info = se)
    expect_same(both$df.residual, one$df.residual, info = se)
  }
  # and the rank is the one-way rank, not one more
  expect_same(
    estimatr:::fe_leverage(list(as.integer(d$g), rep(1L, n)))$rank,
    estimatr:::fe_leverage(list(as.integer(d$g)))$rank
  )
})
