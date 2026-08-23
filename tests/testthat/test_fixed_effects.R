library(estimatr)

# Fixed effects are absorbed via Frisch-Waugh-Lovell (FWL) demeaning.
# This gives bit-identical coefficients and residuals to the dummy-variable
# formulation.
#
# SE type support with fixed_effects. Every se_type is available; what differs
# is whether the fixed effects have to be expanded into dummies to get it.
#
#   Absorbed, no dummy matrix ever built:
#     HC0, HC1, HC2, HC3, classical, stata   (no clusters)
#     CR0, stata                             (with clusters)
#
#   Expands the dummies:
#     CR2                                    (with clusters)
#
# HC2 and HC3 need the hat values of the full [X | FE dummies] design, but
# those decompose exactly, for ANY number of factors: h_ii is the demeaned-X
# hat value plus diag(P_D), and fe_leverage() computes diag(P_D) without
# building D. An earlier version of this comment said HC2, HC3 and CR2 all
# errored under fixed_effects; only the first two were ever expensive, and
# neither is any more.
#
# CR2 is the one left that expands. It needs the whole per-cluster BLOCK of the
# hat matrix rather than its diagonal, and no such decomposition exists, so it
# costs roughly O(J^3) in the number of FE levels, exactly as 1.0.6 did. It is
# built only when CR2 is actually asked for, never to satisfy a default: with
# fixed_effects and clusters the default is CR0, with a warning saying so.
# Weighted CR2 with fixed_effects is refused, as it was in 1.0.6.

dat <- ref_data_fe()
n <- nrow(dat)

# ---- FWL equivalence ----

test_that("FE demeaning gives same coefs as dummy regression", {
  m_fe    <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  m_dummy <- lm_robust(y ~ z + x + factor(bl), data = dat)
  expect_equal(unname(coef(m_fe)["z"]), unname(coef(m_dummy)["z"]), tolerance = 1e-10)
  expect_equal(unname(coef(m_fe)["x"]), unname(coef(m_dummy)["x"]), tolerance = 1e-10)
})

test_that("FE residuals equal dummy regression residuals", {
  m_fe    <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  m_dummy <- lm_robust(y ~ z + x + factor(bl), data = dat)
  resid_fe    <- dat$y - m_fe$fitted.values
  resid_dummy <- dat$y - m_dummy$fitted.values
  expect_equal(resid_fe, resid_dummy, tolerance = 1e-10)
})

test_that("FE df.residual is n - k - (B-1) - 1", {
  m <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  # 2 covariates, 20 blocks (19 dummies + intercept absorbed) → 200 - 2 - 20 = 178
  expect_equal(m$df.residual, n - 2L - 20L)
})

# ---- defaults and se_type behaviour ----

test_that("one-way FE default se_type is HC2 (no clusters)", {
  # Was HC1: HC2 used to be unaffordable with absorbed FE. The leverage
  # identity makes it exact and cheap, so the default is the package's
  # ordinary one again, and it matches estimatr for this call.
  # See test_fe_leverage.R.
  m <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl)
  expect_equal(m$se_type, "HC2")
})

test_that("FE default se_type is CR0 (with clusters)", {
  expect_warning(m <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl, clusters = cl),
                 "`se_type` defaults to")
  expect_equal(m$se_type, "CR0")
})

test_that("HC2 and HC3 now work with one-way FE and match the dummy regression", {
  # These two used to assert an error. The restriction was wrong for any
  # number of FE factors; see test_fe_leverage.R for the identity.
  for (se in c("HC2", "HC3")) {
    fe  <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl, se_type = se)
    dum <- lm_robust(y ~ z + factor(bl), data = dat, se_type = se)
    expect_equal(unname(fe$std.error), unname(dum$std.error["z"]))
  }
})

test_that("HC2 and HC3 with two-way FE match the dummy regression", {
  dat2 <- dat
  dat2$bl2 <- factor(rep(1:4, length.out = nrow(dat2)))
  for (se in c("HC2", "HC3")) {
    fe  <- lm_robust(y ~ z, data = dat2, fixed_effects = ~ bl + bl2, se_type = se)
    dum <- lm_robust(y ~ z + factor(bl) + bl2, data = dat2, se_type = se)
    expect_equal(unname(fe$std.error), unname(dum$std.error["z"]))
  }
})

test_that("CR2 with FE matches the dummy regression", {
  fe  <- lm_robust(y ~ z, data = dat, fixed_effects = ~ bl, clusters = cl,
                   se_type = "CR2")
  dum <- lm_robust(y ~ z + factor(bl), data = dat, clusters = cl, se_type = "CR2")
  expect_equal(unname(fe$std.error), unname(dum$std.error["z"]))
  expect_equal(unname(fe$df), unname(dum$df["z"]))
})

test_that("HC1 with FE does not warn", {
  expect_no_warning(lm_robust(y ~ z, data = dat, fixed_effects = ~bl, se_type = "HC1"))
})

test_that("classical with FE does not warn", {
  expect_no_warning(lm_robust(y ~ z, data = dat, fixed_effects = ~bl, se_type = "classical"))
})

# ---- numerical identity with estimatr ----
#
# HC0/HC1/stata/classical/CR0 match estimatr exactly via FWL: these SE types
# use only residuals (and cluster memberships for CR0), not hat values, so the
# FWL demeaned-X model gives the same result as the full [X|FE dummies] model.
#
# HC2, HC3 and CR2 need hat values from the full [X | FE dummies] design. One
# fixed-effect factor supplies HC2 and HC3 from a leverage vector instead (see
# test_fe_leverage.R); everything else expands the dummies, which costs roughly
# O(g^3) in the number of levels. That is why no *default* reaches for them.

test_that("FE HC1 coefs, SEs, df, and fitted values identical to estimatr", {
  m0 <- ref("fe_HC1")
  mz <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl, se_type = "HC1")
  expect_equal(coef(mz),         coef(m0),         tolerance = REF_TOL)
  expect_equal(mz$std.error,     m0$std.error,     tolerance = REF_TOL)
  expect_equal(mz$df.residual,   m0$df.residual)
  expect_equal(mz$fitted.values, m0$fitted.values, tolerance = REF_TOL)
})

test_that("FE HC0 coefs and SEs identical to estimatr", {
  m0 <- ref("fe_HC0")
  mz <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl, se_type = "HC0")
  expect_equal(coef(mz),     coef(m0),     tolerance = REF_TOL)
  expect_equal(mz$std.error, m0$std.error, tolerance = REF_TOL)
})

test_that("FE stata cluster SEs identical to estimatr", {
  m0 <- ref("fe_cl_stata")
  mz <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl,
                  clusters = cl, se_type = "stata")
  expect_equal(coef(mz),     coef(m0),     tolerance = REF_TOL)
  expect_equal(mz$std.error, m0$std.error, tolerance = REF_TOL)
})

test_that("FE default coefs and fitted values identical to estimatr HC1", {
  m0 <- ref("fe_HC1")
  mz <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  expect_equal(coef(mz),         coef(m0),         tolerance = REF_TOL)
  expect_equal(mz$fitted.values, m0$fitted.values, tolerance = REF_TOL)
  expect_equal(mz$df.residual,   m0$df.residual)
})

test_that("two-way FE HC2 and HC3 identical to estimatr", {
  for (se in c("HC2", "HC3")) {
    m0 <- ref(paste0("fe_2way_", se))
    mz <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~ bl + cl, se_type = se)
    expect_equal(coef(mz),     coef(m0),     tolerance = REF_TOL)
    expect_equal(mz$std.error, m0$std.error, tolerance = REF_TOL)
    expect_equal(mz$df,        m0$df,        tolerance = REF_TOL)
  }
})

test_that("FE clustered CR2 identical to estimatr", {
  m0 <- ref("fe_cl_CR2")
  mz <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~ bl, clusters = cl,
                  se_type = "CR2")
  expect_equal(coef(mz),     coef(m0),     tolerance = REF_TOL)
  expect_equal(mz$std.error, m0$std.error, tolerance = REF_TOL)
  expect_equal(mz$df,        m0$df,        tolerance = REF_TOL)
})

test_that("iv_robust FE clustered CR2 identical to estimatr", {
  m0 <- ref("fe_iv_cl_CR2")
  mz <- iv_robust(y ~ z | iv, data = dat, fixed_effects = ~ bl, clusters = cl,
                  se_type = "CR2")
  expect_equal(coef(mz),     coef(m0),     tolerance = REF_TOL)
  expect_equal(mz$std.error, m0$std.error, tolerance = REF_TOL)
  expect_equal(mz$df,        m0$df,        tolerance = REF_TOL)
})

test_that("FE weighted HC1 coefs and SEs identical to estimatr", {
  dat_w <- ref_data_fe_weighted()
  m0 <- ref("fe_w_HC1")
  mz <- lm_robust(y ~ z + x, data = dat_w, fixed_effects = ~bl,
                  weights = w, se_type = "HC1")
  expect_equal(coef(mz),     coef(m0),     tolerance = REF_TOL)
  expect_equal(mz$std.error, m0$std.error, tolerance = REF_TOL)
})

# ---- return object ----

test_that("FE model has proj_r.squared and r.squared", {
  m <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  expect_true(!is.null(m$proj_r.squared))
  expect_true(!is.null(m$r.squared))
  expect_gte(m$r.squared, 0)
  expect_lte(m$r.squared, 1)
})

test_that("FE fes flag is TRUE", {
  m <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl)
  expect_true(m$fes)
})

test_that("no FE fes flag is FALSE", {
  m <- lm_robust(y ~ z, data = dat)
  expect_false(m$fes)
})

test_that("tidy works on FE model", {
  m  <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  td <- tidy(m)
  expect_s3_class(td, "data.frame")
  expect_true("z" %in% td$term)
})

test_that("summary works on FE model", {
  m <- lm_robust(y ~ z + x, data = dat, fixed_effects = ~bl)
  expect_no_error(summary(m))
})

# ---- multi-way FE ----

test_that("two-way FE converges and gives sensible results", {
  # Two-way FE: block + cluster
  m_2way <- lm_robust(y ~ z, data = dat, fixed_effects = ~bl + cl, se_type = "HC1")
  m_dum  <- lm_robust(y ~ z + factor(bl) + factor(cl), data = dat, se_type = "HC1")
  expect_equal(unname(coef(m_2way)["z"]), unname(coef(m_dum)["z"]), tolerance = 1e-7)
})

# ---- iv_robust with FE ----

test_that("iv_robust with FE returns correct class", {
  m <- iv_robust(y ~ z | iv, data = dat, fixed_effects = ~bl, se_type = "HC1")
  expect_s3_class(m, "iv_robust")
  expect_true(m$fes)
})

test_that("iv_robust with FE coefs match FWL manually", {
  # Demean y, z (endogenous), iv manually then run 2SLS
  y_dm  <- dat$y  - ave(dat$y,  dat$bl, FUN = mean)
  z_dm  <- dat$z  - ave(dat$z,  dat$bl, FUN = mean)
  iv_dm <- dat$iv - ave(dat$iv, dat$bl, FUN = mean)
  dat_dm <- data.frame(y=y_dm, z=z_dm, iv=iv_dm)

  m_fe  <- iv_robust(y ~ z | iv, data = dat, fixed_effects = ~bl, se_type = "HC1")
  m_man <- iv_robust(y ~ z | iv, data = dat_dm, se_type = "HC1")
  expect_equal(unname(coef(m_fe)["z"]), unname(coef(m_man)["z"]), tolerance = 1e-10)
})

test_that("iv_robust FE diagnostics are suppressed with warning", {
  expect_warning(
    iv_robust(y ~ z | iv, data = dat, fixed_effects = ~bl, diagnostics = TRUE,
              se_type = "HC1"),
    "Diagnostics"
  )
})

# ---- The fixed-effects term keeps its name when rows are dropped ----

# With one fixed-effects factor, dropping any row for missingness degrades the
# model frame's one-column code matrix to a vector, and rebuilding it as a
# matrix loses the column name. Three things downstream read that name, and all
# three broke silently: `fe_factors()` handed `model.matrix(~ 0 + .)` an
# unnamed data frame and CR2 died with "'.' in formula", the absorbed effects
# came back as "1", "2", "3" instead of "g1", "g2", "g3", and `predict()`
# rejected its own newdata as containing new levels. The suite saw none of it,
# because every fixed-effects test until now used complete data.
test_that("a single FE factor keeps its term name with missing data", {
  set.seed(20260820)
  n <- 400
  d <- data.frame(
    y  = rnorm(n),
    x1 = rnorm(n),
    g  = factor(sample(6, n, TRUE)),
    cl = sample(5, n, TRUE)
  )
  d_na <- d
  d_na$y[c(1, 200, n)] <- NA

  complete_names <- names(lm_robust(y ~ x1, data = d, fixed_effects = ~ g)$fixed_effects)
  for (se in c("HC1", "HC2", "classical", "stata")) {
    m <- lm_robust(y ~ x1, data = d_na, fixed_effects = ~ g, se_type = se)
    expect_equal(names(m$fixed_effects), complete_names)
    expect_equal(names(m$felevels), "g")
  }

  # CR2 is the case that errored outright: it is the one se_type that still
  # builds the dummy matrix, so it is the only one that reads the name through
  # model.matrix().
  m_cr2 <- lm_robust(y ~ x1, data = d_na, fixed_effects = ~ g,
                     clusters = cl, se_type = "CR2")
  expect_equal(names(m_cr2$fixed_effects), complete_names)

  # Absorbed and explicit-dummy CR2 must still agree, which is what says the
  # restored path computes the right thing rather than merely returning.
  m_dum <- lm_robust(y ~ x1 + g, data = d_na, clusters = cl, se_type = "CR2")
  expect_equal(unname(m_cr2$std.error["x1"]), unname(m_dum$std.error["x1"]),
               tolerance = 1e-12)
  expect_equal(unname(m_cr2$df["x1"]), unname(m_dum$df["x1"]), tolerance = 1e-12)

  # predict() matches newdata against these names, so it fails on the fit's own
  # data when they are lost.
  m <- lm_robust(y ~ x1, data = d_na, fixed_effects = ~ g)
  rows <- c("2", "3", "4", "5")
  expect_equal(
    unname(predict(m, newdata = d_na[rows, ])),
    unname(m$fitted.values[rows]),
    tolerance = 1e-12
  )
})

# A factor can carry a declared level the data never uses. When those levels
# sit BELOW the maximum used code they leave a gap, and clean_model_data()
# renumbers around it. When they sit at the top there is no gap to find, so the
# codes run 1..max and `felevels` keeps every declared level, including the
# unused ones. `absorbed_group_effects()` then has more levels than groups, and
# the level with no rows has no representative row to evaluate at.
#
# It must come back NA, with the vector still one entry per level. The index
# that finds representative rows is built by scatter assignment over the codes,
# and an unused level keeps the 0 that vector was initialised with. Zero is a
# DROP index in R, not a missing one, so leaving it there silently shortens the
# result and `dim<-` fails with "dims [product 10] do not match the length of
# object [6]". This is the only configuration that reaches it.
test_that("a fixed-effects factor with unused trailing levels gives NA effects", {
  set.seed(7)
  n <- 600
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  d$f <- factor(sample(1:6, n, TRUE), levels = 1:10)

  m <- lm_robust(y ~ x, data = d, fixed_effects = ~ f, se_type = "HC1")

  expect_equal(length(m$fixed_effects), 10L)
  expect_equal(sum(is.na(m$fixed_effects)), 4L)
  expect_equal(names(m$fixed_effects), paste0("f", c(1, 10, 2:9)))
  expect_equal(length(m$felevels[[1]]), 10L)

  # The estimates themselves must be untouched by any of this.
  m_dum <- lm_robust(y ~ x + f, data = d, se_type = "HC1")
  expect_equal(unname(m$std.error["x"]), unname(m_dum$std.error["x"]),
               tolerance = 1e-12)
  expect_equal(unname(coef(m)["x"]), unname(coef(m_dum)["x"]), tolerance = 1e-12)
  expect_equal(m$df.residual, m_dum$df.residual)

  # A gap BELOW the maximum is renumbered instead, so no level is left empty.
  d$h <- factor(sample(c(1:4, 9, 10), n, TRUE), levels = 1:10)
  mh <- lm_robust(y ~ x, data = d, fixed_effects = ~ h, se_type = "HC1")
  expect_equal(length(mh$fixed_effects), 6L)
  expect_false(anyNA(mh$fixed_effects))
})


# ---- se_type validation runs before the fixed-effects branches ----

test_that("`fixed_effects` does not put an se_type back on the menu", {
  # check_se_type() used to test the fixed-effects cases first and return
  # before the cluster/no-cluster menu, so `fixed_effects` widened the menu:
  # HC2 with `clusters` was accepted and the clusters were then silently
  # ignored, and CR2 with no `clusters` reached the C++ with a NULL cluster
  # vector. 1.0.6 errored on both. The message must not depend on whether
  # `fixed_effects` was passed.
  set.seed(1); n <- 200
  d <- data.frame(y = rnorm(n), x = rnorm(n),
                  bl = sample(10, n, TRUE), cl = sample(20, n, TRUE))

  for (se in c("HC0", "HC1", "HC2", "HC3", "classical")) {
    expect_error(
      lm_robust(y ~ x, fixed_effects = ~ bl, clusters = cl, data = d, se_type = se),
      "when `clusters` are specified", info = se
    )
  }
  for (se in c("CR0", "CR2")) {
    expect_error(
      lm_robust(y ~ x, fixed_effects = ~ bl, data = d, se_type = se),
      "reserved for a case with clusters", info = se
    )
  }
  expect_error(
    iv_robust(y ~ x | bl, fixed_effects = ~ cl, data = d, se_type = "CR2"),
    "reserved for a case with clusters"
  )

  # every combination that was always legal still is, and the defaults hold
  expect_silent(lm_robust(y ~ x, fixed_effects = ~ bl, data = d, se_type = "HC2"))
  expect_equal(lm_robust(y ~ x, fixed_effects = ~ bl, data = d)$se_type, "HC2")
  cr <- suppressWarnings(
    lm_robust(y ~ x, fixed_effects = ~ bl, clusters = cl, data = d, se_type = "CR2")
  )
  expect_equal(cr$se_type, "CR2")
  expect_equal(
    suppressWarnings(lm_robust(y ~ x, fixed_effects = ~ bl, clusters = cl, data = d))$se_type,
    "CR0"
  )
})

# ---- CR2 with a fixed-effects design short of full rank by two or more ----

test_that("absorbed CR2 matches the dummy expansion when two FE columns drop", {
  # getMeatXtX() compacted the design by removing tossed columns in QR pivot
  # order, which is only correct in descending order: the first left shift
  # renumbers every column after it, so the second removal took the wrong one.
  # One redundant column is exact either way, which is why the nested and
  # disconnected two-factor probes all agreed and this stayed hidden. B3 below
  # is a coarsening of A, so the FE design is short by two.
  skip_if_not_installed("clubSandwich")
  set.seed(3); n <- 400
  d <- data.frame(A = sample(1:12, n, TRUE), C = sample(1:6, n, TRUE))
  d$B3 <- ((d$A - 1) %/% 3) + 1
  d$cl <- sample(1:25, n, TRUE)
  d$x <- rnorm(n)
  d$y <- 0.5 * d$x + d$A * 0.1 + d$C * 0.2 + rnorm(n)

  absorbed <- lm_robust(y ~ x, fixed_effects = ~ A + B3 + C, data = d,
                        clusters = cl, se_type = "CR2")
  dummies <- suppressWarnings(
    lm_robust(y ~ x + factor(A) + factor(B3) + factor(C), data = d,
              clusters = cl, se_type = "CR2")
  )
  # three columns are dropped, so the bug had something to reorder
  expect_gte(sum(is.na(dummies$coefficients)), 2L)

  cs <- clubSandwich::vcovCR(
    lm(y ~ x + factor(A) + factor(B3) + factor(C), data = d),
    cluster = d$cl, type = "CR2"
  )
  # tight, not testthat's default 1.5e-8: these are three routes to one number
  # computed in one session, and the gap the bug left was 3.1%
  expect_equal(absorbed$std.error[["x"]], dummies$std.error[["x"]], tolerance = 1e-9)
  expect_equal(absorbed$std.error[["x"]], sqrt(cs["x", "x"]), tolerance = 1e-9)
  expect_equal(absorbed$df[["x"]], dummies$df[["x"]], tolerance = 1e-9)
})


test_that("B3: multi-way demeaning that runs out of sweeps says so", {
  # Alternating projections are exact in one sweep for a single factor and
  # converge geometrically for several, at a rate set by how well connected the
  # factors are. The cap was silent, so a weakly connected design returned
  # coefficients that were wrong rather than approximate, with no signal.
  set.seed(11); n <- 3000
  worker <- sample(300, n, TRUE)
  d <- data.frame(
    worker = worker,
    firm = pmin(((worker - 1) %/% 10) + 1 + sample(0:1, n, TRUE), 30),
    x = rnorm(n)
  )
  d$y <- 1.45 * d$x + rnorm(n) + d$worker * 0.01 + d$firm * 0.02

  expect_warning(
    lm_robust(y ~ x, fixed_effects = ~ worker + firm, data = d),
    "did not converge"
  )
  a <- suppressWarnings(lm_robust(y ~ x, fixed_effects = ~ worker + firm, data = d))
  b <- lm(y ~ x + factor(worker) + factor(firm), data = d)
  # it really has not converged: the answer differs from the exact fit
  expect_gt(abs(coef(a)[["x"]] - coef(b)[["x"]]), 1e-8)

  # a single factor converges in one sweep and must stay silent
  expect_silent(lm_robust(y ~ x, fixed_effects = ~ worker, data = d))
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
