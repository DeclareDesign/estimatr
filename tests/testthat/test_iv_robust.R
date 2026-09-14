library(estimatr)

# iv_robust's own behaviour: the 2SLS estimate, its residuals, its diagnostics,
# and the regressions filed against it.
#
# HC0 through HC3, CR0 and CR2 against sandwich and clubSandwich are in
# test_vs_sandwich.R and test_vs_clubsandwich.R; Stata's ivregress in
# test_vs_stata.R; reparameterisations and respelled instruments in
# test_invariance.R; underidentified models in test_degenerate.R.

iv_data <- function() {
  set.seed(42)
  n <- 200
  d <- data.frame(
    x = rnorm(n),
    w = rnorm(n),
    inst = rnorm(n),
    inst2 = rnorm(n),
    wt = runif(n, 0.5, 2)
  )
  d$en <- d$inst + 0.5 * d$inst2 + rnorm(n)
  d$y <- 1 + d$en + d$w + rnorm(n)
  d
}

d <- iv_data()

# Measured gaps against AER below are at most 4e-14.
IV_TOL <- 1e-10

test_that("2SLS matches AER::ivreg, with an exogenous covariate and with weights", {
  # Classical, where AER's own vcov() is the same estimator. The robust
  # variances are compared through sandwich in test_vs_sandwich.R.
  skip_if_not_installed("AER")
  fit <- iv_robust(y ~ en + w | inst + w, data = d, se_type = "classical")
  ref <- AER::ivreg(y ~ en + w | inst + w, data = d)
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = IV_TOL)
  expect_equal(unname(fit$vcov), unname(vcov(ref)), tolerance = IV_TOL)

  fit <- iv_robust(y ~ en + w | inst + w, data = d, se_type = "classical", weights = wt)
  ref <- AER::ivreg(y ~ en + w | inst + w, data = d, weights = wt)
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = IV_TOL)
  expect_equal(unname(fit$vcov), unname(vcov(ref)), tolerance = IV_TOL)
})

test_that("just identified with one instrument, the slope is the ratio of covariances", {
  fit <- iv_robust(y ~ en | inst, data = d)
  expect_equal(fit$coefficients[["en"]], cov(d$y, d$inst) / cov(d$en, d$inst),
               tolerance = IV_TOL)
})

test_that("#345: iv_robust returns structural residuals", {
  m <- iv_robust(mpg ~ wt | am, data = mtcars)
  X <- cbind(1, mtcars$wt)
  expect_equal(unname(m$residuals), as.vector(mtcars$mpg - X %*% coef(m)),
               tolerance = 1e-12)
})

# ---- diagnostics ----

test_that("classical diagnostics match AER: weak instruments, Wu-Hausman, Sargan", {
  skip_if_not_installed("AER")
  fit <- iv_robust(y ~ en + w | inst + inst2 + w, data = d, se_type = "classical",
                   diagnostics = TRUE)
  ref <- summary(AER::ivreg(y ~ en + w | inst + inst2 + w, data = d),
                 diagnostics = TRUE)$diagnostics

  weak <- fit$diagnostic_first_stage_fstatistic
  expect_equal(weak[["value"]], ref["Weak instruments", "statistic"], tolerance = IV_TOL)
  expect_equal(weak[["nomdf"]], ref["Weak instruments", "df1"])
  expect_equal(weak[["dendf"]], ref["Weak instruments", "df2"])
  expect_equal(weak[["p.value"]], ref["Weak instruments", "p-value"], tolerance = IV_TOL)

  wu <- fit$diagnostic_endogeneity_test
  expect_equal(wu[["value"]], ref["Wu-Hausman", "statistic"], tolerance = IV_TOL)
  expect_equal(wu[["numdf"]], ref["Wu-Hausman", "df1"])
  expect_equal(wu[["dendf"]], ref["Wu-Hausman", "df2"])
  expect_equal(wu[["p.value"]], ref["Wu-Hausman", "p-value"], tolerance = IV_TOL)

  overid <- fit$diagnostic_overid_test
  expect_equal(overid[["value"]], ref["Sargan", "statistic"], tolerance = IV_TOL)
  expect_equal(overid[["df"]], ref["Sargan", "df1"])
  expect_equal(overid[["p.value"]], ref["Sargan", "p-value"], tolerance = IV_TOL)
})

test_that("robust weak-instrument and Wu-Hausman tests match AER given sandwich's variance", {
  # AER computes both as Wald tests when handed a variance function, which
  # makes it a reference for the robust versions too. It computes no robust
  # over-identification test; Wooldridge's score test is held to Stata's
  # `estat overid` in test_vs_stata.R.
  #
  # This file runs before test_vs_sandwich.R loads the ivreg package, whose
  # summary method then shadows AER's and cannot read an AER fit's hat values.
  skip_if_not_installed("AER")
  skip_if_not_installed("sandwich")
  fml <- y ~ en + w | inst + inst2 + w
  ref_fit <- AER::ivreg(fml, data = d)
  for (ty in c("HC0", "HC1", "HC2", "HC3")) {
    fit <- iv_robust(fml, data = d, se_type = ty, diagnostics = TRUE)
    ref <- summary(ref_fit, vcov. = function(obj) sandwich::vcovHC(obj, type = ty),
                   diagnostics = TRUE)$diagnostics

    weak <- fit$diagnostic_first_stage_fstatistic
    expect_equal(weak[["value"]], ref["Weak instruments", "statistic"], tolerance = IV_TOL,
                 label = paste(ty, "weak-instrument F"))
    expect_equal(weak[["nomdf"]], ref["Weak instruments", "df1"])
    expect_equal(weak[["p.value"]], ref["Weak instruments", "p-value"], tolerance = IV_TOL,
                 label = paste(ty, "weak-instrument p"))

    wu <- fit$diagnostic_endogeneity_test
    expect_equal(wu[["value"]], ref["Wu-Hausman", "statistic"], tolerance = IV_TOL,
                 label = paste(ty, "Wu-Hausman"))
    expect_equal(wu[["numdf"]], ref["Wu-Hausman", "df1"])
    expect_equal(wu[["dendf"]], ref["Wu-Hausman", "df2"])
    expect_equal(wu[["p.value"]], ref["Wu-Hausman", "p-value"], tolerance = IV_TOL,
                 label = paste(ty, "Wu-Hausman p"))
  }
})

test_that("over-identified diagnostics work for every se_type, not just classical", {
  # `first_stage_fits` has its columns renamed `fit_<endog>`, and the rewritten
  # robust branch then indexed it by the bare endogenous names, so every
  # over-identified fit with a non-classical se_type died with "subscript out
  # of bounds". The classical branch took a different function and was fine,
  # which is why the whole path had no test.
  set.seed(2)
  n <- 500
  dd <- data.frame(z1 = rnorm(n), z2 = rnorm(n), z3 = rnorm(n), w = rnorm(n))
  dd$x <- dd$z1 + 0.5 * dd$z2 + 0.3 * dd$z3 + rnorm(n)
  dd$y <- dd$x + dd$w + rnorm(n)
  fml <- y ~ x + w | z1 + z2 + z3 + w

  for (ty in c("classical", "HC0", "HC1", "HC2", "HC3")) {
    m <- iv_robust(fml, data = dd, se_type = ty, diagnostics = TRUE)
    ov <- m$diagnostic_overid_test
    expect_equal(names(ov), c("value", "df", "p.value"), info = ty)
    expect_true(is.finite(ov[["value"]]), info = ty)
    expect_equal(unname(ov[["df"]]), 2, info = ty)   # 3 instruments, 1 endogenous
  }

  # The robust branch is Wooldridge's score test, a different statistic from
  # Sargan's, and every non-classical se_type shares it.
  classical <- iv_robust(fml, data = dd, se_type = "classical", diagnostics = TRUE)
  rb <- vapply(c("HC0", "HC1", "HC2", "HC3"), function(ty) {
    iv_robust(fml, data = dd, se_type = ty, diagnostics = TRUE)$diagnostic_overid_test[["value"]]
  }, numeric(1))
  expect_equal(unname(diff(range(rb))), 0)
  expect_false(isTRUE(all.equal(rb[["HC0"]], classical$diagnostic_overid_test[["value"]])))
})

test_that("a clustered over-identification test sums the score's variance within clusters", {
  # estimatr 1.0.6, and the rewrite until this test, reported the unclustered
  # statistic for a clustered fit. Stata's estat overid computes nothing after
  # vce(cluster), so the reference is the definition: s' S^-1 s, with s the sum
  # of the score contributions and S their outer product summed within clusters.
  set.seed(3)
  n_clusters <- 40
  n <- 10 * n_clusters
  dd <- data.frame(g = rep(seq_len(n_clusters), each = 10), id = seq_len(n),
                   z1 = rnorm(n), z2 = rnorm(n), z3 = rnorm(n), w = rnorm(n))
  dd$x <- dd$z1 + 0.5 * dd$z2 + 0.3 * dd$z3 + rnorm(n)
  dd$y <- dd$x + dd$w + rnorm(n) + rep(rnorm(n_clusters), each = 10)
  fml <- y ~ x + w | z1 + z2 + z3 + w
  overid <- function(...) {
    iv_robust(fml, data = dd, diagnostics = TRUE, ...)$diagnostic_overid_test
  }

  X <- cbind(1, dd$x, dd$w)
  Z <- cbind(1, dd$z1, dd$z2, dd$z3, dd$w)
  xhat <- qr.fitted(qr(Z), X)
  u <- as.vector(dd$y - X %*% qr.coef(qr(xhat), dd$y))
  k <- qr.resid(qr(xhat), cbind(dd$z2, dd$z3)) * u
  s <- colSums(k)
  by_definition <- drop(s %*% solve(crossprod(rowsum(k, dd$g)), s))

  for (ty in c("CR0", "stata", "CR2")) {
    ov <- overid(clusters = g, se_type = ty)
    expect_equal(ov[["value"]], by_definition, tolerance = 1e-10, label = ty)
    expect_equal(unname(ov[["df"]]), 2, label = ty)
  }
  expect_false(isTRUE(all.equal(by_definition, overid(se_type = "HC0")[["value"]])))

  # One observation per cluster is the unclustered statistic.
  expect_equal(overid(clusters = id, se_type = "CR0")[["value"]],
               overid(se_type = "HC0")[["value"]], tolerance = 1e-12)

  # With no more clusters than restrictions S is singular, and the statistic
  # would equal the number of clusters whatever the data.
  dd$g2 <- rep(1:2, each = n / 2)
  expect_warning(ov <- overid(clusters = g2, se_type = "CR0"),
                 "no more clusters than overidentifying restrictions")
  expect_true(is.na(ov[["value"]]))
})

test_that("with weights, each diagnostic is the test on the weighted model under the fit's variance", {
  # estimatr 1.0.6 returned NA for a weighted over-identification test silently.
  # Each statistic is now the unweighted one on the weighted fit's transformed
  # data, sqrt(w) times each column, under the variance the fit uses for its
  # coefficients: Sargan's for classical, the score test for robust. That
  # transformation is the reference.
  set.seed(4)
  n <- 300
  dd <- data.frame(g = rep(1:30, each = 10), z1 = rnorm(n), z2 = rnorm(n),
                   w = rnorm(n), wt = runif(n, 0.5, 2))
  dd$x <- dd$z1 + 0.5 * dd$z2 + rnorm(n)
  dd$y <- dd$x + dd$w + rnorm(n) + rep(rnorm(30), each = 10)
  fml <- y ~ x + w | z1 + z2 + w
  overid <- function(...) {
    iv_robust(fml, data = dd, diagnostics = TRUE, ...)$diagnostic_overid_test[["value"]]
  }

  sw <- sqrt(dd$wt)
  X <- cbind(1, dd$x, dd$w) * sw
  Z <- cbind(1, dd$z1, dd$z2, dd$w) * sw
  xhat <- qr.fitted(qr(Z), X)
  u <- as.vector(dd$y * sw - X %*% qr.coef(qr(xhat), dd$y * sw))
  k <- qr.resid(qr(xhat), cbind(dd$z2 * sw)) * u
  s <- colSums(k)
  by_definition <- drop(s %*% solve(crossprod(k), s))
  by_definition_cl <- drop(s %*% solve(crossprod(rowsum(k, dd$g)), s))

  expect_silent(hc <- iv_robust(fml, data = dd, weights = wt, se_type = "HC1",
                                diagnostics = TRUE))
  expect_equal(hc$diagnostic_overid_test[["value"]], by_definition, tolerance = 1e-10)
  expect_true(is.finite(hc$diagnostic_endogeneity_test[["value"]]))
  expect_equal(overid(weights = wt, clusters = g, se_type = "CR2"), by_definition_cl,
               tolerance = 1e-10)

  # Rescaling the weights changes nothing, and unit weights give the unweighted test.
  dd$wt10 <- 10 * dd$wt
  dd$one <- 1
  expect_equal(overid(weights = wt10, se_type = "HC1"), by_definition, tolerance = 1e-10)
  expect_equal(overid(weights = one, se_type = "HC1"), overid(se_type = "HC1"),
               tolerance = 1e-10)

  # Classical: Sargan's statistic, n times the uncentered R^2 of the transformed
  # residuals on the transformed instruments.
  sargan_by_definition <- n * sum(qr.fitted(qr(Z), u)^2) / sum(u^2)
  expect_silent(cl <- iv_robust(fml, data = dd, weights = wt, se_type = "classical",
                                diagnostics = TRUE))
  expect_equal(cl$diagnostic_overid_test[["value"]], sargan_by_definition, tolerance = 1e-10)
  expect_true(is.finite(cl$diagnostic_endogeneity_test[["value"]]))

  # A zero-weight row is not an observation: both statistics equal the fit with
  # those rows dropped.
  dd$wt0 <- replace(dd$wt, 1:5, 0)
  for (ty in c("classical", "HC1")) {
    expect_equal(overid(weights = wt0, se_type = ty),
                 iv_robust(fml, data = dd[-(1:5), ], weights = wt, se_type = ty,
                           diagnostics = TRUE)$diagnostic_overid_test[["value"]],
                 tolerance = 1e-10, label = ty)
  }
})

test_that("printed diagnostics carry every first stage and its degrees of freedom", {
  # The table read the first-stage df as "numdf", a name the vector does not
  # have, so Df1 printed NA on every fit; with two endogenous regressors the
  # statistic and p-value printed NA too, since those entries are named
  # "<var>:value".
  one <- summary(iv_robust(mpg ~ hp + gear | wt + am + gear, data = mtcars,
                           diagnostics = TRUE))
  tab <- build_ivreg_diagnostics_mat(one)
  fs <- one$diagnostic_first_stage_fstatistic
  expect_equal(rownames(tab), c("Weak instruments", "Wu-Hausman", "Score (robust)"))
  expect_equal(unname(tab["Weak instruments", ]),
               unname(fs[c("value", "nomdf", "dendf", "p.value")]))

  two <- summary(iv_robust(mpg ~ hp + am | wt + gear, data = mtcars,
                           se_type = "classical", diagnostics = TRUE))
  tab2 <- build_ivreg_diagnostics_mat(two)
  fs2 <- two$diagnostic_first_stage_fstatistic
  expect_equal(rownames(tab2),
               c("Weak instruments (hp)", "Weak instruments (am)", "Wu-Hausman", "Sargan"))
  expect_equal(unname(tab2[1:2, "value"]), unname(fs2[c("hp:value", "am:value")]))
  expect_equal(unname(tab2[1:2, "p.value"]), unname(fs2[c("hp:p.value", "am:p.value")]))
  expect_false(anyNA(tab2[1:3, ]))
  expect_output(print(two), "Weak instruments \\(am\\)")
})

test_that("#389: glance() works with multiple endogenous regressors", {
  m <- iv_robust(mpg ~ hp + wt | am + cyl, data = mtcars, diagnostics = TRUE)
  g <- glance(m)
  expect_equal(nrow(g), 1L)
  # reports the weakest of the per-regressor first stages
  fs <- m$diagnostic_first_stage_fstatistic
  expect_equal(g$statistic.weakinst, unname(min(fs[grep("(^|:)value$", names(fs))])))
})

test_that("#397: model.frame() on iv_robust returns the model variables", {
  # Called directly: loading clubSandwich registers a model.frame method for
  # this class too, so bare dispatch depends on what the suite loaded first.
  m <- iv_robust(mpg ~ wt + hp | am + hp, data = mtcars)
  mf <- estimatr:::model.frame.iv_robust(m)
  expect_equal(nrow(mf), nrow(mtcars))
  expect_setequal(names(mf), c("mpg", "wt", "hp", "am"))
})
