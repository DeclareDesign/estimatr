library(estimatrZero)

# These tests confirm numerical identity with estimatr 1.0.6 for all
# overlapping functionality. Where results intentionally differ (weighted
# R-squared bug fix), that divergence is tested explicitly.

skip_if_not_installed("estimatr")

# Fields that must be bit-identical for all model types
INFERENCE_FIELDS <- c(
  "coefficients", "std.error", "statistic", "df",
  "p.value", "conf.low", "conf.high", "vcov"
)

# Additional fields that match for lm-family models
LM_FIELDS <- c(INFERENCE_FIELDS, "fitted.values", "res_var", "fstatistic", "df.residual", "nobs", "rank")

# Compare numeric fields between two model objects.
# tolerance = 1e-12: accepts machine-epsilon accumulation from large design
# matrices and CR2 eigenvalue computations (~5e-15 observed) while catching
# any real numerical discrepancy (would be >= 1e-10 for meaningful differences).
check_identical <- function(e0, ez, fields, label = "") {
  for (nm in fields) {
    a <- e0[[nm]]
    b <- ez[[nm]]
    if (is.null(a) || is.null(b)) next
    info_msg <- if (nzchar(label)) paste0("[", label, "] ", nm) else nm
    expect_equal(a, b, tolerance = 1e-12, label = info_msg)
  }
}

# ---- shared data ----

local({
  set.seed(42)
  n <- 200
  assign("dat", data.frame(
    y  = rnorm(n),
    y2 = rnorm(n),
    x  = rnorm(n),
    z  = rbinom(n, 1, 0.5),
    cl = rep(1:20, 10),
    bl = rep(1:20, each = 10),
    w  = runif(n, 0.5, 2.0),
    iv = rnorm(n) + rnorm(n, 0, 0.3)
  ), envir = parent.env(environment()))
  dat$z_block <<- rep(rep(c(0L, 1L), 5L), 20L)
  dat$cl_z    <<- as.integer(dat$cl %% 2 == 0)
})

# ---- lm_robust: SE types ----

test_that("lm_robust HC0 identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ x + z, data = dat, se_type = "HC0")
  ez <- estimatrZero::lm_robust(y ~ x + z, data = dat, se_type = "HC0")
  check_identical(e0, ez, LM_FIELDS, "HC0")
  check_identical(e0, ez, c("r.squared", "adj.r.squared", "tss"), "HC0-r2")
})

test_that("lm_robust HC1/stata identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ x + z, data = dat, se_type = "HC1")
  ez <- estimatrZero::lm_robust(y ~ x + z, data = dat, se_type = "HC1")
  check_identical(e0, ez, LM_FIELDS, "HC1")
})

test_that("lm_robust HC2 (default) identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ x + z, data = dat)
  ez <- estimatrZero::lm_robust(y ~ x + z, data = dat)
  check_identical(e0, ez, LM_FIELDS, "HC2")
  check_identical(e0, ez, c("r.squared", "adj.r.squared", "tss"), "HC2-r2")
})

test_that("lm_robust HC3 identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ x + z, data = dat, se_type = "HC3")
  ez <- estimatrZero::lm_robust(y ~ x + z, data = dat, se_type = "HC3")
  check_identical(e0, ez, LM_FIELDS, "HC3")
})

test_that("lm_robust classical identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ x + z, data = dat, se_type = "classical")
  ez <- estimatrZero::lm_robust(y ~ x + z, data = dat, se_type = "classical")
  check_identical(e0, ez, LM_FIELDS, "classical")
  check_identical(e0, ez, c("r.squared", "adj.r.squared"), "classical-r2")
})

test_that("lm_robust stata identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ x + z, data = dat, se_type = "stata")
  ez <- estimatrZero::lm_robust(y ~ x + z, data = dat, se_type = "stata")
  check_identical(e0, ez, LM_FIELDS, "stata")
})

# ---- lm_robust: clustered ----

test_that("lm_robust CR2 (default clustered) identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ x + z, data = dat, clusters = cl)
  ez <- estimatrZero::lm_robust(y ~ x + z, data = dat, clusters = cl)
  check_identical(e0, ez, LM_FIELDS, "CR2")
  check_identical(e0, ez, c("r.squared", "adj.r.squared", "nclusters"), "CR2-r2")
})

test_that("lm_robust CR0 identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ x + z, data = dat, clusters = cl, se_type = "CR0")
  ez <- estimatrZero::lm_robust(y ~ x + z, data = dat, clusters = cl, se_type = "CR0")
  check_identical(e0, ez, LM_FIELDS, "CR0")
})

test_that("lm_robust clustered stata identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ x + z, data = dat, clusters = cl, se_type = "stata")
  ez <- estimatrZero::lm_robust(y ~ x + z, data = dat, clusters = cl, se_type = "stata")
  check_identical(e0, ez, LM_FIELDS, "cl-stata")
})

# ---- lm_robust: weighted (coefs/SEs identical; R2 intentionally differs) ----

test_that("lm_robust weighted: coefs and SEs identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ x + z, data = dat, weights = w)
  ez <- estimatrZero::lm_robust(y ~ x + z, data = dat, weights = w)
  check_identical(e0, ez, INFERENCE_FIELDS, "weighted")
  check_identical(e0, ez, c("fitted.values", "res_var", "df.residual", "nobs"), "weighted-aux")
})

test_that("lm_robust weighted: R2 differs from estimatr (bug fix)", {
  # estimatr used weights^2 in TSS; estimatrZero uses weights
  e0 <- estimatr::lm_robust(y ~ x + z, data = dat, weights = w)
  ez <- estimatrZero::lm_robust(y ~ x + z, data = dat, weights = w)
  expect_false(isTRUE(all.equal(e0$tss, ez$tss)))
  expect_false(isTRUE(all.equal(e0$r.squared, ez$r.squared)))
  # estimatrZero R2 is in [0, 1]; estimatr's may not be
  expect_gte(ez$r.squared, 0)
  expect_lte(ez$r.squared, 1.0001)
})

test_that("lm_robust weighted + clustered: coefs and SEs identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ x + z, data = dat, weights = w, clusters = cl)
  ez <- estimatrZero::lm_robust(y ~ x + z, data = dat, weights = w, clusters = cl)
  check_identical(e0, ez, INFERENCE_FIELDS, "weighted-clustered")
})

# ---- lm_robust: multivariate ----

test_that("lm_robust multivariate identical to estimatr", {
  e0 <- estimatr::lm_robust(cbind(y, y2) ~ x + z, data = dat)
  ez <- estimatrZero::lm_robust(cbind(y, y2) ~ x + z, data = dat)
  check_identical(e0, ez, c("coefficients", "std.error", "vcov", "fitted.values", "df.residual"), "mv")
  check_identical(e0, ez, c("r.squared", "adj.r.squared"), "mv-r2")
})

# ---- lm_robust: no-intercept ----

test_that("lm_robust no-intercept identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ 0 + x + z, data = dat)
  ez <- estimatrZero::lm_robust(y ~ 0 + x + z, data = dat)
  check_identical(e0, ez, LM_FIELDS, "no-int")
})

# ---- lm_lin ----

test_that("lm_lin unweighted identical to estimatr", {
  e0 <- estimatr::lm_lin(y ~ z, covariates = ~ x, data = dat)
  ez <- estimatrZero::lm_lin(y ~ z, covariates = ~ x, data = dat)
  check_identical(e0, ez, LM_FIELDS, "lm_lin")
  check_identical(e0, ez, c("r.squared", "adj.r.squared"), "lm_lin-r2")
})

test_that("lm_lin clustered identical to estimatr", {
  e0 <- estimatr::lm_lin(y ~ z, covariates = ~ x, data = dat, clusters = cl)
  ez <- estimatrZero::lm_lin(y ~ z, covariates = ~ x, data = dat, clusters = cl)
  check_identical(e0, ez, INFERENCE_FIELDS, "lm_lin-cl")
})

test_that("lm_lin weighted: coefs and SEs identical, R2 differs", {
  e0 <- estimatr::lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)
  ez <- estimatrZero::lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)
  check_identical(e0, ez, INFERENCE_FIELDS, "lm_lin-w")
  expect_false(isTRUE(all.equal(e0$tss, ez$tss)))
})

test_that("lm_lin scaled_center identical to estimatr", {
  e0 <- estimatr::lm_lin(y ~ z, covariates = ~ x, data = dat)
  ez <- estimatrZero::lm_lin(y ~ z, covariates = ~ x, data = dat)
  expect_equal(e0$scaled_center, ez$scaled_center, tolerance = 0)
})

test_that("lm_lin multi-covariate identical to estimatr", {
  e0 <- estimatr::lm_lin(y ~ z, covariates = ~ x + y2, data = dat)
  ez <- estimatrZero::lm_lin(y ~ z, covariates = ~ x + y2, data = dat)
  check_identical(e0, ez, INFERENCE_FIELDS, "lm_lin-mv")
})

# ---- iv_robust ----

test_that("iv_robust HC2 identical to estimatr", {
  e0 <- estimatr::iv_robust(y ~ z | iv, data = dat)
  ez <- estimatrZero::iv_robust(y ~ z | iv, data = dat)
  check_identical(e0, ez, INFERENCE_FIELDS, "iv")
  check_identical(e0, ez, c("r.squared", "adj.r.squared", "nobs"), "iv-r2")
})

test_that("iv_robust classical identical to estimatr", {
  e0 <- estimatr::iv_robust(y ~ z | iv, data = dat, se_type = "classical")
  ez <- estimatrZero::iv_robust(y ~ z | iv, data = dat, se_type = "classical")
  check_identical(e0, ez, INFERENCE_FIELDS, "iv-classical")
})

test_that("iv_robust clustered CR2 identical to estimatr", {
  e0 <- estimatr::iv_robust(y ~ z | iv, data = dat, clusters = cl)
  ez <- estimatrZero::iv_robust(y ~ z | iv, data = dat, clusters = cl)
  check_identical(e0, ez, INFERENCE_FIELDS, "iv-cr2")
})

test_that("iv_robust with diagnostics identical to estimatr", {
  e0 <- estimatr::iv_robust(y ~ z | iv, data = dat, diagnostics = TRUE)
  ez <- estimatrZero::iv_robust(y ~ z | iv, data = dat, diagnostics = TRUE)
  check_identical(e0, ez, INFERENCE_FIELDS, "iv-diag")
  # Diagnostic test statistics should match
  expect_equal(
    e0$diagnostic_first_stage_fstatistic,
    ez$diagnostic_first_stage_fstatistic,
    tolerance = 0
  )
  expect_equal(
    e0$diagnostic_endogeneity_test,
    ez$diagnostic_endogeneity_test,
    tolerance = 0
  )
})

# ---- difference_in_means ----

DIM_FIELDS <- c("coefficients", "std.error", "df", "statistic", "p.value", "conf.low", "conf.high", "nobs")

test_that("difference_in_means standard (Welch) identical to estimatr", {
  e0 <- estimatr::difference_in_means(y ~ z, data = dat)
  ez <- estimatrZero::difference_in_means(y ~ z, data = dat)
  check_identical(e0, ez, DIM_FIELDS, "dim-std")
  expect_equal(e0$design, ez$design)
})

test_that("difference_in_means clustered (CR2) identical to estimatr", {
  e0 <- estimatr::difference_in_means(y ~ cl_z, clusters = cl, data = dat)
  ez <- estimatrZero::difference_in_means(y ~ cl_z, clusters = cl, data = dat)
  check_identical(e0, ez, DIM_FIELDS, "dim-cl")
  expect_equal(e0$design, ez$design)
})

test_that("difference_in_means blocked identical to estimatr", {
  e0 <- estimatr::difference_in_means(y ~ z_block, blocks = bl, data = dat)
  ez <- estimatrZero::difference_in_means(y ~ z_block, blocks = bl, data = dat)
  check_identical(e0, ez, DIM_FIELDS, "dim-bl")
  expect_equal(e0$design, ez$design)
  expect_equal(e0$nblocks, ez$nblocks)
})

test_that("difference_in_means matched-pairs identical to estimatr", {
  # Paired design: 2 units per block
  set.seed(99)
  n_pairs <- 50
  dp <- data.frame(
    y  = rnorm(n_pairs * 2),
    z  = rep(c(0L, 1L), n_pairs),
    pr = rep(seq_len(n_pairs), each = 2)
  )
  e0 <- estimatr::difference_in_means(y ~ z, blocks = pr, data = dp)
  ez <- estimatrZero::difference_in_means(y ~ z, blocks = pr, data = dp)
  check_identical(e0, ez, DIM_FIELDS, "dim-mp")
  expect_equal(e0$design, ez$design)
})

test_that("difference_in_means weighted: coefs identical, delegates to lm_robust", {
  e0 <- estimatr::difference_in_means(y ~ z, weights = w, data = dat)
  ez <- estimatrZero::difference_in_means(y ~ z, weights = w, data = dat)
  check_identical(e0, ez, c("coefficients", "std.error", "df", "p.value"), "dim-w")
})

# ---- lh_robust ----

test_that("lh_robust lm_robust component identical to estimatr", {
  e0 <- estimatr::lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z + 2*x = 0")
  ez <- estimatrZero::lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z + 2*x = 0")
  check_identical(e0$lm_robust, ez$lm_robust, INFERENCE_FIELDS, "lh-lmr")
})

test_that("lh_robust lh component identical to estimatr", {
  e0 <- estimatr::lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z + 2*x = 0")
  ez <- estimatrZero::lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z + 2*x = 0")
  for (nm in c("coefficients", "std.error", "statistic", "p.value", "conf.low", "conf.high")) {
    expect_equal(e0$lh[[nm]], ez$lh[[nm]], tolerance = 0, label = paste0("lh$", nm))
  }
})

# ---- subset and factor treatment ----

test_that("lm_robust subset identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ x, data = dat, subset = z == 1)
  ez <- estimatrZero::lm_robust(y ~ x, data = dat, subset = z == 1)
  check_identical(e0, ez, LM_FIELDS, "subset")
})

test_that("lm_robust factor covariate identical to estimatr", {
  dat2 <- dat
  dat2$bl_f <- factor(dat2$bl)
  e0 <- estimatr::lm_robust(y ~ z + bl_f, data = dat2)
  ez <- estimatrZero::lm_robust(y ~ z + bl_f, data = dat2)
  check_identical(e0, ez, c("coefficients", "std.error", "p.value"), "factor")
})

# ---- se_type = "none" ----

test_that("lm_robust se_type=none coefs identical to estimatr", {
  e0 <- estimatr::lm_robust(y ~ x + z, data = dat, se_type = "none")
  ez <- estimatrZero::lm_robust(y ~ x + z, data = dat, se_type = "none")
  expect_equal(e0$coefficients, ez$coefficients, tolerance = 0)
  expect_equal(e0$fitted.values, ez$fitted.values, tolerance = 0)
})
