library(estimatr)

# Numerical identity with estimatr 1.0.6 for all overlapping functionality.
#
# The reference is a recording rather than a live `estimatr::` call. Once this
# package takes the estimatr name, calling estimatr here is a dependency on
# itself; and a recording is the stronger test anyway, since it pins the answer
# the released package gave rather than whatever happens to be installed.
# `data-raw/make_estimatr_reference.R` produces it.

# Fields that must be bit-identical for all model types
INFERENCE_FIELDS <- c(
  "coefficients", "std.error", "statistic", "df",
  "p.value", "conf.low", "conf.high", "vcov"
)

# Additional fields that match for lm-family models
LM_FIELDS <- c(INFERENCE_FIELDS, "fitted.values", "res_var", "fstatistic", "df.residual", "nobs", "rank")

# Compare numeric fields between the recorded reference and a fresh fit.
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

dat <- ref_data_vs()

# ---- lm_robust: SE types ----

test_that("lm_robust HC0 identical to estimatr", {
  e0 <- ref("lmr_HC0")
  ez <- lm_robust(y ~ x + z, data = dat, se_type = "HC0")
  check_identical(e0, ez, LM_FIELDS, "HC0")
  check_identical(e0, ez, c("r.squared", "adj.r.squared", "tss"), "HC0-r2")
})

test_that("lm_robust HC1/stata identical to estimatr", {
  e0 <- ref("lmr_HC1")
  ez <- lm_robust(y ~ x + z, data = dat, se_type = "HC1")
  check_identical(e0, ez, LM_FIELDS, "HC1")
})

test_that("lm_robust HC2 (default) identical to estimatr", {
  e0 <- ref("lmr_HC2")
  ez <- lm_robust(y ~ x + z, data = dat)
  check_identical(e0, ez, LM_FIELDS, "HC2")
  check_identical(e0, ez, c("r.squared", "adj.r.squared", "tss"), "HC2-r2")
})

test_that("lm_robust HC3 identical to estimatr", {
  e0 <- ref("lmr_HC3")
  ez <- lm_robust(y ~ x + z, data = dat, se_type = "HC3")
  check_identical(e0, ez, LM_FIELDS, "HC3")
})

test_that("lm_robust classical identical to estimatr", {
  e0 <- ref("lmr_classical")
  ez <- lm_robust(y ~ x + z, data = dat, se_type = "classical")
  check_identical(e0, ez, LM_FIELDS, "classical")
  check_identical(e0, ez, c("r.squared", "adj.r.squared"), "classical-r2")
})

test_that("lm_robust stata identical to estimatr", {
  e0 <- ref("lmr_stata")
  ez <- lm_robust(y ~ x + z, data = dat, se_type = "stata")
  check_identical(e0, ez, LM_FIELDS, "stata")
})

# ---- lm_robust: clustered ----

test_that("lm_robust CR2 (default clustered) identical to estimatr", {
  e0 <- ref("lmr_CR2")
  ez <- lm_robust(y ~ x + z, data = dat, clusters = cl)
  check_identical(e0, ez, LM_FIELDS, "CR2")
  check_identical(e0, ez, c("r.squared", "adj.r.squared", "nclusters"), "CR2-r2")
})

test_that("lm_robust CR0 identical to estimatr", {
  e0 <- ref("lmr_CR0")
  ez <- lm_robust(y ~ x + z, data = dat, clusters = cl, se_type = "CR0")
  check_identical(e0, ez, LM_FIELDS, "CR0")
})

test_that("lm_robust clustered stata identical to estimatr", {
  e0 <- ref("lmr_cl_stata")
  ez <- lm_robust(y ~ x + z, data = dat, clusters = cl, se_type = "stata")
  check_identical(e0, ez, LM_FIELDS, "cl-stata")
})

# ---- lm_robust: weighted ----

test_that("lm_robust weighted: coefs, SEs, and R2 identical to estimatr", {
  e0 <- ref("lmr_w")
  ez <- lm_robust(y ~ x + z, data = dat, weights = w)
  check_identical(e0, ez, INFERENCE_FIELDS, "weighted")
  check_identical(e0, ez, c("fitted.values", "res_var", "df.residual", "nobs"), "weighted-aux")
  check_identical(e0, ez, c("r.squared", "adj.r.squared", "tss"), "weighted-r2")
})

test_that("lm_robust weighted + clustered: coefs and SEs identical to estimatr", {
  e0 <- ref("lmr_w_cl")
  ez <- lm_robust(y ~ x + z, data = dat, weights = w, clusters = cl)
  check_identical(e0, ez, INFERENCE_FIELDS, "weighted-clustered")
})

# ---- lm_robust: multivariate ----

test_that("lm_robust multivariate identical to estimatr", {
  e0 <- ref("lmr_mv")
  ez <- lm_robust(cbind(y, y2) ~ x + z, data = dat)
  check_identical(e0, ez, c("coefficients", "std.error", "vcov", "fitted.values", "df.residual"), "mv")
  check_identical(e0, ez, c("r.squared", "adj.r.squared"), "mv-r2")
})

# ---- lm_robust: no-intercept ----

test_that("lm_robust no-intercept identical to estimatr", {
  e0 <- ref("lmr_noint")
  ez <- lm_robust(y ~ 0 + x + z, data = dat)
  check_identical(e0, ez, LM_FIELDS, "no-int")
})

# ---- lm_lin ----

test_that("lm_lin unweighted identical to estimatr", {
  e0 <- ref("lin")
  ez <- lm_lin(y ~ z, covariates = ~ x, data = dat)
  check_identical(e0, ez, LM_FIELDS, "lm_lin")
  check_identical(e0, ez, c("r.squared", "adj.r.squared"), "lm_lin-r2")
})

test_that("lm_lin clustered identical to estimatr", {
  e0 <- ref("lin_cl")
  ez <- lm_lin(y ~ z, covariates = ~ x, data = dat, clusters = cl)
  check_identical(e0, ez, INFERENCE_FIELDS, "lm_lin-cl")
})

test_that("lm_lin weighted: coefs, SEs, and R2 identical to estimatr", {
  e0 <- ref("lin_w")
  ez <- lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)
  check_identical(e0, ez, INFERENCE_FIELDS, "lm_lin-w")
  check_identical(e0, ez, c("r.squared", "adj.r.squared", "tss"), "lm_lin-w-r2")
})

test_that("lm_lin scaled_center identical to estimatr", {
  e0 <- ref("lin")
  ez <- lm_lin(y ~ z, covariates = ~ x, data = dat)
  expect_equal(e0$scaled_center, ez$scaled_center, tolerance = 0)
})

test_that("lm_lin multi-covariate identical to estimatr", {
  e0 <- ref("lin_mv")
  ez <- lm_lin(y ~ z, covariates = ~ x + y2, data = dat)
  check_identical(e0, ez, INFERENCE_FIELDS, "lm_lin-mv")
})

# ---- iv_robust ----

test_that("iv_robust HC2 identical to estimatr", {
  e0 <- ref("iv_HC2")
  ez <- iv_robust(y ~ z | iv, data = dat)
  check_identical(e0, ez, INFERENCE_FIELDS, "iv")
  check_identical(e0, ez, c("r.squared", "adj.r.squared", "nobs"), "iv-r2")
})

test_that("iv_robust classical identical to estimatr", {
  e0 <- ref("iv_classical")
  ez <- iv_robust(y ~ z | iv, data = dat, se_type = "classical")
  check_identical(e0, ez, INFERENCE_FIELDS, "iv-classical")
})

test_that("iv_robust clustered CR2 identical to estimatr", {
  e0 <- ref("iv_CR2")
  ez <- iv_robust(y ~ z | iv, data = dat, clusters = cl)
  check_identical(e0, ez, INFERENCE_FIELDS, "iv-cr2")
})

test_that("iv_robust with diagnostics identical to estimatr", {
  e0 <- ref("iv_diag")
  ez <- iv_robust(y ~ z | iv, data = dat, diagnostics = TRUE)
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
  e0 <- ref("dim_std")
  ez <- difference_in_means(y ~ z, data = dat)
  check_identical(e0, ez, DIM_FIELDS, "dim-std")
  expect_equal(e0$design, ez$design)
})

test_that("difference_in_means clustered (CR2) identical to estimatr", {
  e0 <- ref("dim_cl")
  ez <- difference_in_means(y ~ cl_z, clusters = cl, data = dat)
  check_identical(e0, ez, DIM_FIELDS, "dim-cl")
  expect_equal(e0$design, ez$design)
})

test_that("difference_in_means blocked identical to estimatr", {
  e0 <- ref("dim_bl")
  ez <- difference_in_means(y ~ z_block, blocks = bl, data = dat)
  check_identical(e0, ez, DIM_FIELDS, "dim-bl")
  expect_equal(e0$design, ez$design)
  expect_equal(e0$nblocks, ez$nblocks)
})

test_that("difference_in_means matched-pairs identical to estimatr", {
  dp <- ref_data_pairs()
  e0 <- ref("dim_mp")
  ez <- difference_in_means(y ~ z, blocks = pr, data = dp)
  check_identical(e0, ez, DIM_FIELDS, "dim-mp")
  expect_equal(e0$design, ez$design)
})

test_that("difference_in_means weighted: coefs identical, delegates to lm_robust", {
  e0 <- ref("dim_w")
  ez <- difference_in_means(y ~ z, weights = w, data = dat)
  check_identical(e0, ez, c("coefficients", "std.error", "df", "p.value"), "dim-w")
})

# ---- lh_robust ----

test_that("lh_robust lm_robust component identical to estimatr", {
  e0 <- ref("lh")
  ez <- lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z + 2*x = 0")
  check_identical(e0$lm_robust, ez$lm_robust, INFERENCE_FIELDS, "lh-lmr")
})

test_that("lh_robust lh component identical to estimatr", {
  e0 <- ref("lh")
  ez <- lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z + 2*x = 0")
  for (nm in c("coefficients", "std.error", "statistic", "p.value", "conf.low", "conf.high")) {
    expect_equal(e0$lh[[nm]], ez$lh[[nm]], tolerance = 0, label = paste0("lh$", nm))
  }
})

# ---- subset and factor treatment ----

test_that("lm_robust subset identical to estimatr", {
  e0 <- ref("lmr_subset")
  ez <- lm_robust(y ~ x, data = dat, subset = z == 1)
  check_identical(e0, ez, LM_FIELDS, "subset")
})

test_that("lm_robust factor covariate identical to estimatr", {
  dat2 <- dat
  dat2$bl_f <- factor(dat2$bl)
  e0 <- ref("lmr_factor")
  ez <- lm_robust(y ~ z + bl_f, data = dat2)
  check_identical(e0, ez, c("coefficients", "std.error", "p.value"), "factor")
})

# ---- se_type = "none" ----

test_that("lm_robust se_type=none coefs identical to estimatr", {
  e0 <- ref("lmr_none")
  ez <- lm_robust(y ~ x + z, data = dat, se_type = "none")
  expect_equal(e0$coefficients, ez$coefficients, tolerance = 0)
  expect_equal(e0$fitted.values, ez$fitted.values, tolerance = 0)
})
