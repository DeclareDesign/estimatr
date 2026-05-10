library(estimatrZero)

# Helper: simulate observational data with a known propensity score model
sim_obs <- function(n, seed, tau = 2, beta_x = 0.5, ps_strength = 0.5) {
  set.seed(seed)
  x <- rnorm(n)
  ps <- plogis(ps_strength * x)
  z <- rbinom(n, 1, ps)
  y <- tau * z + beta_x * x + rnorm(n)
  w <- ifelse(z == 1, 1 / ps, 1 / (1 - ps))
  data.frame(y = y, z = z, x = x, ps = ps, w = w)
}

# ---- centering ----

test_that("lm_lin centers on weighted mean when weights supplied", {
  dat <- sim_obs(500, 1)
  m <- lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)
  expect_equal(
    unname(m$scaled_center),
    weighted.mean(dat$x, dat$w),
    tolerance = 1e-10
  )
})

test_that("lm_lin centers on simple mean when no weights", {
  dat <- sim_obs(500, 1)
  m <- lm_lin(y ~ z, covariates = ~ x, data = dat)
  expect_equal(
    unname(m$scaled_center),
    mean(dat$x),
    tolerance = 1e-10
  )
})

test_that("weighted centering differs from unweighted when weights are informative", {
  dat <- sim_obs(500, 2, ps_strength = 1.0)
  m_w <- lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)
  m_uw <- lm_lin(y ~ z, covariates = ~ x, data = dat)
  # Centering points should differ when propensity scores vary
  expect_false(isTRUE(all.equal(m_w$scaled_center, m_uw$scaled_center, tolerance = 1e-6)))
})

# ---- manual construction equivalence ----

test_that("lm_lin with weights matches manual weighted Lin construction", {
  dat <- sim_obs(300, 3)

  # Build Lin design matrix manually using weighted centering
  x_center <- weighted.mean(dat$x, dat$w)
  x_c <- dat$x - x_center
  dat$x_c <- x_c

  m_auto <- lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)
  m_manual <- lm_robust(y ~ z + x_c + z:x_c, data = dat, weights = w)

  expect_equal(
    unname(coef(m_auto)["z"]),
    unname(coef(m_manual)["z"]),
    tolerance = 1e-10
  )
  expect_equal(
    unname(m_auto$std.error["z"]),
    unname(m_manual$std.error["z"]),
    tolerance = 1e-10
  )
})

test_that("lm_lin without weights matches manual unweighted Lin construction", {
  dat <- sim_obs(300, 4)

  x_center <- mean(dat$x)
  x_c <- dat$x - x_center
  dat$x_c <- x_c

  m_auto <- lm_lin(y ~ z, covariates = ~ x, data = dat)
  m_manual <- lm_robust(y ~ z + x_c + z:x_c, data = dat)

  expect_equal(
    unname(coef(m_auto)["z"]),
    unname(coef(m_manual)["z"]),
    tolerance = 1e-10
  )
  expect_equal(
    unname(m_auto$std.error["z"]),
    unname(m_manual$std.error["z"]),
    tolerance = 1e-10
  )
})

# ---- constant weights ----

test_that("constant IPW weights give same estimate as unweighted", {
  # When all propensity scores are 0.5, all IPW weights are 2 (constant)
  set.seed(5)
  n <- 500
  dat <- data.frame(y = rnorm(n), z = rbinom(n, 1, 0.5), x = rnorm(n))
  dat$w <- rep(2, n)

  m_w <- lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)
  m_uw <- lm_lin(y ~ z, covariates = ~ x, data = dat)

  expect_equal(
    unname(coef(m_w)["z"]),
    unname(coef(m_uw)["z"]),
    tolerance = 1e-10
  )
})

test_that("constant weights give same SEs as unweighted", {
  set.seed(6)
  n <- 500
  dat <- data.frame(y = rnorm(n), z = rbinom(n, 1, 0.5), x = rnorm(n))
  dat$w <- rep(3, n)

  m_w <- lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)
  m_uw <- lm_lin(y ~ z, covariates = ~ x, data = dat)

  expect_equal(
    unname(m_w$std.error["z"]),
    unname(m_uw$std.error["z"]),
    tolerance = 1e-10
  )
})

# ---- diagnostics ----

test_that("IPW lm_lin SEs are finite and positive", {
  dat <- sim_obs(500, 7)
  m <- lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)
  expect_true(all(is.finite(m$std.error)))
  expect_true(all(m$std.error > 0))
})

test_that("IPW lm_lin R-squared is in [0, 1]", {
  dat <- sim_obs(500, 8)
  m <- lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)
  expect_gte(m$r.squared, 0)
  expect_lte(m$r.squared, 1.0001)
})

test_that("IPW lm_lin p-values are in (0, 1)", {
  dat <- sim_obs(1000, 9, tau = 2)
  m <- lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)
  expect_true(all(m$p.value > 0, na.rm = TRUE))
  expect_true(all(m$p.value < 1, na.rm = TRUE))
  # With tau = 2 and n = 1000, treatment should be significant
  expect_lt(m$p.value["z"], 0.05)
})

test_that("IPW lm_lin confidence intervals cover true ATE in large samples", {
  # With n = 5000, CI should comfortably cover true ATE of 2
  dat <- sim_obs(5000, 10, tau = 2.0, ps_strength = 0.5)
  m <- lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)
  expect_gte(m$conf.low["z"], 1.7)
  expect_lte(m$conf.high["z"], 2.3)
})

# ---- weights + clusters ----

test_that("lm_lin with weights and clusters runs without error", {
  set.seed(11)
  n <- 200
  dat <- data.frame(
    y = rnorm(n),
    z = rbinom(n, 1, 0.5),
    x = rnorm(n),
    cl = rep(1:20, 10)
  )
  dat$w <- runif(n, 0.5, 2)

  m <- lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w, clusters = cl)
  expect_s3_class(m, "lm_robust")
  expect_equal(m$se_type, "CR2")
  expect_true(all(is.finite(m$std.error)))
})

# ---- multiple covariates ----

test_that("IPW lm_lin works with multiple covariates", {
  set.seed(12)
  n <- 500
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  ps <- plogis(0.3 * x1 + 0.3 * x2)
  z <- rbinom(n, 1, ps)
  y <- 2 * z + 0.5 * x1 + 0.3 * x2 + rnorm(n)
  w <- ifelse(z == 1, 1 / ps, 1 / (1 - ps))
  dat <- data.frame(y = y, z = z, x1 = x1, x2 = x2, w = w)

  m <- lm_lin(y ~ z, covariates = ~ x1 + x2, data = dat, weights = w)
  expect_equal(length(m$scaled_center), 2L)
  expect_equal(
    unname(m$scaled_center["x1"]),
    weighted.mean(x1, w),
    tolerance = 1e-10
  )
  expect_equal(
    unname(m$scaled_center["x2"]),
    weighted.mean(x2, w),
    tolerance = 1e-10
  )
  expect_true(all(is.finite(m$std.error)))
})

# ---- weighted R2 is always in [0, 1] ----

test_that("weighted R2 is in [0, 1] when weights are positively correlated with Y", {
  set.seed(13)
  n <- 200
  x <- 1:n / n
  y <- x + rnorm(n, 0, 0.1)
  w_ipw <- x + 0.1
  dat <- data.frame(y = y, x = x, z = rbinom(n, 1, 0.5), w = w_ipw)

  m <- lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)
  expect_gte(m$r.squared, 0)
  expect_lte(m$r.squared, 1)
})
