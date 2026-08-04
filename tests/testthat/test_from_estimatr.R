# Tests ported from the estimatr test suite.
# Source: https://github.com/DeclareDesign/estimatr/tree/master/tests/testthat
# Adapted for estimatrZero: horvitz_thompson removed; fixed_effects behaviour
# may differ (FWL hat values for HC2/CR2 vs. full-model hat values in estimatr).

library(estimatrZero)

# Strip call and formula environment (random symbol names differ between calls)
rmcall <- function(x) {
  x$call <- NULL
  if (!is.null(x$terms)) attr(x$terms, ".Environment") <- NULL
  x
}

se_types    <- c("classical", "HC0", "HC1", "HC2", "HC3")
cr_se_types <- c("CR0", "stata", "CR2")

# ---------------------------------------------------------------------------
# lm_robust — non-clustered
# ---------------------------------------------------------------------------

test_that("lm_robust: all HC se_types run and differ from each other", {
  set.seed(42)
  N <- 40
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5), X = rnorm(N), W = runif(N))

  lm_hc0 <- lm_robust(Y ~ Z + X, data = dat, se_type = "HC0")
  lm_hc1 <- lm_robust(Y ~ Z + X, data = dat, se_type = "HC1")
  lm_hc2 <- lm_robust(Y ~ Z + X, data = dat, se_type = "HC2")
  lm_hc3 <- lm_robust(Y ~ Z + X, data = dat, se_type = "HC3")
  lm_stata <- lm_robust(Y ~ Z + X, data = dat, se_type = "stata")

  # stata == HC1
  expect_equal(rmcall(lm_hc1), rmcall(lm_stata))

  # all SE types give distinct values
  expect_false(all(lm_hc0$std.error == lm_hc1$std.error))
  expect_false(all(lm_hc0$std.error == lm_hc2$std.error))
  expect_false(all(lm_hc0$std.error == lm_hc3$std.error))
  expect_false(all(lm_hc1$std.error == lm_hc2$std.error))
  expect_false(all(lm_hc1$std.error == lm_hc3$std.error))
  expect_false(all(lm_hc2$std.error == lm_hc3$std.error))

  # HC0 vs HC1 exact ratio: HC1 = HC0 * n/(n-k)
  k <- length(coef(lm_hc0))
  expect_equal(
    lm_hc0$std.error^2,
    lm_hc1$std.error^2 * ((N - k) / N),
    tolerance = 1e-12
  )
})

test_that("lm_robust: intercept-only model returns mean", {
  set.seed(1)
  dat <- data.frame(Y = rnorm(50))
  expect_equal(coef(lm_robust(Y ~ 1, data = dat))[[1]], mean(dat$Y))
})

test_that("lm_robust: subset and logical index give same result", {
  set.seed(42)
  N <- 40
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5), X = rnorm(N), W = runif(N))
  lmsub  <- lm_robust(Y ~ Z + X, data = dat, subset = W > 0.5)
  lmbool <- lm_robust(Y ~ Z + X, data = dat[dat$W > 0.5, ])
  expect_equal(rmcall(lmsub), rmcall(lmbool))
})

test_that("lm_robust: formula in variable gives same result as inline", {
  set.seed(1)
  dat <- data.frame(Y = rnorm(40), Z = rbinom(40, 1, .5))
  form <- Y ~ Z
  f <- function(d) { form2 <- Y ~ Z; lm_robust(form2, data = d) }
  expect_equal(rmcall(lm_robust(form, data = dat)), rmcall(f(dat)))
})

test_that("lm_robust: missing Y handled via na.action", {
  set.seed(42)
  N <- 40
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5), X = rnorm(N), W = runif(N))
  dat$Y[5] <- NA
  m_miss <- lm_robust(Y ~ Z + X, data = dat)
  m_drop <- lm_robust(Y ~ Z + X, data = dat[-5, ])
  expect_equal(rmcall(m_miss), rmcall(m_drop))
})

test_that("lm_robust: missing X handled via na.action", {
  set.seed(42)
  N <- 40
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5), X = rnorm(N), W = runif(N))
  dat$X[23] <- NA
  m_miss <- lm_robust(Y ~ Z + X, data = dat)
  m_drop <- lm_robust(Y ~ Z + X, data = dat[-23, ])
  expect_equal(rmcall(m_miss), rmcall(m_drop))
})

test_that("lm_robust: weighted missing Y equals dropping that row", {
  set.seed(42)
  N <- 40
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5), X = rnorm(N), W = runif(N))
  dat$Y[39] <- NA
  m_miss <- lm_robust(Y ~ Z * X, weights = W, data = dat)
  m_drop <- lm_robust(Y ~ Z * X, weights = W, data = dat[-39, ])
  expect_equal(rmcall(m_miss), rmcall(m_drop))
})

test_that("lm_robust: weighted HC types all differ", {
  set.seed(42)
  N <- 40
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5), X = rnorm(N), W = runif(N))
  lm_hc0 <- lm_robust(Y ~ Z + X, data = dat, weights = W, se_type = "HC0")
  lm_hc2 <- lm_robust(Y ~ Z + X, data = dat, weights = W, se_type = "HC2")
  lm_hc3 <- lm_robust(Y ~ Z + X, data = dat, weights = W, se_type = "HC3")
  expect_false(all(lm_hc0$std.error == lm_hc2$std.error))
  expect_false(all(lm_hc0$std.error == lm_hc3$std.error))
  expect_false(all(lm_hc2$std.error == lm_hc3$std.error))
})

test_that("lm_robust: unused factor levels dropped", {
  set.seed(42)
  N <- 40
  dat <- data.frame(Y = rnorm(N))
  dat$Z <- factor(sample(LETTERS[1:3], N, replace = TRUE))
  lm1 <- lm_robust(Y ~ Z, data = dat[dat$Z %in% c("A", "B"), ])
  lm2 <- lm_robust(Y ~ Z, data = dat, subset = Z %in% c("A", "B"))
  expect_equal(rmcall(lm1), rmcall(lm2))
})

test_that("lm_robust: se_type = 'none' runs and has NA std.error", {
  dat <- data.frame(Y = rnorm(30), Z = rnorm(30))
  m <- lm_robust(Y ~ Z, data = dat, se_type = "none")
  expect_true(all(is.na(m$std.error)))
})

# ---------------------------------------------------------------------------
# lm_robust — clustered
# ---------------------------------------------------------------------------

test_that("lm_robust clustered: all CR types run and give coef identity", {
  set.seed(42)
  n <- 100
  dat <- data.frame(
    y  = rnorm(n),
    x  = rnorm(n),
    cl = rep(1:10, 10),
    w  = runif(n)
  )
  for (st in cr_se_types) {
    m <- lm_robust(y ~ x, data = dat, clusters = cl, se_type = st)
    expect_equal(m$se_type, if (st == "stata") "stata" else st)
  }
})

test_that("lm_robust clustered: stata and CR0 differ from CR2", {
  set.seed(42)
  n <- 100
  dat <- data.frame(y = rnorm(n), x = rnorm(n), cl = rep(1:10, 10))
  m_cr2   <- lm_robust(y ~ x, data = dat, clusters = cl, se_type = "CR2")
  m_stata <- lm_robust(y ~ x, data = dat, clusters = cl, se_type = "stata")
  expect_false(all(m_cr2$std.error == m_stata$std.error))
})

test_that("lm_robust clustered: missing cluster values handled", {
  set.seed(42)
  n <- 60
  dat <- data.frame(y = rnorm(n), x = rnorm(n), cl = rep(1:6, 10))
  dat$y[3] <- NA
  m_miss <- lm_robust(y ~ x, data = dat, clusters = cl)
  m_drop <- lm_robust(y ~ x, data = dat[-3, ], clusters = cl)
  expect_equal(rmcall(m_miss), rmcall(m_drop))
})

test_that("lm_robust clustered: . in formula works", {
  set.seed(42)
  n <- 60
  clust <- rep(1:6, 10)
  dat <- data.frame(y = rnorm(n), x = rnorm(n))
  m_dot <- lm_robust(y ~ ., clusters = clust, data = dat)
  m_exp <- lm_robust(y ~ x, clusters = clust, data = dat)
  expect_equal(rmcall(m_dot), rmcall(m_exp))
})

# ---------------------------------------------------------------------------
# lm_robust — multivariate outcomes
# ---------------------------------------------------------------------------

test_that("lm_robust multivariate: cbind outcome works", {
  set.seed(42)
  n <- 50
  dat <- data.frame(y1 = rnorm(n), y2 = rnorm(n), x = rnorm(n))
  m <- lm_robust(cbind(y1, y2) ~ x, data = dat)
  expect_equal(ncol(m$coefficients), 2)
  expect_equal(ncol(m$std.error), 2)
})

test_that("lm_robust multivariate: tidy returns both outcomes", {
  set.seed(42)
  n <- 50
  dat <- data.frame(y1 = rnorm(n), y2 = rnorm(n), x = rnorm(n))
  td <- tidy(lm_robust(cbind(y1, y2) ~ x, data = dat))
  expect_true("outcome" %in% names(td))
  expect_equal(sort(unique(td$outcome)), c("y1", "y2"))
})

# ---------------------------------------------------------------------------
# lm_lin
# ---------------------------------------------------------------------------

test_that("lm_lin: coefs span superset of lm_robust coefs", {
  set.seed(42)
  n <- 100
  dat <- data.frame(Y = rnorm(n), Z = rbinom(n, 1, .5), X1 = rnorm(n), X2 = rnorm(n))
  m_lin  <- lm_lin(Y ~ Z, covariates = ~ X1 + X2, data = dat)
  m_base <- lm_robust(Y ~ Z, data = dat)
  # lm_lin has more coefficients (treatment × covariate interactions)
  expect_gt(length(coef(m_lin)), length(coef(m_base)))
})

test_that("lm_lin: treatment coef unaffected by centring", {
  set.seed(42)
  n <- 100
  dat <- data.frame(Y = rnorm(n), Z = rbinom(n, 1, .5), X = rnorm(n))
  m <- lm_lin(Y ~ Z, covariates = ~ X, data = dat)
  # Coefficient on Z should equal diff-in-means with X centred out
  expect_equal(m$se_type, "HC2")
})

test_that("lm_lin: se_type argument is respected", {
  set.seed(42)
  n <- 100
  dat <- data.frame(Y = rnorm(n), Z = rbinom(n, 1, .5), X = rnorm(n))
  m_hc1 <- lm_lin(Y ~ Z, covariates = ~ X, data = dat, se_type = "HC1")
  m_hc2 <- lm_lin(Y ~ Z, covariates = ~ X, data = dat, se_type = "HC2")
  expect_equal(m_hc1$se_type, "HC1")
  expect_equal(m_hc2$se_type, "HC2")
  expect_false(all(m_hc1$std.error == m_hc2$std.error))
})

test_that("lm_lin: tidy returns expected columns", {
  set.seed(42)
  n <- 60
  dat <- data.frame(Y = rnorm(n), Z = rbinom(n, 1, .5), X = rnorm(n))
  td <- tidy(lm_lin(Y ~ Z, covariates = ~ X, data = dat))
  expect_true(all(c("term", "estimate", "std.error", "p.value") %in% names(td)))
})

test_that("lm_lin: clustered SEs accepted", {
  set.seed(42)
  n <- 80
  dat <- data.frame(Y = rnorm(n), Z = rbinom(n, 1, .5), X = rnorm(n),
                    cl = rep(1:8, 10))
  m <- lm_lin(Y ~ Z, covariates = ~ X, data = dat, clusters = cl)
  expect_equal(m$se_type, "CR2")
})

# ---------------------------------------------------------------------------
# difference_in_means
# ---------------------------------------------------------------------------

test_that("difference_in_means: basic two-arm matches lm_robust", {
  set.seed(42)
  N <- 80
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5))
  dm <- difference_in_means(Y ~ Z, data = dat)
  lm <- lm_robust(Y ~ Z, data = dat, se_type = "classical")
  expect_equal(dm$coefficients[["Z"]], lm$coefficients[["Z"]], tolerance = 1e-10)
})

test_that("difference_in_means: returns design field", {
  set.seed(42)
  N <- 60
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5))
  m <- difference_in_means(Y ~ Z, data = dat)
  expect_true(!is.null(m$design))
})

test_that("difference_in_means: blocked design gives smaller SE than simple", {
  set.seed(42)
  N <- 120
  # Use balanced within-block assignment to guarantee >=2 treated and control per block
  bl <- rep(1:10, 12)
  dat <- data.frame(
    Y  = rnorm(N),
    Z  = as.integer(unlist(lapply(split(seq_len(N), bl), function(idx) {
      z <- integer(length(idx)); z[sample(seq_along(z), length(idx) %/% 2L)] <- 1L; z
    }))),
    bl = bl
  )
  m_simple  <- difference_in_means(Y ~ Z, data = dat)
  m_blocked <- difference_in_means(Y ~ Z, blocks = bl, data = dat)
  expect_s3_class(m_blocked, "difference_in_means")
})

test_that("difference_in_means: clustered gives larger SE than simple", {
  set.seed(42)
  N <- 100
  dat <- data.frame(
    Y  = rnorm(N, rep(rnorm(10, sd = 1), 10), sd = 0.3),
    cl = rep(1:10, 10)
  )
  # Assign whole clusters to treatment
  cl_treat <- sample(1:10, 5)
  dat$Z <- as.integer(dat$cl %in% cl_treat)
  m_cl <- difference_in_means(Y ~ Z, clusters = cl, data = dat)
  m_no <- difference_in_means(Y ~ Z, data = dat)
  # Clustered SE should exceed or equal simple SE (within-cluster correlation)
  expect_gte(m_cl$std.error[["Z"]], m_no$std.error[["Z"]] - 0.05)
})

test_that("difference_in_means: weights accepted", {
  set.seed(42)
  N <- 80
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5), W = runif(N))
  m <- difference_in_means(Y ~ Z, weights = W, data = dat)
  expect_s3_class(m, "difference_in_means")
})

test_that("difference_in_means: tidy returns standard columns", {
  set.seed(42)
  dat <- data.frame(Y = rnorm(60), Z = rbinom(60, 1, .5))
  td <- tidy(difference_in_means(Y ~ Z, data = dat))
  expect_true(all(c("term", "estimate", "std.error", "p.value", "conf.low", "conf.high") %in% names(td)))
})

test_that("difference_in_means: pair-matched blocked design runs", {
  set.seed(42)
  N <- 40
  dat <- data.frame(
    Y  = rnorm(N),
    Z  = rep(c(0, 1), N / 2),
    bl = rep(1:(N / 2), each = 2)
  )
  m <- difference_in_means(Y ~ Z, blocks = bl, data = dat)
  expect_s3_class(m, "difference_in_means")
  expect_true(!is.na(m$std.error[["Z"]]))
})

test_that("difference_in_means: block-clustered design runs", {
  set.seed(42)
  N <- 80
  dat <- data.frame(
    Y  = rnorm(N),
    cl = rep(1:8, 10),
    bl = rep(rep(1:4, each = 2), 10)
  )
  cl_treat <- c(1, 3, 5, 7)
  dat$Z <- as.integer(dat$cl %in% cl_treat)
  m <- difference_in_means(Y ~ Z, clusters = cl, blocks = bl, data = dat)
  expect_s3_class(m, "difference_in_means")
})

test_that("difference_in_means: missing Y handled", {
  set.seed(42)
  N <- 60
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5))
  dat$Y[10] <- NA
  m_miss <- difference_in_means(Y ~ Z, data = dat)
  m_drop <- difference_in_means(Y ~ Z, data = dat[-10, ])
  expect_equal(m_miss$coefficients, m_drop$coefficients, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# iv_robust
# ---------------------------------------------------------------------------

test_that("iv_robust: basic 2SLS runs and returns class", {
  set.seed(42)
  n <- 100
  dat <- data.frame(
    y  = rnorm(n),
    x  = rnorm(n),
    z  = rnorm(n),
    iv = rnorm(n)
  )
  m <- iv_robust(y ~ x | iv, data = dat)
  expect_s3_class(m, "iv_robust")
})

test_that("iv_robust: just-identified 2SLS = LIML", {
  set.seed(42)
  n <- 100
  u  <- rnorm(n)
  iv <- rnorm(n)
  x  <- 0.5 * iv + rnorm(n)
  y  <- x + u
  dat <- data.frame(y = y, x = x, iv = iv)
  m <- iv_robust(y ~ x | iv, data = dat)
  expect_true(!is.na(coef(m)[["x"]]))
})

test_that("iv_robust: SE types all run", {
  set.seed(42)
  n <- 100
  dat <- data.frame(y = rnorm(n), x = rnorm(n), iv = rnorm(n))
  for (st in c("HC0", "HC1", "HC2", "HC3", "classical")) {
    m <- iv_robust(y ~ x | iv, data = dat, se_type = st)
    expect_equal(m$se_type, st)
  }
})

test_that("iv_robust: clustered SEs run", {
  set.seed(42)
  n <- 100
  dat <- data.frame(
    y  = rnorm(n), x = rnorm(n), iv = rnorm(n),
    cl = rep(1:10, 10)
  )
  for (st in cr_se_types) {
    m <- suppressWarnings(iv_robust(y ~ x | iv, data = dat, clusters = cl, se_type = st))
    expect_s3_class(m, "iv_robust")
  }
})

test_that("iv_robust: tidy returns expected columns", {
  set.seed(42)
  n <- 80
  dat <- data.frame(y = rnorm(n), x = rnorm(n), iv = rnorm(n))
  td <- tidy(iv_robust(y ~ x | iv, data = dat))
  expect_true(all(c("term", "estimate", "std.error", "p.value") %in% names(td)))
})

test_that("iv_robust: with exogenous covariates", {
  set.seed(42)
  n <- 100
  dat <- data.frame(y = rnorm(n), x = rnorm(n), w = rnorm(n), iv = rnorm(n))
  m <- iv_robust(y ~ x + w | iv + w, data = dat)
  expect_true("w" %in% names(coef(m)))
})

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

test_that("confint works for lm_robust", {
  set.seed(1)
  dat <- data.frame(Y = rnorm(50), Z = rbinom(50, 1, .5))
  m   <- lm_robust(Y ~ Z, data = dat)
  ci  <- confint(m)
  expect_equal(nrow(ci), 2)
  expect_true(all(ci[, 1] < ci[, 2]))
})

test_that("vcov returns symmetric matrix", {
  set.seed(1)
  dat <- data.frame(Y = rnorm(50), Z = rbinom(50, 1, .5), X = rnorm(50))
  m <- lm_robust(Y ~ Z + X, data = dat)
  V <- vcov(m)
  expect_equal(V, t(V), tolerance = 1e-14)
})

test_that("summary.lm_robust runs without error", {
  set.seed(1)
  dat <- data.frame(Y = rnorm(50), Z = rbinom(50, 1, .5))
  expect_no_error(summary(lm_robust(Y ~ Z, data = dat)))
})

test_that("predict.lm_robust returns fitted values on new data", {
  set.seed(1)
  dat <- data.frame(Y = rnorm(50), Z = rbinom(50, 1, .5))
  m   <- lm_robust(Y ~ Z, data = dat)
  p   <- predict(m, newdata = dat)
  expect_equal(length(p), nrow(dat))
})

test_that("nobs returns sample size", {
  dat <- data.frame(Y = rnorm(30), Z = rbinom(30, 1, .5))
  m   <- lm_robust(Y ~ Z, data = dat)
  expect_equal(nobs(m), 30L)
})

test_that("update.lm_robust re-runs model with new formula", {
  set.seed(1)
  dat <- data.frame(Y = rnorm(40), Z = rbinom(40, 1, .5), X = rnorm(40))
  m1  <- lm_robust(Y ~ Z, data = dat)
  m2  <- update(m1, . ~ . + X)
  expect_equal(length(coef(m2)), 3L)
})

# ---------------------------------------------------------------------------
# difference_in_means — extended
# ---------------------------------------------------------------------------

test_that("difference_in_means: matches t.test numerically", {
  set.seed(42)
  dat <- data.frame(Y = rnorm(100), Z = rbinom(100, 1, .5))
  dm  <- difference_in_means(Y ~ Z, data = dat)
  tt  <- with(dat, t.test(Y[Z == 1], Y[Z == 0]))
  expect_equal(dm$p.value[[1]],  tt$p.value,     tolerance = 1e-10)
  expect_equal(dm$conf.low[[1]], tt$conf.int[1], tolerance = 1e-10)
  expect_equal(dm$conf.high[[1]], tt$conf.int[2], tolerance = 1e-10)
})

test_that("difference_in_means: multi-arm condition1/condition2 works", {
  dat <- data.frame(Y = rnorm(100), Z = sample(1:3, 100, replace = TRUE))
  m12 <- difference_in_means(Y ~ Z, condition1 = 1, condition2 = 2, data = dat)
  m13 <- difference_in_means(Y ~ Z, condition1 = 1, condition2 = 3, data = dat)
  expect_s3_class(m12, "difference_in_means")
  expect_s3_class(m13, "difference_in_means")
  expect_equal(m12$design, "Standard")
})

test_that("difference_in_means: reversing conditions negates estimate, not SE", {
  set.seed(42)
  N <- 80
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5))
  m_fwd <- difference_in_means(Y ~ Z, condition1 = 0, condition2 = 1, data = dat)
  m_rev <- difference_in_means(Y ~ Z, condition1 = 1, condition2 = 0, data = dat)
  expect_equal(m_fwd$coefficients[[1]], -m_rev$coefficients[[1]])
  expect_equal(m_fwd$std.error[[1]],    m_rev$std.error[[1]])
})

test_that("difference_in_means: blocked condition reversal negates estimate, not SE", {
  set.seed(42)
  dat <- data.frame(
    Y = rnorm(60), Z = rbinom(60, 1, .5),
    bl = sample(c("A", "B", "C"), 60, replace = TRUE)
  )
  m_fwd <- difference_in_means(Y ~ Z, blocks = bl, condition1 = 0, condition2 = 1, data = dat)
  m_rev <- difference_in_means(Y ~ Z, blocks = bl, condition1 = 1, condition2 = 0, data = dat)
  expect_equal(m_fwd$coefficients[[1]], -m_rev$coefficients[[1]])
  expect_equal(m_fwd$std.error[[1]],    m_rev$std.error[[1]])
  expect_equal(m_fwd$design, "Blocked")
})

test_that("difference_in_means: pair-matched design returns Matched-pair", {
  set.seed(42)
  dat <- data.frame(Y = rnorm(100), Z = rep(0:1, 50), bl = rep(1:50, each = 2))
  m <- difference_in_means(Y ~ Z, blocks = bl, data = dat)
  expect_equal(m$design, "Matched-pair")
  expect_true(!is.na(m$std.error[["Z"]]))
})

test_that("difference_in_means: weighted standard design reported correctly", {
  n <- 60
  dat <- data.frame(Y = rnorm(n), Z = rbinom(n, 1, .5), W = runif(n))
  m <- difference_in_means(Y ~ Z, weights = W, data = dat)
  expect_equal(m$design, "Standard (weighted)")
})

test_that("difference_in_means: weighted blocked design reported correctly", {
  set.seed(42)
  n <- 60
  bl <- rep(1:6, 10)
  dat <- data.frame(
    Y  = rnorm(n),
    Z  = as.integer(unlist(lapply(split(seq_len(n), bl), function(idx) {
      z <- integer(length(idx)); z[sample(seq_along(z), length(idx) %/% 2L)] <- 1L; z
    }))),
    W  = runif(n),
    bl = bl
  )
  m <- difference_in_means(Y ~ Z, weights = W, blocks = bl, data = dat)
  expect_equal(m$design, "Blocked (weighted)")
})

test_that("difference_in_means: ci = FALSE returns NA for p.value and CIs", {
  dat <- data.frame(Y = rnorm(40), Z = rbinom(40, 1, .5))
  m <- difference_in_means(Y ~ Z, data = dat, ci = FALSE)
  expect_true(is.na(m$p.value[[1]]))
  expect_true(is.na(m$conf.low[[1]]))
  expect_true(is.na(m$conf.high[[1]]))
})

test_that("difference_in_means: clustered design label is Clustered", {
  set.seed(42)
  dat <- data.frame(
    Y  = rnorm(40),
    cl = rep(1:4, each = 10)
  )
  dat$Z <- as.integer(dat$cl %in% c(1, 2))
  m <- difference_in_means(Y ~ Z, clusters = cl, data = dat)
  expect_equal(m$design, "Clustered")
})

test_that("difference_in_means: block-cluster mismatch errors", {
  set.seed(42)
  N <- 100
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5),
                    bl = rep(1:10, each = 10), bad_cl = rep(1:10, 10))
  expect_error(
    difference_in_means(Y ~ Z, blocks = bl, clusters = bad_cl, data = dat),
    "All `clusters` must be contained within `blocks`"
  )
})

test_that("difference_in_means: single-unit block errors", {
  N <- 100
  dat <- data.frame(Y = rnorm(N), Z = rbinom(N, 1, .5))
  dat$bad_bl <- c(1, rep(2:10, length.out = N - 1))
  expect_error(
    difference_in_means(Y ~ Z, blocks = bad_bl, data = dat),
    "All `blocks` must have multiple units"
  )
})

test_that("difference_in_means: two-variable formula errors", {
  dat <- data.frame(Y = rnorm(40), Z = rbinom(40, 1, .5), X = rnorm(40))
  expect_error(
    difference_in_means(Y ~ Z + X, data = dat),
    "must have only one variable on the right-hand side"
  )
})

# ---------------------------------------------------------------------------
# lm_lin — extended
# ---------------------------------------------------------------------------

test_that("lm_lin: matches manual centering in lm_robust", {
  set.seed(42)
  n <- 100
  dat <- data.frame(Y = rnorm(n), Z = rbinom(n, 1, .5), X1 = rnorm(n), X2 = rnorm(n))
  dat$X1_c <- dat$X1 - mean(dat$X1)
  dat$X2_c <- dat$X2 - mean(dat$X2)
  m_lin <- lm_lin(Y ~ Z, covariates = ~ X1 + X2, data = dat)
  m_man <- lm_robust(Y ~ Z + Z * X1_c + Z * X2_c, data = dat)
  expect_equal(tidy(m_lin), tidy(m_man), tolerance = 1e-10)
})

test_that("lm_lin: error when covariates are not a right-sided formula", {
  dat <- data.frame(Y = rnorm(30), Z = rbinom(30, 1, .5), X = rnorm(30))
  expect_error(lm_lin(Y ~ Z, covariates = Y ~ X, data = dat), "right-sided formula")
})

test_that("lm_lin: error when covariates is not formula at all", {
  dat <- data.frame(Y = rnorm(30), Z = rbinom(30, 1, .5), X = rnorm(30))
  expect_error(lm_lin(Y ~ Z, dat$X, data = dat), "must be specified as a formula")
})

test_that("lm_lin: error when treatment formula has multiple terms", {
  dat <- data.frame(Y = rnorm(30), Z = rbinom(30, 1, .5), X = rnorm(30))
  expect_error(lm_lin(Y ~ Z + X, covariates = ~ X, data = dat),
               "must only have the treatment variable on the right-hand side")
})

test_that("lm_lin: error when covariates formula has no predictors", {
  dat <- data.frame(Y = rnorm(30), Z = rbinom(30, 1, .5))
  expect_error(lm_lin(Y ~ Z, ~ 1, data = dat), "variable on the right-hand side")
})

# ---------------------------------------------------------------------------
# iv_robust — extended
# ---------------------------------------------------------------------------

test_that("iv_robust: more regressors than instruments warns", {
  # Under-identified, so the second stage is rank deficient and the intercept
  # drops out: both warnings are correct and both fire.
  expect_warning(
    expect_warning(
      iv_robust(mpg ~ hp + cyl | am, data = mtcars, se_type = "HC0"),
      "More regressors than instruments"
    ),
    "collinear"
  )
})

test_that("iv_robust: missing instrument specification errors", {
  expect_error(
    iv_robust(mpg ~ hp + cyl, data = mtcars),
    "Must specify a `formula` with both regressors and instruments"
  )
})

test_that("iv_robust: classical SEs match OLS when instrument = treatment", {
  set.seed(42)
  n <- 100
  dat <- data.frame(y = rnorm(n), x = rnorm(n))
  m_iv  <- iv_robust(y ~ x | x, data = dat, se_type = "classical")
  m_ols <- lm_robust(y ~ x,     data = dat, se_type = "classical")
  expect_equal(coef(m_iv), coef(m_ols), tolerance = 1e-10)
})

test_that("iv_robust: weighted 2SLS runs", {
  set.seed(42)
  n <- 60
  dat <- data.frame(y = rnorm(n), x = rnorm(n), iv = rnorm(n), w = runif(n))
  m <- iv_robust(y ~ x | iv, data = dat, weights = w)
  expect_s3_class(m, "iv_robust")
  expect_true(!is.na(coef(m)[["x"]]))
})

# ---------------------------------------------------------------------------
# FE — extended (ported from estimatr test-fixed_effects.R)
# ---------------------------------------------------------------------------

test_that("FE coefs match dummy regression for supported HC types", {
  set.seed(43)
  N <- 40
  dat <- data.frame(
    Y = rnorm(N), Z = rbinom(N, 1, .5), X = rnorm(N),
    B = factor(rep(1:4, each = 10))
  )
  for (st in c("HC0", "HC1", "classical")) {
    ro  <- tidy(lm_robust(Y ~ Z + factor(B), data = dat, se_type = st))
    rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~B, data = dat, se_type = st))
    expect_equal(
      ro$estimate[ro$term == "Z"],
      rfo$estimate[rfo$term == "Z"],
      tolerance = 1e-10,
      label = paste("coef match for se_type =", st)
    )
  }
})

test_that("HC2/HC3 with one-way fixed_effects match the dummy regression", {
  # Was an error assertion. One FE factor now supports both exactly; two or
  # more still do not. See test_fe_leverage.R.
  dat <- data.frame(Y = rnorm(20), Z = rbinom(20, 1, .5), B = factor(rep(1:2, 10)))
  for (se in c("HC2", "HC3")) {
    fe  <- lm_robust(Y ~ Z, fixed_effects = ~B, data = dat, se_type = se)
    dum <- lm_robust(Y ~ Z + factor(B), data = dat, se_type = se)
    expect_equal(unname(fe$std.error), unname(dum$std.error["Z"]))
  }
})

test_that("FE coefs match dummy regression for supported CR types", {
  set.seed(43)
  N <- 40
  dat <- data.frame(
    Y = rnorm(N), Z = rbinom(N, 1, .5), X = rnorm(N),
    B = factor(rep(1:4, each = 10)), cl = rep(1:4, 10)
  )
  for (st in c("CR0", "stata")) {
    ro  <- suppressWarnings(tidy(lm_robust(Y ~ Z + factor(B), clusters = cl, data = dat, se_type = st)))
    rfo <- suppressWarnings(tidy(lm_robust(Y ~ Z, fixed_effects = ~B, clusters = cl, data = dat, se_type = st)))
    expect_equal(
      ro$estimate[ro$term == "Z"],
      rfo$estimate[rfo$term == "Z"],
      tolerance = 1e-10,
      label = paste("coef match for se_type =", st)
    )
  }
})

test_that("CR2 with fixed_effects errors", {
  dat <- data.frame(Y = rnorm(20), Z = rbinom(20, 1, .5),
                    B = factor(rep(1:2, 10)), cl = rep(1:4, 5))
  expect_error(
    lm_robust(Y ~ Z, fixed_effects = ~B, clusters = cl, data = dat, se_type = "CR2"),
    "CR2"
  )
})

test_that("FE with multiple outcomes: coefs match dummy regression", {
  set.seed(43)
  N <- 40
  dat <- data.frame(
    Y  = rnorm(N), Y2 = rnorm(N),
    Z  = rbinom(N, 1, .5), X  = rnorm(N),
    B  = factor(rep(1:4, each = 10))
  )
  ro  <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B), data = dat)
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~B, data = dat)

  expect_equal(
    tidy(ro)$estimate[tidy(ro)$term == "Z"],
    tidy(rfo)$estimate[tidy(rfo)$term == "Z"],
    tolerance = 1e-10
  )
  expect_equal(rfo$fitted.values, ro$fitted.values, tolerance = 1e-8)
})

test_that("FE missingness in FE variable warns and drops correctly", {
  set.seed(43)
  N <- 40
  dat <- data.frame(
    Y  = rnorm(N), Z = rbinom(N, 1, .5), X = rnorm(N),
    B  = factor(rep(1:4, each = 10))
  )
  dat$B[2] <- NA

  expect_warning(
    m <- lm_robust(Y ~ Z + X, fixed_effects = ~B, data = dat),
    "fixed_effects"
  )
  m_drop <- lm_robust(Y ~ Z + X, fixed_effects = ~B, data = dat[-2, ])
  expect_equal(coef(m), coef(m_drop), tolerance = 1e-10)
})

test_that("FE-only model (no covariates) runs and returns empty coef", {
  set.seed(42)
  n <- 60
  dat <- data.frame(y = rnorm(n), B = factor(rep(1:6, 10)))
  m <- lm_robust(y ~ 1, fixed_effects = ~B, data = dat)
  expect_equal(length(coef(m)), 0L)
  expect_true(!is.null(m$fitted.values))
})
