library(estimatr)

test_that("horvitz_thompson aligns matrix and data.frame condition_prs after NA drops", {
  set.seed(1)
  n <- 10
  y <- rnorm(n)
  z <- rep(0:1, each = 5)
  p1 <- seq(0.2, 0.8, length.out = n)
  p_mat <- cbind("0" = 1 - p1, "1" = p1)
  p_df <- as.data.frame(p_mat)
  p_tbl <- tibble::as_tibble(p_df)
  dat <- data.frame(y = y, z = z)
  dat_na <- dat
  dat_na$y[1] <- NA

  ref <- horvitz_thompson(y ~ z, data = dat[-1, ], condition_prs = p_mat[-1, ])

  for (cpr in list(p_mat, p_df, p_tbl)) {
    fit <- horvitz_thompson(y ~ z, data = dat_na, condition_prs = cpr)
    expect_equal(fit$coefficients, ref$coefficients, tolerance = 1e-12)
    expect_equal(fit$std.error, ref$std.error, tolerance = 1e-12)
    expect_equal(fit$nobs, n - 1L)
  }
})

test_that("horvitz_thompson aligns matrix and data.frame condition_prs with >2 arms", {
  dat3 <- data.frame(
    y = c(10, 20, 30, 40, 50, 60),
    z = c(2, 0, 1, 2, 0, 1)
  )
  p_mat3 <- cbind(
    "0" = c(0.1, 0.3, 0.4, 0.1, 0.5, 0.2),
    "1" = c(0.8, 0.3, 0.4, 0.8, 0.3, 0.6)
  )

  # Units in condition2 = 1 are rows 3 (y=30, p=0.4) and 6 (y=60, p=0.6).
  # Units in condition1 = 0 are rows 2 (y=20, p=0.3) and 5 (y=50, p=0.5).
  expected_est <- ((30 / 0.4 + 60 / 0.6) - (20 / 0.3 + 50 / 0.5)) / 6
  expected_se <- sqrt(
    ((30 / 0.4)^2 + (60 / 0.6)^2 + (20 / 0.3)^2 + (50 / 0.5)^2) / 6^2
  )

  fit_mat <- horvitz_thompson(
    y ~ z, data = dat3, condition_prs = p_mat3, condition1 = 0, condition2 = 1
  )
  fit_df <- horvitz_thompson(
    y ~ z, data = dat3, condition_prs = as.data.frame(p_mat3),
    condition1 = 0, condition2 = 1
  )

  expect_equal(unname(fit_mat$coefficients), expected_est, tolerance = 1e-12)
  expect_equal(unname(fit_mat$std.error), expected_se, tolerance = 1e-12)
  expect_equal(fit_df$coefficients, fit_mat$coefficients, tolerance = 1e-12)
  expect_equal(fit_df$std.error, fit_mat$std.error, tolerance = 1e-12)
})

test_that("horvitz_thompson aligns per-unit numeric vector condition_prs after NA drops", {
  set.seed(2)
  n <- 10
  y <- rnorm(n)
  z <- rep(0:1, each = 5)
  p1 <- seq(0.25, 0.75, length.out = n)
  dat <- data.frame(y = y, z = z)
  dat_na <- dat
  dat_na$y[c(1, 6)] <- NA

  ref <- horvitz_thompson(y ~ z, data = dat[-c(1, 6), ], condition_prs = p1[-c(1, 6)])
  fit <- horvitz_thompson(y ~ z, data = dat_na, condition_prs = p1)

  expect_equal(fit$coefficients, ref$coefficients, tolerance = 1e-12)
  expect_equal(fit$std.error, ref$std.error, tolerance = 1e-12)
  expect_equal(fit$nobs, n - 2L)
})

test_that("horvitz_thompson refuses matrix/data.frame condition_prs with wrong row count", {
  dat <- data.frame(y = rnorm(10), z = rep(0:1, each = 5))
  bad_mat <- cbind("0" = rep(0.5, 8), "1" = rep(0.5, 8))

  expect_error(
    horvitz_thompson(y ~ z, data = dat, condition_prs = bad_mat),
    "has 8 rows but the data has 10 rows"
  )
  expect_error(
    horvitz_thompson(y ~ z, data = dat, condition_prs = as.data.frame(bad_mat)),
    "has 8 rows but the data has 10 rows"
  )
})
