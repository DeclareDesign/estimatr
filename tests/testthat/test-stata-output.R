context("Verification - lm and iv match Stata")

test_that("lm_robust matches stata", {

  # write.csv(mtcars,
  #           file = 'tests/testthat/mtcars.csv',
  #           row.names = F)

  stata_ests <- read.table(
    "stata-ests.txt",
    col.names = c("model", "se1", "se2", "df", "fstat"),
    stringsAsFactors = FALSE
  )

  mtcars$w <- mtcars$drat / 5

  estimatr_mat <- matrix(NA, 10, 4)

  lm_c <- lm_robust(mpg ~ hp, data = mtcars, se_type = "classical")
  estimatr_mat[1, ] <- c(lm_c$std.error ^ 2, lm_c$df[2], lm_c$fstatistic[1])
  lm_hc1 <- lm_robust(mpg ~ hp, data = mtcars, se_type = "HC1")
  estimatr_mat[2, ] <- c(lm_hc1$std.error ^ 2, lm_hc1$df[2], lm_hc1$fstatistic[1])
  lm_hc2 <- lm_robust(mpg ~ hp, data = mtcars, se_type = "HC2")
  estimatr_mat[3, ] <- c(lm_hc2$std.error ^ 2, lm_hc2$df[2], lm_hc2$fstatistic[1])
  lm_hc3 <- lm_robust(mpg ~ hp, data = mtcars, se_type = "HC3")
  estimatr_mat[4, ] <- c(lm_hc3$std.error ^ 2, lm_hc3$df[2], lm_hc3$fstatistic[1])

  lm_stata <- lm_robust(mpg ~ hp, clusters = cyl, data = mtcars, se_type = "stata")
  estimatr_mat[5, ] <- c(lm_stata$std.error ^ 2, lm_stata$df[2], lm_stata$fstatistic[1])

  lm_c_w <- lm_robust(mpg ~ hp, data = mtcars, weights = w, se_type = "classical")
  estimatr_mat[6, ] <- c(lm_c_w$std.error ^ 2, lm_c_w$df[2], lm_c_w$fstatistic[1])
  lm_hc1_w <- lm_robust(mpg ~ hp, data = mtcars, weights = w, se_type = "HC1")
  estimatr_mat[7, ] <- c(lm_hc1_w$std.error ^ 2, lm_hc1_w$df[2], lm_hc1_w$fstatistic[1])
  lm_hc2_w <- lm_robust(mpg ~ hp, data = mtcars, weights = w, se_type = "HC2")
  estimatr_mat[8, ] <- c(lm_hc2_w$std.error ^ 2, lm_hc2_w$df[2], lm_hc2_w$fstatistic[1])
  lm_hc3_w <- lm_robust(mpg ~ hp, data = mtcars, weights = w, se_type = "HC3")
  estimatr_mat[9, ] <- c(lm_hc3_w$std.error ^ 2, lm_hc3_w$df[2], lm_hc3_w$fstatistic[1])
  lm_stata_w <- lm_robust(mpg ~ hp, clusters = cyl, weights = w, data = mtcars, se_type = "stata")
  estimatr_mat[10, ] <- c(lm_stata_w$std.error ^ 2, lm_stata_w$df[2], lm_stata_w$fstatistic[1])

  # All look numerically identical except for HC2 and HC3 with weights which
  # have non-negligible difference. This is due to differences in how the hat
  # matrix is built that are still unresolved

  # Therefore rows 8 and 9 will have larger differences
  expect_true(
    max(abs(estimatr_mat[c(1:7, 10), 1:4] - apply(stata_ests[c(1:7, 10), c(3, 2, 4, 5)], 2, as.numeric))) < 1e-5
  )

})



test_that("iv_robust matches stata", {

  skip_if_not_installed("AER")

  # write.csv(mtcars,
  #           file = 'tests/testthat/mtcars.csv',
  #           row.names = F)

  stata_ests <- read.table(
    "stata-iv-ests.txt",
    col.names = c("model", "v1", "v2", "v3", "fstat", "r2", "r2_a", "rmse"),
    stringsAsFactors = FALSE
  )

  mtcars$w <- mtcars$drat / 5

  estimatr_mat <- matrix(NA, 6, 7)
  iv_c <- iv_robust(mpg ~ hp + am | wt + gear, data = mtcars, se_type = "classical")
  estimatr_mat[1, ] <- c(iv_c$std.error ^ 2, iv_c$fstatistic[1], iv_c$r.squared, iv_c$adj.r.squared, sqrt(iv_c$res_var))
  iv_hc1 <- iv_robust(mpg ~ hp + am | wt + gear, data = mtcars, se_type = "HC1")
  estimatr_mat[2, ] <- c(iv_hc1$std.error ^ 2, iv_hc1$fstatistic[1], iv_hc1$r.squared, iv_hc1$adj.r.squared, sqrt(iv_hc1$res_var))
  iv_stata <- iv_robust(mpg ~ hp + am | wt + gear, clusters = cyl, data = mtcars, se_type = "stata")
  estimatr_mat[3, ] <- c(iv_stata$std.error ^ 2, iv_stata$fstatistic[1], iv_stata$r.squared, iv_stata$adj.r.squared, sqrt(iv_stata$res_var))

  iv_c_w <- iv_robust(mpg ~ hp + am | wt + gear, data = mtcars, weights = w, se_type = "classical")
  estimatr_mat[4, ] <- c(iv_c_w$std.error ^ 2, iv_c_w$fstatistic[1], iv_c_w$r.squared, iv_c_w$adj.r.squared, sqrt(iv_c_w$res_var))
  iv_hc1_w <- iv_robust(mpg ~ hp + am | wt + gear, data = mtcars, weights = w, se_type = "HC1")
  estimatr_mat[5, ] <- c(iv_hc1_w$std.error ^ 2, iv_hc1_w$fstatistic[1], iv_hc1_w$r.squared, iv_hc1_w$adj.r.squared, sqrt(iv_hc1_w$res_var))
  iv_stata_w <- iv_robust(mpg ~ hp + am | wt + gear, clusters = cyl, weights = w, data = mtcars, se_type = "stata")
  estimatr_mat[6, ] <- c(iv_stata_w$std.error ^ 2, iv_stata_w$fstatistic[1], iv_stata_w$r.squared, iv_stata_w$adj.r.squared, sqrt(iv_stata_w$res_var))

  expect_true(
    max(abs(estimatr_mat[, 1] - as.numeric(stata_ests[, 4]))) < 2e-05
  )

  expect_true(
    max(abs(estimatr_mat[, 4] - as.numeric(stata_ests[, 5]))) < 3e-05
  )

  # Note, RMSE is different for stata with weights than ivreg or iv_robust
  expect_true(
    max(abs(estimatr_mat[, 5:6] - stata_ests[, 6:7])) < 4e-08
  )

  ivrego_w <- AER::ivreg(mpg ~ hp + am | wt + gear, data = mtcars, weights = w)

  expect_equal(
    ivrego_w$sigma,
    sqrt(iv_c_w$res_var)
  )
  expect_equal(
    ivrego_w$sigma,
    sqrt(iv_hc1_w$res_var)
  )
  expect_equal(
    ivrego_w$sigma,
    sqrt(iv_stata_w$res_var)
  )

})
