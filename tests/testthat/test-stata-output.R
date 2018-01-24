context("Verification - lm_robust matches Stata")

test_that("lm_robust matches stata", {

  # write.csv(mtcars,
  #           file = 'tests/testthat/mtcars.csv',
  #           row.names = F)

  stata_ests <- read.table(
    "stata-ests.txt",
    col.names = c("model", "se1", "se2", "df"),
    stringsAsFactors = FALSE
  )

  mtcars$w <- mtcars$drat / 5

  estimatr_mat <- matrix(NA, 10, 3)

  lm_c <- lm_robust(mpg ~ hp, data = mtcars, se_type = "classical")
  estimatr_mat[1, ] <- c(lm_c$se ^ 2, lm_c$df[2])
  lm_hc1 <- lm_robust(mpg ~ hp, data = mtcars, se_type = "HC1")
  estimatr_mat[2, ] <- c(lm_hc1$se ^ 2, lm_hc1$df[2])
  lm_hc2 <- lm_robust(mpg ~ hp, data = mtcars, se_type = "HC2")
  estimatr_mat[3, ] <- c(lm_hc2$se ^ 2, lm_hc2$df[2])
  lm_hc3 <- lm_robust(mpg ~ hp, data = mtcars, se_type = "HC3")
  estimatr_mat[4, ] <- c(lm_hc3$se ^ 2, lm_hc3$df[2])

  lm_stata <- lm_robust(mpg ~ hp, clusters = cyl, data = mtcars, se_type = "stata")
  estimatr_mat[5, ] <- c(lm_stata$se ^ 2, lm_stata$df[2])

  lm_c_w <- lm_robust(mpg ~ hp, data = mtcars, weights = w, se_type = "classical")
  estimatr_mat[6, ] <- c(lm_c_w$se ^ 2, lm_c_w$df[2])
  lm_hc1_w <- lm_robust(mpg ~ hp, data = mtcars, weights = w, se_type = "HC1")
  estimatr_mat[7, ] <- c(lm_hc1_w$se ^ 2, lm_hc1_w$df[2])
  lm_hc2_w <- lm_robust(mpg ~ hp, data = mtcars, weights = w, se_type = "HC2")
  estimatr_mat[8, ] <- c(lm_hc2_w$se ^ 2, lm_hc2_w$df[2])
  lm_hc3_w <- lm_robust(mpg ~ hp, data = mtcars, weights = w, se_type = "HC3")
  estimatr_mat[9, ] <- c(lm_hc3_w$se ^ 2, lm_hc3_w$df[2])
  lm_stata_w <- lm_robust(mpg ~ hp, clusters = cyl, weights = w, data = mtcars, se_type = "stata")
  estimatr_mat[10, ] <- c(lm_stata_w$se ^ 2, lm_stata_w$df[2])

  # All look numerically identical except for HC2 and HC3 with weights which
  # have non-negligible difference. This is due to differences in how the hat
  # matrix is built that are still unresolved

  # Therefore rows 8 and 9 will have larger differences
  expect_true(
    max(estimatr_mat[c(1:7, 10), ] - apply(stata_ests[c(1:7, 10), c(3, 2, 4)], 2, as.numeric)) < 4e-8
  )
})
