context("statprep")

test_that("stargazer works", {

  Y <- rnorm(100)
  X <- rbinom(100, 1, .5)
  fit_1 <- lm(Y ~ X)
  fit_2 <- lm(Y ~ X)

  fit_1_r <- lm_robust(Y ~ X)
  fit_2_r <- lm_robust(Y ~ X)

  library(stargazer)
  stargazer(fit_1, fit_2,
                       se = starprep(fit_1_r, fit_2_r),
                       p = starprep(fit_1_r, fit_2_r, stat = "p"))

})
