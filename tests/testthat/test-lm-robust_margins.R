context("Helper - lm_robust margins")

skip_if_not_installed("margins")

# install.packages(c("margins", "prediction"))
# devtools::install_github("leeper/prediction")
# devtools::install_github("leeper/margins")
library(margins)

mv <- c("AME", "SE", "z", "p")

test_that("lm robust can work with margins",{

  x <- lm(mpg ~ cyl * hp + wt, data = mtcars)
  lmr <- lm_robust(mpg ~ cyl * hp + wt, data = mtcars)

  # Note old package vce defaults to delta
  # new margins on github defaults to none with our obj
  lm_sum_marg <- summary(
    margins::margins(
      x,
      vcov = sandwich::vcovHC(x, type = "HC2"),
      vce = "delta"
    )
  )

  lmr_sum_marg <- summary(margins::margins(lmr, vce = "delta"))

  # Close enough with HC2?
  expect_true(
    max(lm_sum_marg[, mv] - lmr_sum_marg[, mv]) < 2e-05
  )

  # Close with classical
  lmr_class <- lm_robust(mpg ~ cyl * hp + wt, data = mtcars, se_type = "classical")
  lmrc <- summary(margins(lmr_class, vce = "delta"))
  lmc <- summary(margins(x, vce = "delta"))
  expect_true(
    max(lmc[, mv] - lmrc[, mv]) < 2e-05
  )

  # Works with other vce
  set.seed(42)
  lmrc <- summary(margins(lmr_class, vce = "bootstrap", iterations = 10L))
  expect_true(!any(is.na(lmrc)))
  lmrc <- summary(margins(lmr_class, vce = "simulation", iterations = 10L))
  expect_true(!any(is.na(lmrc)))
  lmrc <- summary(margins(lmr_class, vce = "simulation", iterations = 10L))
  expect_true(!any(is.na(lmrc)))

})

test_that("lm robust + weights can work with margins",{

  x <- lm(mpg ~ cyl * hp, data = mtcars, weights = wt)
  x2 <- lm_robust(mpg ~ cyl * hp, data = mtcars, weights = wt, se_type = "classical")
  expect_equal(marginal_effects(x), marginal_effects(x2))


  expect_warning(
    lmc <- round(summary(margins(x, vce = "delta"))[, mv], 3),
    "'weights'"
  )
  expect_warning(
    lmr <- round(summary(margins(x2, vce = "delta"))[, mv], 3),
    "'weights'"
  )
  # Have to round quite a bit!
  expect_equal(lmc, lmr)

})

test_that("lm robust + cluster can work with margins",{

  # works but throws a lot of warnings
  x <- lm(mpg ~ cyl * hp + wt, data = mtcars)
  x2 <- lm_robust(mpg ~ cyl * hp + wt, data = mtcars, clusters = am)

  lmc <- round(summary(margins(x, vce = "delta"))[, mv], 8)

  expect_warning(
    lmr <- round(summary(margins(x2, vce = "delta"))[, mv], 8),
    "cluster"
  )

  # With rounding
  expect_equal(lmc[, 1], lmr[, 1])
  expect_true(
    !any(lmc[, 2] == lmr[, 2])
  )

})


test_that("lm lin can work with margins",{

  data("alo_star_men")
  lml <- lm_lin(GPA_year1 ~ ssp, ~  gpa0, data = alo_star_men, se_type = "classical")

  alo_star_men$gpa0_tilde <- alo_star_men$gpa0 - mean(alo_star_men$gpa0)

  lmo <- lm(GPA_year1 ~ ssp * gpa0_tilde, data = alo_star_men)

  lml_sum <- summary(margins(lml, vce = "delta"))
  lmo_sum <- summary(margins(lmo, vce = "delta"))

  expect_equal(
    round(lml_sum[, 4], 5),
    round(lmo_sum[, 4], 5)
  )


})

