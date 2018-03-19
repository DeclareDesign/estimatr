# This file is ignored by .Rbuildignore to keep from suggesting broom

context("S3 - tidy and broom compatability")

test_that("estimatr::tidy works loaded before or after after broom", {
  detach("package:estimatr", unload = TRUE)

  library(broom)
  library(estimatr)

  model <- lm(extra~group, sleep)
  model2 <- lm_robust(extra~group, sleep, clusters = ID)

  expect_identical(
    environmentName(environment(get("tidy", envir = globalenv(), inherit = TRUE))),
    "estimatr"
  )

  estimatr_tidy_lm <- tidy(model)
  estimatr_tidy_lm_robust <- tidy(model2)

  detach("package:estimatr", unload = TRUE)
  detach("package:broom", unload = TRUE)

  library(estimatr)
  library(broom)

  model <- lm(extra~group, sleep)
  model2 <- lm_robust(extra~group, sleep, clusters = ID)

  expect_identical(
    environmentName(environment(get("tidy", envir = globalenv(), inherit = TRUE))),
    "broom"
  )

  broom_tidy_lm <- tidy(model)
  broom_tidy_lm_robust <- tidy(model2)

  expect_identical(
    broom_tidy_lm,
    estimatr_tidy_lm
  )

  expect_identical(
    broom_tidy_lm_robust,
    estimatr_tidy_lm_robust
  )
})
