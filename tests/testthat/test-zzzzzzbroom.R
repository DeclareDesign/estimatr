# This file is ignored by .Rbuildignore to keep from suggesting broom

context("S3 - tidy and broom compatability")

test_that("estimatr::tidy works loaded before broom", {

  library(broom)

  model <- lm(extra~group, sleep)
  model2 <- lm_robust(extra~group, sleep, clusters = ID)

  expect_identical(
    environmentName(environment(get("tidy", envir = globalenv(), inherit = TRUE))),
    "broom"
  )

  expect_error(
    broom_tidy_lm <- tidy(model),
    NA
  )

  expect_error(
    broom_tidy_lm_robust <- tidy(model2),
    NA
  )

})

test_that("broom::tidy works if estimatr loaded after", {

  skip("Skip as detach messes up coveralls")
  detach("package:estimatr", unload = TRUE)

  library(broom)
  library(estimatr)

  model <- lm(extra~group, sleep)
  model2 <- lm_robust(extra~group, sleep, clusters = ID)

  expect_identical(
    environmentName(environment(get("tidy", envir = globalenv(), inherit = TRUE))),
    "estimatr"
  )

  expect_error(
    estimatr_tidy_lm <- tidy(model),
    NA
  )

  expect_error(
    estimatr_tidy_lm_robust <- tidy(model2),
    NA
  )

})

# Can't test unless both can be run...
# expect_identical(
#   broom_tidy_lm,
#   estimatr_tidy_lm
# )
#
# expect_identical(
#   broom_tidy_lm_robust,
#   estimatr_tidy_lm_robust
# )
