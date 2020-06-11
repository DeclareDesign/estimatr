context("Estimator - Arg checking fails as expected.")

test_that("#349 Early fail when formula is a string", {
  expect_error(
    estimatr::lm_robust("mpg~hp", data = mtcars, cluster = wt),
    "formula"
  )
  expect_length(
    estimatr::lm_robust(mpg~hp, data = mtcars, cluster = wt),
    29
  )
})
