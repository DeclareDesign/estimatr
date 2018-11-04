context("zzz.R - .onLoad")

test_that("onLoad makes generics if texreg is present", {

  e <- environment(.onLoad)

  environment(.onLoad) <- new.env(parent = parent.env(e))
  environment(.onLoad)$extract.lm_robust <- e$extract.lm_robust

  expect_null(.onLoad("estimatr", "estimatr"))
  environment(.onLoad) <- e
})

test_that(".onLoad does not message if new version of 'broom' is loaded", {
  skip_if_not_installed("broom", minimum_version = "0.5.0.9000")
  library(broom)
  expect_silent(.onLoad("estimatr", "estimatr"))

  expect_is(
    tidy(lm_robust(mpg ~ hp, mtcars)),
    "data.frame"
  )
})


test_that(".onLoad message if old version of 'broom' is installed", {
  skip_if_not_installed("broom")
  skip_if(packageVersion("broom") > "0.5.0")
  library(broom)
  expect_message(
    .onLoad("estimatr", "estimatr"),
    "the `broom` package version 0.5.0 or earlier is loaded"
  )
})
