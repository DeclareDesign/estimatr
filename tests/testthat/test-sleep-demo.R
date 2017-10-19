context("sleep_demo")


test_that("lm_robust works when called from global environment", {
  expect_output_file(
    demo("sleep", package = "estimatr"),
    file = "sleep-demo.txt"
  )
})
