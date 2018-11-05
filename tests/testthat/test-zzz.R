context("zzz.R - .onLoad")

test_that("onLoad makes generics if texreg is present", {

  e <- environment(.onLoad)

  environment(.onLoad) <- new.env(parent = parent.env(e))
  environment(.onLoad)$extract.lm_robust <- e$extract.lm_robust

  expect_null(.onLoad("estimatr", "estimatr"))
  environment(.onLoad) <- e
})
