context("zzz.R - .onLoad")

test_that("onLoad makes generics if texreg is present", {

  # .onLoad() calls setGeneric(), which writes a binding into the package
  # namespace. That namespace is locked once the package is loaded, and it is
  # always locked under covr, so re-running .onLoad() here aborts with
  # "cannot add bindings to a locked environment". The test only ever passed
  # where the namespace happened to be unlocked; skip where it cannot run
  # rather than fail.
  skip_if(environmentIsLocked(asNamespace("estimatr")))

  e <- environment(.onLoad)

  environment(.onLoad) <- new.env(parent = parent.env(e))
  environment(.onLoad)$extract.lm_robust <- e$extract.lm_robust

  expect_null(.onLoad("estimatr", "estimatr"))
  environment(.onLoad) <- e
})
