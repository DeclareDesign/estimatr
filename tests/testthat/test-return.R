context("Output - test similiarity across estimators")

test_that("Structure of output is the same", {

  n <- 40
  dat <- data.frame(
    y = rnorm(n),
    z = c(0, 0, rep(0:1, times = 9)),
    x = rnorm(n),
    bl = rep(1:4, each = 10),
    cl = rep(1:20, each = 2)
  )

  # Should be in all estimator returns
  in_return <-
    c("est",
      "se",
      "df",
      "p",
      "ci_lower",
      "ci_upper",
      "outcome",
      "alpha",
      "N"
    )

  lmr_o <- lm_robust(y ~ z, data = dat)
  lmr_cl_o <- lm_robust(y ~ z, data = dat, clusters = cl)
  lml_o <- lm_lin(y ~ z, ~ x, data = dat)
  lml_cl_o <- lm_lin(y ~ z, ~ x, data = dat)
  # Major branching for diff estimators is for blocks
  ht_o <- horvitz_thompson(y ~ z, data = dat)
  ht_bl_o <- horvitz_thompson(y ~ z, blocks = bl, data = dat)
  dim_o <- difference_in_means(y ~ z, data = dat)
  dim_bl_o <- difference_in_means(y ~ z, blocks = bl, data = dat)

  expect_true( all(in_return %in% names(lmr_o)) )
  expect_true( all(in_return %in% names(lmr_cl_o)) )
  expect_true( all(in_return %in% names(lml_o)) )
  expect_true( all(in_return %in% names(lml_cl_o)) )
  expect_true( all(in_return %in% names(dim_o)) )
  expect_true( all(in_return %in% names(dim_bl_o)) )
  expect_true( all(in_return %in% names(ht_o)) )
  expect_true( all(in_return %in% names(ht_bl_o)) )

  expect_equal(
    colnames(tidy(lmr_o)),
    colnames(tidy(lmr_cl_o)),
    colnames(tidy(lml_o)),
    colnames(tidy(lml_cl_o)),
    colnames(tidy(ht_o)),
    colnames(tidy(ht_bl_o)),
    colnames(tidy(dim_o)),
    colnames(tidy(dim_bl_o))
  )

})

test_that("Warns properly if df is negative or 0", {
  dat = data.frame(y = 1, z = 1, p = .5)
  expect_warning(
    horvitz_thompson(y ~ z, data = dat, condition_prs = p),
    "Estimated negative or zero degrees of freedom"
  )
})
