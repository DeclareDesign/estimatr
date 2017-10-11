context('Test vcov methods')

test_that('vcov works', {
  n <- 10
  dat <- data.frame(x = rbinom(n, size = 1, prob = 0.5),
                    p = 0.5,
                    z = rnorm(n),
                    y = rnorm(n))

  # not identical due to < 1e-15 diffs
  expect_equal(
    vcov(lm_robust(y ~ x, data = dat, se_type = 'classical')),
    vcov(lm(y ~ x, data = dat))
  )

  expect_error(
    vcov(lm_lin(y ~ x, ~ z, data = dat)),
    NA
  )

  expect_error(
    vcov(horvitz_thompson(y ~ x, condition_pr_variable_name = p, data = dat)),
    "supported|horvitz_thompson"
  )


  expect_error(
    vcov(difference_in_means(y ~ x, data = dat)),
    "supported|difference_in_means"
  )

})
