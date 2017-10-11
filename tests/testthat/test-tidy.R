context('Test tidy methods')

test_that('tidy works', {
  n <- 10
  dat <- data.frame(x = rbinom(n, size = 1, prob = 0.5),
                    p = 0.5,
                    z = rnorm(n),
                    y = rnorm(n))

  expect_is(
    tidy(lm_robust(y ~ x, data = dat, se_type = 'classical')),
    'data.frame'
  )

  expect_equal(
    nrow(tidy(lm_robust(y ~ x, coefficient_name = 'x', data = dat, se_type = 'classical'))),
    1
  )

  expect_is(
    tidy(lm_lin(y ~ x, ~ z, data = dat)),
    'data.frame'
  )

  expect_is(
    tidy(horvitz_thompson(y ~ x, condition_pr_variable_name = p, data = dat)),
    "data.frame"
  )

  expect_is(
    tidy(difference_in_means(y ~ x, data = dat)),
    "data.frame"
  )

})
