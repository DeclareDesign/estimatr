context('Test S3 methods')

test_that('tidy, summary, and print work', {
  n <- 10
  dat <- data.frame(x = rep(0:1, times = 5),
                    p = 0.5,
                    z = rnorm(n),
                    y = rnorm(n))

  ## lm_robust
  lmo <- lm_robust(y ~ x, data = dat, se_type = 'classical')

  summary(lmo)$coefficients

  expect_is(
    tidy(lmo),
    'data.frame'
  )

  expect_equal(
    nrow(tidy(lm_robust(y ~ x, coefficient_name = 'x', data = dat, se_type = 'classical'))),
    1
  )

  expect_equivalent(
    tidy(lmo),
    print(lmo)
  )


  ## lm_lin
  lmlo <- lm_lin(y ~ x, ~ z, data = dat)
  expect_is(
    tidy(lmlo),
    'data.frame'
  )

  expect_equivalent(
    tidy(lmlo),
    print(lmlo)
  )

  ## horvitz_thompson
  ht <- horvitz_thompson(y ~ x, condition_prs = p, data = dat)
  expect_is(
    tidy(ht),
    "data.frame"
  )

  expect_equivalent(
    tidy(ht),
    print(ht)
  )

  ## difference_in_means
  dim <- difference_in_means(y ~ x, data = dat)
  expect_is(
    tidy(dim),
    "data.frame"
  )

  summary(dim)$coefficients

  expect_equivalent(
    tidy(dim),
    print(dim)
  )

})


test_that('vcov works', {
  n <- 10
  dat <- data.frame(x = rep(0:1, times = 5),
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
    vcov(horvitz_thompson(y ~ x, condition_prs = p, data = dat)),
    "supported|horvitz_thompson"
  )


  expect_error(
    vcov(difference_in_means(y ~ x, data = dat)),
    "supported|difference_in_means"
  )

})


test_that('coef and confint work', {
  n <- 10
  dat <- data.frame(x = rep(0:1, times = 5),
                    p = 0.5,
                    z = rnorm(n),
                    y = rnorm(n))

  lmo <- lm_robust(y ~ x, data = dat)
  expect_equivalent(
    coef(lmo),
    lmo$est
  )

  expect_equivalent(
    confint(lmo),
    cbind(lmo$ci_lower, lmo$ci_upper)
  )

  expect_equivalent(
    confint(lmo, parm = 'x', level = 0.15),
    with(lm_robust(y ~ x, data = dat, coefficient_name = 'x', alpha = 0.15),
         cbind(ci_lower[2], ci_upper[2]))
  )

  lmlo <- lm_lin(y ~ x, ~ z, data = dat, se_type = 'HC3')
  expect_equivalent(
    confint(lmlo),
    cbind(lmlo$ci_lower, lmlo$ci_upper)
  )

  dim <- difference_in_means(y ~ x, data = dat)
  expect_equivalent(
    confint(dim),
    cbind(dim$ci_lower, dim$ci_upper)
  )

  ht <- horvitz_thompson(y ~ x, condition_prs = p, data = dat)
  expect_equivalent(
    confint(ht),
    cbind(ht$ci_lower, ht$ci_upper)
  )

})

