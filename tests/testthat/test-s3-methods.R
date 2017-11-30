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

test_that('predict works', {
  set.seed(42)
  n <- 10
  dat <- data.frame(x = rep(0:1, times = 5),
                    w = runif(n),
                    z = rnorm(n),
                    cl = as.factor(sample(letters[1:2], size = n, replace = T)),
                    y = rnorm(n))

  lm_out <- lm(y ~ z * x + cl, data = dat)
  lmr_out <- lm_robust(y ~ z * x + cl, data = dat, se_type = 'classical')

  expect_equivalent(
    predict(lm_out, dat),
    predict(lmr_out, dat)
  )

  # various specifications
  expect_equivalent(
    predict(lm_out, dat, se.fit = T, interval = 'confidence')[c(1, 2)],
    predict(lmr_out, dat, se.fit = T, interval = 'confidence')[c(1, 2)]
  )

  expect_equivalent(
    predict(lm_out, dat, se.fit = T, interval = 'confidence')[c(1, 2)],
    predict(lmr_out, dat, se.fit = T, interval = 'confidence')[c(1, 2)]
  )

  # missingness
  n <- 11
  new_dat <- data.frame(x = rep(0:1, times = c(5, 6)),
                    w = runif(n),
                    z = rnorm(n),
                    cl = as.factor(sample(letters[1:2], size = n, replace = T)),
                    y = rnorm(n))
  # remove one level
  new_dat <- new_dat[new_dat$cl == 'a', ]
  new_dat[1, 'x'] <- NA

  expect_equivalent(
    predict(lm_out, new_dat),
    predict(lmr_out, new_dat)
  )

  # various specifications
  expect_equivalent(
    predict(lm_out, new_dat, se.fit = T, interval = 'confidence')[c(1, 2)],
    predict(lmr_out, new_dat, se.fit = T, interval = 'confidence')[c(1, 2)]
  )

  expect_equivalent(
    predict(lm_out, new_dat, se.fit = T, interval = 'confidence')[c(1, 2)],
    predict(lmr_out, new_dat, se.fit = T, interval = 'confidence')[c(1, 2)]
  )

  # weights
  lm_out <- lm(y ~ z * x + cl, data = dat, weights = w)
  lmr_out <- lm_robust(y ~ z * x + cl, data = dat, se_type = 'classical', weights = w)

  expect_equivalent(
    predict(lm_out, dat),
    predict(lmr_out, dat)
  )

  expect_equivalent(
    predict(lm_out, dat, se.fit = T, interval = 'confidence')[c(1, 2)],
    predict(lmr_out, dat, se.fit = T, interval = 'confidence')[c(1, 2)]
  )

  expect_warning(
    plmo <- predict(lm_out, dat, se.fit = T, interval = 'prediction')[c(1, 2)],
    'Assuming constant prediction variance'
  )
  expect_warning(
    plmro <- predict(lmr_out, dat, se.fit = T, interval = 'prediction')[c(1, 2)],
    'Assuming constant prediction variance'
  )
  expect_equivalent(
    plmo,
    plmro
  )

  # Now with missingness and newdat

  expect_equivalent(
    predict(lm_out, new_dat),
    predict(lmr_out, new_dat)
  )


  # mimic lm behavior with missing weights (can't get prediction intervals)
  new_dat$w[3] <- NA
  # various specifications
  expect_equivalent(
    predict(lm_out, new_dat, se.fit = T, interval = 'confidence', weights = ~ w)[c(1, 2)],
    predict(lmr_out, new_dat, se.fit = T, interval = 'confidence', weights = w)[c(1, 2)]
  )

  expect_equivalent(
    predict(lm_out, new_dat, se.fit = T, interval = 'prediction', weights = ~w)[c(1, 2)],
    predict(lmr_out, new_dat, se.fit = T, interval = 'prediction', weights = w)[c(1, 2)]
  )

  # other arguments
  expect_equivalent(
    predict(lm_out, new_dat, se.fit = T, interval = 'prediction', pred.var = 2.3)[c(1, 2)],
    predict(lmr_out, new_dat, se.fit = T, interval = 'prediction', pred.var = 2.3)[c(1, 2)]
  )

  # lm_lin

  # TODO works with rank deficient X
})

