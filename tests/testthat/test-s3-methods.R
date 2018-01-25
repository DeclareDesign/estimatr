context("S3")

test_that("tidy, summary, and print work", {
  n <- 10
  dat <- data.frame(
    x = rep(0:1, times = 5),
    p = 0.5,
    z = rnorm(n),
    y = rnorm(n)
  )

  expect_is(
    tidy(NULL),
    "data.frame"
  )

  # Might need to force this test of tidy.default

  d <- date()
  expect_warning(
    estimatr:::tidy.default(d),
    "using as.data.frame"
  )

  ## lm_robust
  lmo <- lm_robust(y ~ x, data = dat, se_type = "classical")

  capture_output(
    summary(lmo)
  )

  expect_is(
    tidy(lmo),
    "data.frame"
  )

  expect_equal(
    nrow(tidy(lm_robust(y ~ x, data = dat, se_type = "classical"))),
    2
  )

  capture_output(
    expect_equivalent(
      summary(lmo)$coefficients,
      print(lmo)
    )
  )


  ## lm_lin
  lmlo <- lm_lin(y ~ x, ~ z, data = dat)
  expect_is(
    tidy(lmlo),
    "data.frame"
  )


  capture_output(
    expect_equivalent(
      summary(lmlo)$coefficients,
      print(lmlo)
    )
  )

  ## horvitz_thompson
  ht <- horvitz_thompson(y ~ x, condition_prs = p, data = dat)
  expect_is(
    tidy(ht),
    "data.frame"
  )

  expect_equivalent(
    as.matrix(tidy(ht)[, c("coefficients", "se", "p", "ci_lower", "ci_upper", "df")]),
    summary(ht)$coefficients
  )


  capture_output(
    expect_equivalent(
      summary(ht)$coefficients,
      print(ht)
    )
  )

  ## difference_in_means
  dim <- difference_in_means(y ~ x, data = dat)
  expect_is(
    tidy(dim),
    "data.frame"
  )

  summary(dim)$coefficients

  capture_output(
    expect_equivalent(
      summary(dim)$coefficients,
      print(dim)
    )
  )

  # rank deficient
  dat$z2 <- dat$z
  lmro <- lm_robust(y ~ z + z2 + x, data = dat)
  tidy(lmro)
})


test_that("vcov works", {
  n <- 10
  dat <- data.frame(
    x = rep(0:1, times = 5),
    p = 0.5,
    z = rnorm(n),
    y = rnorm(n)
  )

  # not identical due to < 1e-15 diffs
  expect_equal(
    vcov(lm_robust(y ~ x, data = dat, se_type = "classical")),
    vcov(lm(y ~ x, data = dat))
  )

  expect_error(
    vcov(lm_lin(y ~ x, ~ z, data = dat)),
    NA
  )

  expect_error(
    vcov(lm_robust(y ~ x, data = dat, return_vcov = FALSE)),
    "return_vcov = TRUE"
  )

  expect_error(
    vcov(horvitz_thompson(y ~ x, condition_prs = p, data = dat)),
    "supported|horvitz_thompson"
  )


  expect_error(
    vcov(difference_in_means(y ~ x, data = dat)),
    "supported|difference_in_means"
  )

  # rank deficient
  dat$z2 <- dat$z
  lmro <- lm_robust(y ~ z + z2 + x, data = dat)
  expect_equivalent(
    dim(vcov(lmro)),
    c(3, 3)
  )
})


test_that("coef and confint work", {
  n <- 10
  dat <- data.frame(
    x = rep(0:1, times = 5),
    p = 0.5,
    z = rnorm(n),
    y = rnorm(n)
  )

  lmo <- lm_robust(y ~ x, data = dat)
  expect_equivalent(
    coef(lmo),
    lmo$coefficients
  )

  expect_equivalent(
    confint(lmo),
    cbind(lmo$ci_lower, lmo$ci_upper)
  )

  lm2o <- lm_robust(y ~ x + z, data = dat)
  expect_equivalent(
    coef(lm2o)[2],
    lm2o$coefficients["x"]
  )

  expect_equivalent(
    confint(lm2o)[2, , drop = FALSE],
    confint(lm2o, parm = "x")
  )

  expect_equivalent(
    confint(lmo, parm = "x", level = 0.90),
    with(
      lm_robust(y ~ x, data = dat, alpha = 0.10),
      cbind(ci_lower[2], ci_upper[2])
    )
  )

  lmlo <- lm_lin(y ~ x, ~ z, data = dat, se_type = "HC3")
  expect_equivalent(
    confint(lmlo),
    cbind(lmlo$ci_lower, lmlo$ci_upper)
  )

  dim <- difference_in_means(y ~ x, data = dat)
  expect_equivalent(
    coef(dim),
    dim$coefficients
  )
  expect_equivalent(
    confint(dim),
    cbind(dim$ci_lower, dim$ci_upper)
  )

  ht <- horvitz_thompson(y ~ x, condition_prs = p, data = dat)
  expect_equivalent(
    coef(ht),
    ht$coefficients
  )
  expect_equivalent(
    confint(ht),
    cbind(ht$ci_lower, ht$ci_upper)
  )

  # rank deficient
  dat$z2 <- dat$z
  lmro <- lm_robust(y ~ z + z2 + x, data = dat)
  confint(lmro)
  coef(lmro)
  capture.output(
    expect_equal(
      nobs(lmro),
      nobs(summary(lmro))
    )
  )

})

test_that("predict works", {
  set.seed(42)
  n <- 10
  dat <- data.frame(
    x = rep(0:1, times = 5),
    w = runif(n),
    z = rnorm(n),
    cl = as.factor(sample(letters[1:3], size = n, replace = T)),
    y = rnorm(n)
  )

  lm_out <- lm(y ~ z * x + cl, data = dat)
  lmr_out <- lm_robust(y ~ z * x + cl, data = dat, se_type = "classical")

  expect_equivalent(
    predict(lm_out, dat),
    predict(lmr_out, dat)
  )

  # various specifications
  expect_equivalent(
    predict(lm_out, dat, se.fit = T, interval = "confidence")[c(1, 2)],
    predict(lmr_out, dat, se.fit = T, interval = "confidence")[c(1, 2)]
  )

  expect_equivalent(
    predict(lm_out, dat, se.fit = T, interval = "prediction")[c(1, 2)],
    predict(lmr_out, dat, se.fit = T, interval = "prediction")[c(1, 2)]
  )

  # missingness
  n <- 11
  new_dat <- data.frame(
    x = rep(0:1, times = c(5, 6)),
    w = runif(n),
    z = rnorm(n),
    cl = as.factor(sample(letters[1:3], size = n, replace = T)),
    y = rnorm(n)
  )
  # remove one level to make sure works with missing levels
  new_dat <- new_dat[new_dat$cl == "a", ]
  new_dat[1, "x"] <- NA

  expect_equivalent(
    predict(lm_out, new_dat),
    predict(lmr_out, new_dat)
  )

  # various specifications
  expect_equivalent(
    predict(lm_out, new_dat, se.fit = T, interval = "confidence")[c(1, 2)],
    predict(lmr_out, new_dat, se.fit = T, interval = "confidence")[c(1, 2)]
  )

  expect_equivalent(
    predict(lm_out, new_dat, se.fit = T, interval = "prediction")[c(1, 2)],
    predict(lmr_out, new_dat, se.fit = T, interval = "prediction")[c(1, 2)]
  )

  # weights
  lm_out <- lm(y ~ z * x + cl, data = dat, weights = w)
  lmr_out <- lm_robust(y ~ z * x + cl, data = dat, se_type = "classical", weights = w)

  expect_equivalent(
    predict(lm_out, dat),
    predict(lmr_out, dat)
  )

  expect_equivalent(
    predict(lm_out, dat, se.fit = T, interval = "confidence")[c(1, 2)],
    predict(lmr_out, dat, se.fit = T, interval = "confidence")[c(1, 2)]
  )

  expect_warning(
    plmo <- predict(lm_out, dat, se.fit = T, interval = "prediction")[c(1, 2)],
    "Assuming constant prediction variance"
  )
  expect_warning(
    plmro <- predict(lmr_out, dat, se.fit = T, interval = "prediction")[c(1, 2)],
    "Assuming constant prediction variance"
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
    predict(lm_out, new_dat, se.fit = T, interval = "confidence", weights = ~ w)[c(1, 2)],
    predict(lmr_out, new_dat, se.fit = T, interval = "confidence", weights = w)[c(1, 2)]
  )

  expect_equivalent(
    predict(lm_out, new_dat, se.fit = T, interval = "prediction", weights = ~w)[c(1, 2)],
    predict(lmr_out, new_dat, se.fit = T, interval = "prediction", weights = w)[c(1, 2)]
  )

  # other arguments
  expect_equivalent(
    predict(lm_out, new_dat, se.fit = T, interval = "prediction", pred.var = 2.3)[c(1, 2)],
    predict(lmr_out, new_dat, se.fit = T, interval = "prediction", pred.var = 2.3)[c(1, 2)]
  )

  # lm_lin
  n <- 11
  new_dat <- data.frame(
    x = rep(0:1, times = c(5, 6)),
    z = rnorm(n),
    cl = as.factor(sample(letters[1:3], size = n, replace = TRUE)),
    y = rnorm(n)
  )

  lml_out <- lm_lin(y ~ x, covariates = ~ z + cl, data = dat, se_type = "classical")

  dat$z_bar <- dat$z - mean(dat$z)
  dat$clb <- as.numeric(dat$cl == "b")
  dat$clc <- as.numeric(dat$cl == "c")
  dat$clb_bar <- dat$clb - mean(dat$clb)
  dat$clc_bar <- dat$clc - mean(dat$clc)

  lm_int_out <- lm(y ~ x + x * z_bar + x * clb_bar + x * clc_bar, data = dat)

  # have to scale new data by old mean values!
  # now predict does this for you! ok.emoji

  new_dat$z_bar <- new_dat$z - mean(dat$z)
  new_dat$clb <- as.numeric(new_dat$cl == "b")
  new_dat$clc <- as.numeric(new_dat$cl == "c")
  new_dat$clb_bar <- new_dat$clb - mean(dat$clb)
  new_dat$clc_bar <- new_dat$clc - mean(dat$clc)

  # not identical due to some numerical difference, presumably due to the way I save the means from lm_lin
  expect_equivalent(
    predict(lml_out, new_dat, se.fit = TRUE, interval = "confidence")[c(1, 2)],
    predict(lm_int_out, new_dat, se.fit = TRUE, interval = "confidence")[c(1, 2)]
  )

  # working with rank deficient X
  head(dat)
  dat$z2 <- dat$z

  lm_out <- lm(y ~ z * x + z2 + cl, data = dat)
  lmr_out <- lm_robust(y ~ z * x + z2 + cl + z, data = dat, se_type = "classical")

  suppressWarnings({
    expect_equivalent(
      predict(lm_out, dat),
      predict(lmr_out, dat)
    )

    # various specifications
    expect_equivalent(
      predict(lm_out, dat, se.fit = T, interval = "confidence")[c(1, 2)],
      predict(lmr_out, dat, se.fit = T, interval = "confidence")[c(1, 2)]
    )

    expect_equivalent(
      predict(lm_out, dat, se.fit = T, interval = "prediction")[c(1, 2)],
      predict(lmr_out, dat, se.fit = T, interval = "prediction")[c(1, 2)]
    )
  })
})
