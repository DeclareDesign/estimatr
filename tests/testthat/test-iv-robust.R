context("Estimator - iv_robust")

N <- 20
dat <- data.frame(
  y = rnorm(N),
  x = rnorm(N),
  x2 = rnorm(N),
  z2 = rnorm(N),
  w = runif(N),
  clust = sample(letters[1:3], size = N, replace = TRUE)
)
dat$z <- dat$x * 0.5 + rnorm(N)

test_that("iv_robust warnings and errors are correct", {
  expect_warning(
    ivro <- iv_robust(mpg ~ hp + cyl | am, data = mtcars, se_type = "HC0"),
    "More regressors than instruments"
  )

  expect_error(
    iv_robust(mpg ~ hp + cyl, data = mtcars),
    "Must specify a `formula` with both regressors and instruments."
  )
})

test_that("iv_robust matches AER + ivpack", {

  ivco <- iv_robust(y ~ x | z, data = dat, se_type = "classical")
  ivfit <- AER::ivreg(y ~ x | z, data = dat)
  ivo <- summary(ivfit)

  expect_equivalent(
    as.matrix(tidy(ivco)[, c("estimate", "std.error", "statistic", "p.value")]),
    coef(ivo)[, 1:4]
  )
  # Same as stata if you specify `small` as a stata option
  # which applies the N / N-k finite sample correction

  expect_equivalent(
    ivfit$fitted.values,
    ivco$fitted.values
  )

  # Stata defaults to HC0 as well, but does HC1 with `small`
  ivro <- iv_robust(y ~ x | z, data = dat, se_type = "HC0")
  capture_output(ivpackrob <- ivpack::robust.se(ivfit))

  expect_equivalent(
    as.matrix(tidy(ivro)[, c("estimate", "std.error", "statistic", "p.value")]),
    ivpackrob[, 1:4]
  )

  expect_equivalent(
    ivfit$fitted.values,
    ivro$fitted.values
  )

  # "Stata" clustered SEs are CR0, but they are the same as below with `small`
  ivclusto <- iv_robust(y ~ x | z, data = dat, se_type = "stata", clusters = clust)
  capture_output(ivpackclust <- ivpack::cluster.robust.se(ivfit, dat$clust))

  # Our p-values are bigger (ivpack is be using less conservative DF, we use J - 1 which
  # is what stata uses for clusters w/ `small` and in OLS)
  expect_equivalent(
    as.matrix(tidy(ivclusto)[, c("estimate", "std.error")]),
    ivpackclust[, c(1, 2)]
  )

  expect_equivalent(
    ivfit$fitted.values,
    ivclusto$fitted.values
  )

  # CR2
  ivcr2o <- iv_robust(y ~ x | z, data = dat, clusters = clust)
  clubsando <- clubSandwich::coef_test(ivfit, vcov = "CR2", cluster = dat$clust)

  expect_equivalent(
    as.matrix(tidy(ivcr2o)[, c("estimate", "std.error", "df", "p.value")]),
    as.matrix(clubsando)
  )

  expect_equivalent(
    ivfit$fitted.values,
    ivcr2o$fitted.values
  )

  # CR0
  ivcr0o <- iv_robust(y ~ x | z, data = dat, clusters = clust, se_type = "CR0")
  clubsandCR0o <- clubSandwich::coef_test(ivfit, vcov = "CR0", cluster = dat$clust, test = "naive-t")

  expect_equivalent(
    as.matrix(tidy(ivcr0o)[, c("estimate", "std.error", "p.value")]),
    as.matrix(clubsandCR0o)
  )

  # Weighting classical
  ivcw <- iv_robust(y ~ x | z, data = dat, weights = w, se_type = "classical")
  ivw <- AER::ivreg(y ~ x | z, weights = w, data = dat)
  ivregsum <- summary(ivw)

  expect_equivalent(
    as.matrix(tidy(ivcw)[, c("estimate", "std.error", "statistic", "p.value")]),
    coef(ivregsum)[, 1:4]
  )

  expect_equivalent(
    ivw$fitted.values,
    ivcw$fitted.values
  )

  # HC0 weighted
  ivrw <- iv_robust(y ~ x | z, data = dat, weights = w, se_type = "HC0")
  capture_output(ivpackrobw <- ivpack::robust.se(ivw))

  expect_equivalent(
    as.matrix(tidy(ivrw)[, c("estimate", "std.error", "statistic", "p.value")]),
    ivpackrobw[, 1:4]
  )

  expect_equivalent(
    ivrw$fitted.values,
    ivcw$fitted.values
  )

  # CR2 weighted
  ivcr2wo <- iv_robust(y ~ x | z, data = dat, clusters = clust, weights = w)
  clubsandwo <- clubSandwich::coef_test(ivw, vcov = "CR2", cluster = dat$clust)

  expect_equivalent(
    as.matrix(tidy(ivcr2wo)[, c("estimate", "std.error", "df", "p.value")]),
    as.matrix(clubsandwo)
  )

  expect_equivalent(
    ivcr2wo$fitted.values,
    ivcw$fitted.values
  )

  # CR0 weighted
  ivcr0wo <- iv_robust(y ~ x | z, data = dat, clusters = clust, weights = w, se_type = "CR0")
  clubsandCR0wo <- clubSandwich::coef_test(ivw, vcov = "CR0", cluster = dat$clust, test = "naive-t")

  expect_equivalent(
    as.matrix(tidy(ivcr0wo)[, c("estimate", "std.error", "p.value")]),
    as.matrix(clubsandCR0wo)
  )

  expect_equivalent(
    ivcr0wo$fitted.values,
    ivcw$fitted.values
  )

  # Rank-deficiency
  # HC0
  dat$x1_c <- dat$x
  ivdefr <- iv_robust(y ~ x + x1_c| z + z2, data = dat, se_type = "HC0")
  ivdef <- AER::ivreg(y ~ x + x1_c| z + z2, data = dat)
  capture_output(ivdefse <- ivpack::robust.se(ivdef))

  expect_equal(
    coef(ivdefr),
    coef(ivdef)
  )

  expect_equivalent(
    as.matrix(tidy(ivdefr)[1:2, c("estimate", "std.error", "statistic", "p.value")]),
    ivdefse[, 1:4]
  )

  expect_equivalent(
    ivdefr$fitted.values,
    ivdef$fitted.values
  )

  # # Does not work if instrument is collinear with other instrument
  # ivdefri <- iv_robust(y ~ z + z2| x + x1_c, data = dat, se_type = "HC0")
  # ivdefi <- AER::ivreg(y ~ z + z2| x + x1_c, data = dat)
  # ivdefsei <- ivpack::robust.se(ivdefi)
  #
  # # No longer equal!
  # expect_equal(
  #   coef(ivdefri),
  #   coef(ivdefi)
  # )

  # expect_equivalent(
  #   as.matrix(tidy(ivdefri)[1:2, c("estimate", "std.error", "statistic", "p.value")]),
  #   ivdefsei[, 1:4]
  # )

  # Stata
  ivdefclr <- iv_robust(y ~ x + x1_c | z + z2, data = dat, clusters = clust, se_type = "stata")
  ivdefcl <- AER::ivreg(y ~ x + x1_c | z + z2, data = dat)
  capture_output(ivdefclse <- ivpack::cluster.robust.se(ivdefcl, clusterid = dat$clust))

  expect_equal(
    coef(ivdefclr),
    coef(ivdefcl)
  )

  expect_equivalent(
    as.matrix(tidy(ivdefclr)[1:2, c("estimate", "std.error")]),
    ivdefclse[, c(1, 2)]
  )


  expect_equivalent(
    ivdefclr$fitted.values,
    ivdefcl$fitted.values
  )

  # CR2
  ivdefcl2r <- iv_robust(y ~ x + x1_c | z + z2, data = dat, clusters = clust, se_type = "CR2")
  ivdefcl2 <- AER::ivreg(y ~ x + x1_c | z + z2, data = dat)
  ivdefcl2se <- clubSandwich::coef_test(ivdefcl2, vcov = "CR2", cluster = dat$clust)


  expect_equivalent(
    na.omit(as.matrix(tidy(ivdefcl2r)[, c("estimate", "std.error", "df", "p.value")])),
    na.omit(as.matrix(ivdefcl2se))
  )

  expect_equivalent(
    ivdefcl2r$fitted.values,
    ivdefcl2$fitted.values
  )

  # HC0 Weighted
  ivdefrw <- iv_robust(y ~ x + x1_c| z + z2, weights = w, data = dat, se_type = "HC0")
  ivdefw <- AER::ivreg(y ~ x + x1_c| z + z2, weights = w, data = dat)
  capture_output(ivdefsew <- ivpack::robust.se(ivdefw))

  expect_equal(
    coef(ivdefrw),
    coef(ivdefw)
  )

  expect_equivalent(
    as.matrix(tidy(ivdefrw)[1:2, c("estimate", "std.error", "statistic", "p.value")]),
    ivdefsew[, 1:4]
  )

  expect_equivalent(
    ivdefrw$fitted.values,
    ivdefw$fitted.values
  )

  # CR2 Weighted
  ivdefclrw <- iv_robust(y ~ x + x1_c | z + z2, data = dat, clusters = clust, weights = w, se_type = "CR2")
  ivdefclw <- AER::ivreg(y ~ x + x1_c | z + z2, weights = w, data = dat)
  ivdefclsew <- clubSandwich::coef_test(ivdefclw, vcov = "CR2", cluster = dat$clust)

  expect_equivalent(
    na.omit(as.matrix(tidy(ivdefclrw)[, c("estimate", "std.error", "df", "p.value")])),
    na.omit(as.matrix(ivdefclsew)[, 1:4])
  )

  expect_equivalent(
    ivdefclrw$fitted.values,
    ivdefclw$fitted.values
  )

  ivdef2clrw <- iv_robust(y ~ x + z | x + x1_c, data = dat, clusters = clust, weights = w, se_type = "CR2")
  ivdef2clw <- AER::ivreg(y ~ x + z | x + x1_c, weights = w, data = dat)
  ivdef2clsew <- clubSandwich::coef_test(ivdef2clw, vcov = "CR2", cluster = dat$clust)

  expect_equivalent(
    na.omit(as.matrix(tidy(ivdef2clrw)[, c("estimate", "std.error", "df", "p.value")])),
    na.omit(as.matrix(ivdef2clsew)[, 1:4])
  )

  expect_equivalent(
    ivdef2clrw$fitted.values,
    ivdef2clw$fitted.values
  )

  # F-stat fails properly with blocks of size 1
  set.seed(42)
  N <- 20
  dat <- data.frame(y = rnorm(N), x = rnorm(N), z = rnorm(N), bl = sample(letters, size = N, replace = T))
  ivr <- iv_robust(y ~ bl + x | bl + z, data = dat, se_type = "stata")
  expect_equivalent(
    ivr$fstatistic[1],
    NA_integer_
  )

})


test_that("iv_robust different specifications work", {
  # More instruments than endog. regressors
  ivro <- iv_robust(mpg ~ wt | hp + cyl, data = mtcars, se_type = "HC0")
  ivo <- AER::ivreg(mpg ~ wt | hp + cyl, data = mtcars)
  capture_output(ivpo <- ivpack::robust.se(ivo))
  expect_equivalent(
    as.matrix(tidy(ivro)[, c("estimate", "std.error", "statistic", "p.value")]),
    ivpo[, 1:4]
  )

  # . notation for multiple exog, doesnt work!
  # ivro <- iv_robust(mpg ~ wt + hp + vs | . - vs + cyl, data = mtcars, se_type = "HC0")
  # ivo <- AER::ivreg(mpg ~ wt + hp + vs | . - vs + cyl, data = mtcars)
  # ivpo <- ivpack::robust.se(ivo)
  # expect_equivalent(
  #   as.matrix(tidy(ivro)[, c("estimate", "std.error", "statistic", "p.value")]),
  #   ivpo[, 1:4]
  # )

  # . notation in general
  ivro <- iv_robust(mpg ~ .| ., data = mtcars, se_type = "HC0")
  ivo <- AER::ivreg(mpg ~ . | ., data = mtcars)
  capture_output(ivpo <- ivpack::robust.se(ivo))

  expect_equivalent(
    as.matrix(tidy(ivro)[, c("estimate", "std.error", "statistic", "p.value")]),
    ivpo[, 1:4]
  )

})

test_that("S3 methods", {

  ivo <- AER::ivreg(mpg ~ hp + cyl | wt + gear, data = mtcars)
  ivro <- iv_robust(mpg ~ hp + cyl | wt + gear, data = mtcars, se_type = "classical")

  expect_equal(
    vcov(ivro),
    vcov(ivo)
  )

  expect_is(
    tidy(ivro),
    "data.frame"
  )

  expect_equal(
    nrow(tidy(ivro)),
    3
  )

  summary(ivro)

  siv <- capture_output(
    summary(ivro),
    print = TRUE
  )

  expect_true(
    grepl(
      "iv\\_robust\\(formula = mpg \\~ hp \\+ cyl \\| wt \\+ gear, data = mtcars,",
      siv
    )
  )

  expect_true(
    grepl(
      "F\\-statistic\\: 33\\.73 on 2 and 29 DF,  p\\-value\\: 2\\.706e\\-08",
      siv
    )
  )

  capture_output(
    expect_equivalent(
      coef(summary(ivro)),
      print(ivro)
    )
  )

  expect_equivalent(
    ivro$fstatistic,
    summary(ivo)$waldtest[-2]
  )

  expect_equal(
    predict(ivro, newdata = mtcars),
    predict(ivo)
  )

  # no intercept
  ivo <- AER::ivreg(mpg ~ hp + cyl + 0 | wt + gear, data = mtcars)
  ivro <- iv_robust(mpg ~ hp + cyl + 0 | wt + gear, data = mtcars, se_type = "classical")

  expect_equivalent(
    ivro$fstatistic,
    summary(ivo)$waldtest[-2]
  )

})

