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

  skip_if_not_installed("AER")
  skip_if_not_installed("ivpack")
  skip_if_not_installed("clubSandwich")

  library(AER)
  library(ivpack)

  ivco <- iv_robust(y ~ x | z, data = dat, se_type = "classical")
  ivfit <- ivreg(y ~ x | z, data = dat)
  ivo <- summary(ivfit)

  expect_equivalent(
    as.matrix(tidy(ivco)[, c("coefficients", "se", "p")]),
    ivo$coefficients[, c(1, 2, 4)]
  )
  # Same as stata if you specify `small` as a stata option
  # which applies the N / N-k finite sample correction

  # Stata defaults to HC0 as well, but does HC1 with `small`
  ivro <- iv_robust(y ~ x | z, data = dat, se_type = "HC0")
  ivpackrob <- robust.se(ivfit)

  expect_equivalent(
    as.matrix(tidy(ivro)[, c("coefficients", "se", "p")]),
    ivpackrob[, c(1, 2, 4)]
  )

  # "Stata" clustered SEs are CR0, but they are the same as below with `small`
  ivclusto <- iv_robust(y ~ x | z, data = dat, se_type = "stata", clusters = clust)
  ivpackclust <- cluster.robust.se(ivfit, dat$clust)

  # Our p-values are bigger (ivpack is be using less conservative DF, we use J - 1 which
  # is what stata uses for clusters w/ `small` and in OLS)
  expect_equivalent(
    as.matrix(tidy(ivclusto)[, c("coefficients", "se")]),
    ivpackclust[, c(1, 2)]
  )

  # CR2
  ivcr2o <- iv_robust(y ~ x | z, data = dat, clusters = clust)
  clubsando <- clubSandwich::coef_test(ivfit, vcov = "CR2", cluster = dat$clust)

  expect_equivalent(
    as.matrix(tidy(ivcr2o)[, c("coefficients", "se", "df", "p")]),
    as.matrix(clubsando)
  )
  # CR0
  ivcr0o <- iv_robust(y ~ x | z, data = dat, clusters = clust, se_type = "CR0")
  clubsandCR0o <- clubSandwich::coef_test(ivfit, vcov = "CR0", cluster = dat$clust, test = "naive-t")

  expect_equivalent(
    as.matrix(tidy(ivcr0o)[, c("coefficients", "se", "p")]),
    as.matrix(clubsandCR0o)
  )

  # Weighting classical
  ivcw <- iv_robust(y ~ x | z, data = dat, weights = w, se_type = "classical")
  ivw <- ivreg(y ~ x | z, weights = w, data = dat)
  ivregsum <- summary(ivw)

  expect_equivalent(
    as.matrix(tidy(ivcw)[, c("coefficients", "se", "p")]),
    ivregsum$coefficients[, c(1, 2, 4)]
  )

  # HC1 weighted
  ivrw <- iv_robust(y ~ x | z, data = dat, weights = w, se_type = "HC0")
  ivpackrobw <- robust.se(ivw)

  expect_equivalent(
    as.matrix(tidy(ivrw)[, c("coefficients", "se", "p")]),
    ivpackrobw[, c(1, 2, 4)]
  )

  # CR2 weighted
  ivcr2wo <- iv_robust(y ~ x | z, data = dat, clusters = clust, weights = w)
  clubsandwo <- clubSandwich::coef_test(ivw, vcov = "CR2", cluster = dat$clust)

  expect_equivalent(
    as.matrix(tidy(ivcr2wo)[, c("coefficients", "se", "df", "p")]),
    as.matrix(clubsandwo)
  )

  # CR0 weighted
  ivcr0wo <- iv_robust(y ~ x | z, data = dat, clusters = clust, weights = w, se_type = "CR0")
  clubsandCR0wo <- clubSandwich::coef_test(ivw, vcov = "CR0", cluster = dat$clust, test = "naive-t")

  expect_equivalent(
    as.matrix(tidy(ivcr0wo)[, c("coefficients", "se", "p")]),
    as.matrix(clubsandCR0wo)
  )

  # Rank-deficiency
  # HC0
  dat$x1_c <- dat$x
  ivdefr <- iv_robust(y ~ x + x1_c| z + z2, data = dat, se_type = "HC0")
  ivdef <- ivreg(y ~ x + x1_c| z + z2, data = dat)
  ivdefse <- robust.se(ivdef)

  expect_equal(
    ivdefr$coefficients,
    ivdef$coefficients
  )

  expect_equivalent(
    as.matrix(tidy(ivdefr)[1:2, c("coefficients", "se", "p")]),
    ivdefse[, c(1, 2, 4)]
  )

  # Also works as instrument
  ivdefri <- iv_robust(y ~ z + z2| x + x1_c, data = dat, se_type = "HC0")
  ivdefi <- ivreg(y ~ z + z2| x + x1_c, data = dat)
  ivdefsei <- robust.se(ivdefi)

  expect_equal(
    ivdefri$coefficients,
    ivdefi$coefficients
  )

  expect_equivalent(
    as.matrix(tidy(ivdefri)[1:2, c("coefficients", "se", "p")]),
    ivdefsei[, c(1, 2, 4)]
  )

  # Stata
  ivdefclr <- iv_robust(y ~ x + x1_c | z + z2, data = dat, clusters = clust, se_type = "stata")
  ivdefcl <- ivreg(y ~ x + x1_c | z + z2, data = dat)
  ivdefclse <- cluster.robust.se(ivdefcl, clusterid = dat$clust)

  expect_equal(
    ivdefclr$coefficients,
    ivdefcl$coefficients
  )

  expect_equivalent(
    as.matrix(tidy(ivdefclr)[1:2, c("coefficients", "se")]),
    ivdefclse[, c(1, 2)]
  )

  # CR2
  ivdefcl2r <- iv_robust(y ~ x + x1_c | z + z2, data = dat, clusters = clust, se_type = "CR2")
  ivdefcl2 <- ivreg(y ~ x + x1_c | z + z2, data = dat)
  ivdefcl2se <- clubSandwich::coef_test(ivdefcl2, vcov = "CR2", cluster = dat$clust)


  expect_equivalent(
    as.matrix(tidy(ivdefcl2r)[1:2, c("coefficients", "se", "df", "p")]),
    as.matrix(ivdefcl2se)
  )

  # HC0 Weighted
  ivdefrw <- iv_robust(y ~ x + x1_c| z + z2, weights = w, data = dat, se_type = "HC0")
  ivdefw <- ivreg(y ~ x + x1_c| z + z2, weights = w, data = dat)
  ivdefsew <- robust.se(ivdefw)

  expect_equal(
    ivdefrw$coefficients,
    ivdefw$coefficients
  )

  expect_equivalent(
    as.matrix(tidy(ivdefrw)[1:2, c("coefficients", "se", "p")]),
    ivdefsew[, c(1, 2, 4)]
  )

  # CR2 Weighted
  ivdefclrw <- iv_robust(y ~ x + x1_c | z + z2, data = dat, clusters = clust, weights = w, se_type = "CR2")
  ivdefclw <- ivreg(y ~ x + x1_c | z + z2, weights = w, data = dat)
  ivdefclsew <- clubSandwich::coef_test(ivdefclw, vcov = "CR2", cluster = dat$clust)

  expect_equivalent(
    as.matrix(tidy(ivdefclrw)[1:2, c("coefficients", "se", "p")]),
    as.matrix(ivdefclsew)[, c(1, 2, 4)]
  )
})


test_that("iv_robust different specifications work", {
  skip_if_not_installed("AER")
  skip_if_not_installed("ivpack")
  library(AER)

  # More instruments than endog. regressors
  ivro <- iv_robust(mpg ~ wt | hp + cyl, data = mtcars, se_type = "HC0")
  ivo <- ivreg(mpg ~ wt | hp + cyl, data = mtcars)
  ivpo <- robust.se(ivo)
  expect_equivalent(
    as.matrix(tidy(ivro)[, c("coefficients", "se", "p")]),
    ivpo[, c(1, 2, 4)]
  )

  # . notation for multiple exog, doesnt work!
  # ivro <- iv_robust(mpg ~ wt + hp + vs | . - vs + cyl, data = mtcars, se_type = "HC0")
  # ivo <- ivreg(mpg ~ wt + hp + vs | . - vs + cyl, data = mtcars)
  # ivpo <- robust.se(ivo)
  # expect_equivalent(
  #   as.matrix(tidy(ivro)[, c("coefficients", "se", "p")]),
  #   ivpo[, c(1, 2, 4)]
  # )

  # . notation in general
  ivro <- iv_robust(mpg ~ .| ., data = mtcars, se_type = "HC0")
  ivo <- ivreg(mpg ~ . | ., data = mtcars)
  ivpo <- robust.se(ivo)
  expect_equivalent(
    as.matrix(tidy(ivro)[, c("coefficients", "se", "p")]),
    ivpo[, c(1, 2, 4)]
  )

})

test_that("S3 methods", {

  skip_if_not_installed("AER")

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

  capture_output(
    expect_equivalent(
      summary(ivro)$coefficients,
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
  ivo <- AER::ivreg(mpg ~ hp + cyl +0 | wt + gear, data = mtcars)
  ivro <- iv_robust(mpg ~ hp + cyl +0| wt + gear, data = mtcars, se_type = "classical")

  expect_equivalent(
    ivro$fstatistic,
    summary(ivo)$waldtest[-2]
  )

})
