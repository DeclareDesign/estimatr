context("Estimator - iv_robust")

skip_if_not_installed("AER")
skip_if_not_installed("ivpack")
skip_if_not_installed("clubSandwich")

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
# xhaven::write_dta(dat, path = "~/test.dta", version = 13)

library(AER)
library(ivpack)

test_that("iv_robust matches AER + ivpack", {
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
  # Weighting
  ivcw <- iv_robust(y ~ x | z, data = dat, weights = w, se_type = "classical")
  ivw <- ivreg(y ~ x | z, weights = w, data = dat)
  ivpackrobw <- robust.se(ivw)

  expect_equivalent(
    as.matrix(tidy(ivrw)[, c("coefficients", "se", "p")]),
    ivpackrobw[, c(1, 2, 4)]
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

})
