context("Estimator - iv_robust")

N <- 20
dat <- data.frame(
  y = rnorm(N),
  x = rnorm(N),
  x2 = rnorm(N),
  z2 = rnorm(N),
  clust = sample(letters[1:3], size = N, replace = TRUE)
)
dat$z <- dat$x * 0.5 + rnorm(N)
haven::write_dta(dat, path = "~/test.dta", version = 13)

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

  # Stata defaults to HC0 as well
  ivro <- iv_robust(y ~ x | z, data = dat, se_type = "HC0")
  ivpackrob <- robust.se(ivfit)

  expect_equivalent(
    as.matrix(tidy(ivro)[, c("coefficients", "se", "p")]),
    ivpackrob[, c(1, 2, 4)]
  )

  # TODO "Stata" clustered SEs are actually the same as ivpack
  # but Stata actually uses CR0, have to fix and discuss options
  ivclusto <- iv_robust(y ~ x | z, data = dat, se_type = "stata", clusters = clust)
  ivpackclust <- cluster.robust.se(ivfit, dat$clust)

  # Our p-values are bigger (must be using less conservative DF, we use J - 1 which
  # is what stata uses for clusters in OLS)
  expect_equivalent(
    as.matrix(tidy(ivclusto)[, c("coefficients", "se")]),
    ivpackclust[, c(1, 2)]
  )
})
