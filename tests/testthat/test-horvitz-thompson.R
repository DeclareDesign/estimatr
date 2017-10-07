context("Horvitz Thompson")


test_that("Horvitz thompson matches d-i-m under certain conditions", {

  n <- 4
  dat <- data.frame(y0 = rnorm(n),
                    z = rep(0:1, each = n/2),
                    ps = rep(0.5, n))

  dat$y1 <- dat$y0 + 0.43
  dat$y <- ifelse(dat$z, dat$y1, dat$y0)

  horvitz_thompson(y ~ z,
                   ps,
                   data = dat)
  difference_in_means(y ~ z,
                      data = dat)
  # tests for missingness
  # tests for equality across declare vs. args

})
