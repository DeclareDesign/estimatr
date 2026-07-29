# Extracted from test_lm_robust.R:346

# prequel ----------------------------------------------------------------------
library(estimatrZero)
set.seed(42)
n <- 200
dat <- data.frame(
  y = rnorm(n),
  x = rnorm(n),
  z = rbinom(n, 1, 0.5),
  cl = rep(1:20, 10),
  block = rep(1:20, each = 10),
  w = runif(n, 0.5, 2)
)
dat$z_block <- rep(rep(c(0L, 1L), 5L), 20L)

# test -------------------------------------------------------------------------
m <- lm_robust(y ~ x, data = dat, fixed_effects = ~ block)
nd <- dat[1, ]
nd$block <- 999
expect_error(predict(m, nd), "new levels")
