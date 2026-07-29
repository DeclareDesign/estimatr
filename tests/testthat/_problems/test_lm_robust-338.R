# Extracted from test_lm_robust.R:338

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
d <- dat
d$f <- factor(rep(letters[1:4], length.out = n))
m <- lm_robust(y ~ x + f, data = d, fixed_effects = ~ block)
nd <- d[1:5, ]
expect_equal(unname(predict(m, nd)),
               unname(predict(lm(y ~ x + f + factor(block), data = d), nd)),
               tolerance = 1e-8)
