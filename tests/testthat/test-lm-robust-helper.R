# N <- 10000
# dat <- data.frame(y = rnorm(N), x1 = rnorm(N), x2 = rnorm(N))
# X <- model.matrix.default(~ x1 + x2, data = dat)
# y <- dat$y
# fit <- lm(y ~ x1 + x2, data = dat)
#
# # Speed Tests -------------------------------------------------------------
#
# system.time(replicate(500, lm_robust_helper(X = X, y = dat$y, type = "classical")))
#
# system.time(replicate(500, lm_robust_helper(X = X, y = dat$y, type = "HC0")))
# system.time(replicate(500, lm_robust_helper(X = X, y = dat$y, type = "HC1")))
# system.time(replicate(500, lm_robust_helper(X = X, y = dat$y, type = "HC2")))
# system.time(replicate(500, lm_robust_helper(X = X, y = dat$y, type = "HC3")))
#
#
