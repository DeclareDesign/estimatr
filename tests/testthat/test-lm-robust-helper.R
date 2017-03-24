library(sandwich)
N <- 10000
dat <- data.frame(y = rnorm(N), x1 = rnorm(N), x2 = rnorm(N))
X <- model.matrix.default(~ x1 + x2, data = dat)
y <- dat$y
fit <- lm(y ~ x1 + x2, data = dat)



# confirm calcs -----------------------------------------------------------

# Classical
lm_robust_helper(y = y, X = X, type = "classical")
vcov(fit)

# HC0
vcovHC(fit, type = "HC0")
lm_robust_helper(y = y, X = X, type = "HC0")

# HC1
vcovHC(fit, type = "HC1")
lm_robust_helper(y = y, X = X, type = "HC1")

#HC2
vcovHC(fit, type = "HC2")
lm_robust_helper(y = y, X = X, type = "HC2")

#HC3
vcovHC(fit, type = "HC3")
lm_robust_helper(y = y, X = X, type = "HC3")

# Speed Tests -------------------------------------------------------------

system.time(replicate(500, RcppArmadillo::fastLm(X = X, y = dat$y)))
system.time(replicate(500, RcppArmadillo::fastLmPure(X = X, y = dat$y)))
system.time(replicate(500, lm_robust_helper(X = X, y = dat$y, type = "classical")))

system.time(replicate(500, lm_robust_helper(X = X, y = dat$y, type = "HC0")))
system.time(replicate(500, lm_robust_helper(X = X, y = dat$y, type = "HC1")))
system.time(replicate(500, lm_robust_helper(X = X, y = dat$y, type = "HC2")))
system.time(replicate(500, lm_robust_helper(X = X, y = dat$y, type = "HC3")))






