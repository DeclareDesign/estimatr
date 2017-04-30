
df <- data.frame(Y = rnorm(100), Z = rbinom(100, 1, .5), X = rnorm(100))

lm_robust_se(Y ~ Z, data = df)

lm_robust_se(Y ~ Z + X, data = df)
lm_robust_se(Y ~ Z + X, coefficient_name = "X", data = df)
lm_robust_se(Y ~ Z + X, coefficient_name = c("Z", "X"), data = df)
lm_robust_se(Y ~ Z + X, coefficient_name = c("(Intercept)", "Z", "X"), data = df)


lm_robust_se(Y ~ Z*X, coefficient_name = "Z:X", data = df)
