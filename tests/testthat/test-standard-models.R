
df <- data.frame(Y = rnorm(100), Z = rbinom(100, 1, .5), X = rnorm(100))

use_model(Y ~ Z, data = df, coefficient_name = c("(Intercept)", "Z"))
