context("Use Model")

test_that("Use model",{


df <- data.frame(Y = rnorm(100), Z = rbinom(100, 1, .5), X = rnorm(100))

use_model(Y ~ Z, data = df, coefficient_name = c("(Intercept)", "Z"))
use_model(Y ~ Z, data = df, model = glm, coefficient_name = c("(Intercept)", "Z"))
use_model(Z ~ X, data = df, model = glm, family = binomial, coefficient_name = c("(Intercept)", "X"))
use_model(Z ~ X, data = df, model = glm, coefficient_name = c("(Intercept)", "X"))

})
