
context("Difference in means")


test_that("DIM", {

  df <- data.frame(Y = rnorm(100), Z = rbinom(100, 1, .5), X = rnorm(100))

  difference_in_means(Y ~ Z, data = df)
  difference_in_means(Y ~ Z, condition1 = 0, condition2 = 1, data = df)
  difference_in_means(Y ~ Z, condition1 = 1, condition2 = 0, data = df)

  difference_in_means(Y ~ Z, alpha = .05, data = df)
  difference_in_means(Y ~ Z, alpha = .10, data = df)

  df <- data.frame(Y = rnorm(100), Z = sample(1:3, 100, replace = TRUE), X = rnorm(100))

  difference_in_means(Y ~ Z, data = df)
  difference_in_means(Y ~ Z, condition1 = 1, condition2 = 2, data = df)
  difference_in_means(Y ~ Z, condition1 = 2, condition2 = 1, data = df)
  difference_in_means(Y ~ Z, condition1 = 3, condition2 = 1, data = df)
  difference_in_means(Y ~ Z, condition1 = 3, condition2 = 2, data = df)

})

test_that("DIM Blocked", {

  df <- data.frame(Y = rnorm(100),
                   Z = rbinom(100, 1, .5),
                   block = sample(c("A", "B", "C"), 100, replace = TRUE))
  difference_in_means(Y ~ Z, block_variable_name = block, data = df)
  difference_in_means(Y ~ Z, condition1 = 0, condition2 = 1, block_variable_name = block, data = df)
  difference_in_means(Y ~ Z, condition1 = 1, condition2 = 0, block_variable_name = block, data = df)

  difference_in_means(Y ~ Z, alpha = .05, block_variable_name = block, data = df)
  difference_in_means(Y ~ Z, alpha = .10, block_variable_name = block, data = df)

})

test_that("DIM Weighted", {

  df <- data.frame(Y = rnorm(100),
                   Z = rbinom(100, 1, .5),
                   weights = runif(100),
                   weights2 = 1)

  difference_in_means(Y ~ Z, alpha = .10, weights = weights, data = df)
  difference_in_means(Y ~ Z, alpha = .10, weights = weights2, data = df)
  difference_in_means(Y ~ Z, alpha = .10, data = df)

})


test_that("DIM Clustered", {

  df <- data.frame(weights = runif(100),
                   weights2 = 1,
                   J = rep(1:4, each = 25))

  df$Y <- rnorm(100, mean = rep(rnorm(4, sd = sqrt(0.1)), each = 25), sd = sqrt(0.9))
  df$Z <- as.numeric(df$J %in% 1:2)

  difference_in_means(Y ~ Z, alpha = .10, data = df)
  difference_in_means(Y ~ Z, alpha = .10, cluster_variable_name = J, data = df)

})

test_that("DIM Pair Matched", {

  df <- data.frame(Y = rnorm(100),
                   Z = rbinom(100, 1, .5),
                   weights = runif(100),
                   weights2 = 1,
                   block = rep(1:50, each = 2))

  difference_in_means(Y ~ Z, alpha = .05, block_variable_name = block, data = df)

})
