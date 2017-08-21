
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

  expect_equal(
    difference_in_means(Y ~ Z, alpha = .10, weights = weights2, data = df),
    difference_in_means(Y ~ Z, alpha = .10, data = df)
  )

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

  expect_error(
    difference_in_means(Y ~ Z, alpha = .05, block_variable_name = block, data = df),
    'both treatment'
  )

  df$Z <- rep(0:1, 50)
  difference_in_means(Y ~ Z, alpha = .05, block_variable_name = block, data = df)

})


test_that("DIM Matched Pair Cluster Randomization", {

  df <- data.frame(Y = rnorm(100),
                   block = rep(1:25, each = 4),
                   cluster = as.character(rep(1:50, each = 2)),
                   Z = rep(0:1, times = 50))

  expect_error(
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      block_variable_name = block,
      cluster_variable_name = cluster,
      data = df
    ),
    'same treatment condition'
  )

  df$Z <- c(rep(rep(0:1, each = 4), 12), rep(0, 4))
  expect_error(
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      block_variable_name = block,
      cluster_variable_name = cluster,
      data = df
    ),
    'both treatment conditions'
  )

  df$Z <- rep(rep(0:1, each = 2), 25)
  difference_in_means(
    Y ~ Z,
    alpha = .05,
    block_variable_name = block,
    cluster_variable_name = cluster,
    data = df
  )


})

test_that("DIM Matched Pair Cluster Randomization = Matched Pair when cluster size is 1", {
  df <- data.frame(Y = rnorm(100),
                   block = rep(1:25, each = 4),
                   cluster = 1:100,
                   Z = rep(c(0,0,1,1), times = 25))

  expect_equal(
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      block_variable_name = block,
      cluster_variable_name = cluster,
      data = df
    ),
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      block_variable_name = block,
      data = df
    )
  )

})

test_that("DIM works with missingness", {
  df <- data.frame(Y = rnorm(100),
                   block = rep(1:2, each = 50),
                   cluster = 1:100,
                   Z = rep(c(0,0,1,1), times = 25))

  ## Missingness on treatment
  df$Z[23] <- NA

  expect_error(
    estimatr_dim_out <- difference_in_means(
      Y ~ Z,
      alpha = .05,
      block_variable_name = block,
      data = df
    ),
    NA
  )

  expect_equal(
    estimatr_dim_out,
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      block_variable_name = block,
      data = df[-23, ]
    )
  )

  ## Missingness on block
  df$block[35] <- NA

  expect_warning(
    estimatr_missblock_dim <- difference_in_means(
      Y ~ Z,
      alpha = .05,
      block_variable_name = block,
      data = df
    ),
    'missingness in the block'
  )

  expect_equal(
    estimatr_missblock_dim,
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      block_variable_name = block,
      data = df[-c(23, 35), ]
    )
  )

  ## Missingness on cluster
  df$cluster[1] <- NA

  expect_warning(
    estimatr_missclust_dim <- difference_in_means(
      Y ~ Z,
      alpha = .05,
      cluster_variable_name = cluster,
      data = df
    ),
    'missingness in the cluster'
  )

  expect_equal(
    estimatr_missclust_dim,
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      cluster_variable_name = cluster,
      data = df[-c(1, 23), ]
    )
  )


})
