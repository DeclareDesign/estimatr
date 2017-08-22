
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

test_that("DIM unbiased", {

  dat <- data.frame(i = 1:10,
                    Y0 = c(2.1, 3.5, -131.2, -1.3, -4,
                           0.1, 8.1, -1.3, 1.1, 9.1),
                    Y1 = c(2.6, 3.0, -132, -0.7, -3.3,
                           0.5, 24.3, -1, 1.6, 0.3))

  # True SATE = 0.91
  trueSATE <- mean(dat$Y1) - mean(dat$Y0)

  ## Complete Randomization
  # True se(SATE_hat)
  true_seSATE <- sqrt( (var(dat$Y0) + var(dat$Y1) + 2 * cov(dat$Y0, dat$Y1)) / (10 - 1))
  declaration <- declare_ra(N = nrow(dat))
  treatment_perms <- obtain_permutation_matrix(declaration)
  ests <- matrix(NA,
                 nrow = ncol(treatment_perms),
                 ncol = 2)

  for (i in 1:nrow(ests)) {
    dat$Z <- treatment_perms[, i]
    dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
    dim <- difference_in_means(Y ~ Z, data = dat)
    ests[i, ] <- c(dim$est, dim$se)
  }

  expect_equivalent(
    trueSATE,
    mean(ests[, 1])
  )

  ## cluster randomized design, 5 blocks of 2
  dat$cluster <- rep(1:5, each = 2)
  declaration <- declare_ra(N = nrow(dat),
                            clust_var = dat$cluster)
  treatment_perms <- obtain_permutation_matrix(declaration)
  ests <- numeric(ncol(treatment_perms))

  for (i in 1:length(ests)) {
    dat$Z <- treatment_perms[, i]
    dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
    dim <- difference_in_means(Y ~ Z,
                               cluster_variable_name = cluster,
                               data = dat)
    ests[i] <- dim$est
  }

  expect_equivalent(
    trueSATE,
    mean(ests)
  )

  ## Matched pair design, 5 blocks of 2
  dat$blocks <- rep(1:5, each = 2)
  declaration <- declare_ra(N = nrow(dat),
                            block_var = dat$blocks,
                            block_m = rep(1, 5))
  treatment_perms <- obtain_permutation_matrix(declaration)
  ests <- numeric(ncol(treatment_perms))

  for (i in 1:length(ests)) {
    dat$Z <- treatment_perms[, i]
    dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
    dim <- difference_in_means(Y ~ Z,
                               block_variable_name = blocks,
                               data = dat)
    ests[i] <- dim$est
  }

  expect_equivalent(
    trueSATE,
    mean(ests)
  )

  ## block randomized design, 2 blocks of 5
  dat$blocks <- rep(1:2, each = 5)
  declaration <- declare_ra(N = nrow(dat),
                            block_var = dat$blocks,
                            block_m = c(3, 3))
  treatment_perms <- obtain_permutation_matrix(declaration)
  ests <- numeric(ncol(treatment_perms))

  for (i in 1:length(ests)) {
    dat$Z <- treatment_perms[, i]
    dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
    dim <- difference_in_means(Y ~ Z,
                               block_variable_name = blocks,
                               data = dat)
    ests[i] <- dim$est
  }

  expect_equivalent(
    trueSATE,
    mean(ests)
  )

  ## cluster matched pair, 2 blocks of 5
  dat$blocks <- rep(1:3, times = c(4, 4, 2))
  dat$clusters <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 6)
  declaration <- declare_ra(N = nrow(dat),
                            block_var = dat$blocks,
                            clust_var = dat$clusters)
  treatment_perms <- obtain_permutation_matrix(declaration)
  ests <- numeric(ncol(treatment_perms))

  for (i in 1:length(ests)) {
    dat$Z <- treatment_perms[, i]
    dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
    dim <- difference_in_means(Y ~ Z,
                               block_variable_name = blocks,
                               cluster_variable_name = clusters,
                               data = dat)
    ests[i] <- dim$est
  }

  expect_equivalent(
    trueSATE,
    mean(ests)
  )


})
