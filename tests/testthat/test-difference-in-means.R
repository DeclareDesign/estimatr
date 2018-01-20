
context("Difference in means")


test_that("DIM", {

  dat <- data.frame(Y = rnorm(100), Z = rbinom(100, 1, .5), X = rnorm(100))

  difference_in_means(Y ~ Z, data = dat)
  difference_in_means(Y ~ Z, condition1 = 0, condition2 = 1, data = dat)
  difference_in_means(Y ~ Z, condition1 = 1, condition2 = 0, data = dat)

  difference_in_means(Y ~ Z, alpha = .05, data = dat)
  difference_in_means(Y ~ Z, alpha = .10, data = dat)

  dat <- data.frame(Y = rnorm(100), Z = sample(1:3, 100, replace = TRUE), X = rnorm(100))

  difference_in_means(Y ~ Z, condition1 = 1, condition2 = 2, data = dat)
  difference_in_means(Y ~ Z, condition1 = 2, condition2 = 1, data = dat)
  difference_in_means(Y ~ Z, condition1 = 3, condition2 = 1, data = dat)
  difference_in_means(Y ~ Z, condition1 = 3, condition2 = 2, data = dat)

})

test_that("DIM Blocked", {

  dat <- data.frame(Y = rnorm(100),
                   Z = rbinom(100, 1, .5),
                   block = sample(c("A", "B", "C"), 100, replace = TRUE))

  difference_in_means(Y ~ Z, blocks = block, data = dat)
  difference_in_means(Y ~ Z, condition1 = 0, condition2 = 1, blocks = block, data = dat)
  difference_in_means(Y ~ Z, condition1 = 1, condition2 = 0, blocks = block, data = dat)

  difference_in_means(Y ~ Z, alpha = .05, blocks = block, data = dat)
  difference_in_means(Y ~ Z, alpha = .10, blocks = block, data = dat)

})

test_that("DIM same as t.test", {

  # test df correction
  dat <- data.frame(Y = rnorm(100), Z = rbinom(100, 1, .5), X = rnorm(100))

  expect_equal(
    unlist(difference_in_means(Y ~ Z, data = dat)[c('p', 'ci_lower', 'ci_upper', 'df')],
           F,
           F),
    unlist(with(dat, t.test(Y[Z==1], Y[Z==0]))[c('p.value', 'conf.int', 'parameter')],
           F,
           F)
  )


})


test_that("DIM Weighted", {

  dat <- data.frame(Y = rnorm(100),
                   Z = rbinom(100, 1, .5),
                   weights = runif(100),
                   weights2 = 1)

  difference_in_means(Y ~ Z, alpha = .10, weights = weights, data = dat)

  expect_equal(
    difference_in_means(Y ~ Z, alpha = .10, weights = weights2, data = dat),
    difference_in_means(Y ~ Z, alpha = .10, data = dat)
  )

})


test_that("DIM Clustered", {

  dat <- data.frame(weights = runif(100),
                   weights2 = 1,
                   J = rep(1:4, each = 25))

  dat$Y <- rnorm(100, mean = rep(rnorm(4, sd = sqrt(0.1)), each = 25), sd = sqrt(0.9))
  dat$Z <- as.numeric(dat$J %in% 1:2)

  difference_in_means(Y ~ Z, alpha = .10, data = dat)
  difference_in_means(Y ~ Z, alpha = .10, clusters = J, data = dat)

})

test_that("DIM Pair Matched", {

  dat <- data.frame(Y = rnorm(100),
                   Z = rbinom(100, 1, .5),
                   weights = runif(100),
                   weights2 = 1,
                   block = rep(1:50, each = 2))

  expect_error(
    difference_in_means(Y ~ Z, alpha = .05, blocks = block, data = dat),
    'both treatment'
  )

  dat$Z <- rep(0:1, 50)
  difference_in_means(Y ~ Z, alpha = .05, blocks = block, data = dat)

})


test_that("DIM Matched Pair Cluster Randomization", {

  dat <- data.frame(Y = rnorm(100),
                   block = rep(1:25, each = 4),
                   cluster = as.character(rep(1:50, each = 2)),
                   Z = rep(0:1, times = 50))

  expect_error(
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      blocks = block,
      clusters = cluster,
      data = dat
    ),
    'same treatment condition'
  )

  dat$Z <- c(rep(rep(0:1, each = 4), 12), rep(0, 4))
  expect_error(
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      blocks = block,
      clusters = cluster,
      data = dat
    ),
    'both treatment conditions'
  )

  dat$Z <- rep(rep(0:1, each = 2), 25)
  difference_in_means(
    Y ~ Z,
    alpha = .05,
    blocks = block,
    clusters = cluster,
    data = dat
  )


})

test_that("DIM Matched Pair Cluster Randomization = Matched Pair when cluster size is 1", {
  dat <- data.frame(Y = rnorm(100),
                   block = rep(1:25, each = 4),
                   cluster = 1:100,
                   Z = rep(c(0,0,1,1), times = 25))

  expect_equal(
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      blocks = block,
      clusters = cluster,
      data = dat
    ),
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      blocks = block,
      data = dat
    )
  )

})

test_that("DIM works with missingness", {
  dat <- data.frame(Y = rnorm(100),
                   block = rep(1:2, each = 50),
                   cluster = 1:100,
                   Z = rep(c(0,0,1,1), times = 25))

  ## Missingness on treatment
  dat$Z[23] <- NA

  expect_error(
    estimatr_dim_out <- difference_in_means(
      Y ~ Z,
      alpha = .05,
      blocks = block,
      data = dat
    ),
    NA
  )

  expect_equal(
    estimatr_dim_out,
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      blocks = block,
      data = dat[-23, ]
    )
  )

  ## Missingness on block
  dat$block[35] <- NA

  expect_warning(
    estimatr_missblock_dim <- difference_in_means(
      Y ~ Z,
      alpha = .05,
      blocks = block,
      data = dat
    ),
    'missingness in the block'
  )

  expect_equal(
    estimatr_missblock_dim,
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      blocks = block,
      data = dat[-c(23, 35), ]
    )
  )

  ## Missingness on cluster
  dat$cluster[1] <- NA

  expect_warning(
    estimatr_missclust_dim <- difference_in_means(
      Y ~ Z,
      alpha = .05,
      clusters = cluster,
      data = dat
    ),
    'missingness in the cluster'
  )

  expect_equal(
    estimatr_missclust_dim,
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      clusters = cluster,
      data = dat[-c(1, 23), ]
    )
  )


})


test_that("DIM works with character args", {
  dat <- data.frame(Y = rnorm(100),
                   block = rep(1:25, each = 4),
                   cluster = 1:100,
                   Z = rep(c(0,0,1,1), times = 25))

  dim_unquote <- difference_in_means(
    Y ~ Z,
    alpha = .05,
    blocks = block,
    clusters = cluster,
    data = dat
  )

  dim_quote <- difference_in_means(
    Y ~ Z,
    alpha = .05,
    blocks = "block",
    clusters = "cluster",
    data = dat
  )

  expect_equal(
    dim_unquote,
    dim_quote
  )

  expect_identical(
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      weights = cluster,
      data = dat
    ),
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      weights = "cluster",
      data = dat
    )
  )

})

test_that("DIM unbiased", {

  # Skip these tests until randomizr::obtain_permutation_matrix is pushed to CRAN
  # They have been passing locally
  check_bias <- FALSE

  if(check_bias) {
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
    declaration <- randomizr::declare_ra(N = nrow(dat))
    treatment_perms <- randomizr::obtain_permutation_matrix(declaration)

    ests <- apply(treatment_perms,
                  2,
                  function(x) {
                    dat$Z <- x
                    dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
                    dim <- difference_in_means(Y ~ Z, data = dat)
                    dim$est
                  }
    )

    expect_equivalent(
      trueSATE,
      mean(ests)
    )

    ## cluster randomized design, 5 blocks of 2
    dat$cluster <- rep(1:5, each = 2)
    declaration <- randomizr::declare_ra(N = nrow(dat),
                                         clust_var = dat$cluster)
    treatment_perms <- randomizr::obtain_permutation_matrix(declaration)

    ests <- apply(treatment_perms,
                  2,
                  function(x) {
                    dat$Z <- x
                    dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
                    dim <- difference_in_means(Y ~ Z,
                                               clusters = cluster,
                                               data = dat)
                    dim$est
                  }
    )

    expect_equivalent(
      trueSATE,
      mean(ests)
    )

    ## Matched pair design, 5 blocks of 2
    dat$blocks <- rep(1:5, each = 2)
    declaration <- randomizr::declare_ra(N = nrow(dat),
                                         block_var = dat$blocks,
                                         block_m = rep(1, 5))
    treatment_perms <- randomizr::obtain_permutation_matrix(declaration)

    ests <- apply(treatment_perms,
                  2,
                  function(x) {
                    dat$Z <- x
                    dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
                    dim <- difference_in_means(Y ~ Z,
                                               blocks = blocks,
                                               data = dat)
                    dim$est
                  }
    )

    expect_equivalent(
      trueSATE,
      mean(ests)
    )

    ## block randomized design, 2 blocks of 5
    dat$blocks <- rep(1:2, each = 5)
    declaration <- randomizr::declare_ra(N = nrow(dat),
                                         block_var = dat$blocks,
                                         block_m = c(3, 3))
    treatment_perms <- randomizr::obtain_permutation_matrix(declaration)

    ests <- apply(treatment_perms,
                  2,
                  function(x) {
                    dat$Z <- x
                    dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
                    dim <- difference_in_means(Y ~ Z,
                                               blocks = blocks,
                                               data = dat)
                    dim$est
                  }
    )

    expect_equivalent(
      trueSATE,
      mean(ests)
    )

    ## cluster matched pair, different sized blocks
    dat$blocks <- rep(1:3, times = c(4, 4, 2))
    dat$clusters <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 6)
    declaration <- randomizr::declare_ra(N = nrow(dat),
                                         block_var = dat$blocks,
                                         clust_var = dat$clusters)
    treatment_perms <- randomizr::obtain_permutation_matrix(declaration)

    ests <- apply(treatment_perms,
                  2,
                  function(x) {
                    dat$Z <- x
                    dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
                    dim <- difference_in_means(Y ~ Z,
                                               blocks = blocks,
                                               clusters = clusters,
                                               data = dat)
                    dim$est
                  }
    )

    expect_equivalent(
      trueSATE,
      mean(ests)
    )
  }

})

test_that("DIM matches lm_robust under certain conditions", {

  dat <- data.frame(Y = rnorm(400))
  ## DIM and lm_robust agree without clustering except for DoF because DIM uses Satterthwaite approx
  dat$z <- c(0, 1)
  lm_o <- lm_robust(Y ~ z, data = dat, coefficient_name = "z")
  dim_o <- difference_in_means(Y ~ z, data = dat)
  expect_equivalent(
    tidy(lm_o)[, 1:3],
    tidy(dim_o)[, 1:3]
  )

  ## DIM and lm_robust agree with clustering becuase DIM just uses lm_robust w/ CR2!
  dat$cl_diff_size <- sample(100, size = 400, replace = TRUE)
  dat$z_clustered <- as.numeric(dat$cl_diff_size <= 50)

  lm_cl_o <- lm_robust(Y ~ z_clustered, clusters = cl_diff_size, data = dat, coefficient_name = "z_clustered")
  dim_cl_o <- difference_in_means(Y ~ z_clustered, clusters = cl_diff_size, data = dat)
  expect_equivalent(
    tidy(lm_cl_o),
    tidy(dim_cl_o)
  )

  ## Blocked design is equivalent to lm_lin
  dat$bl <- rep(1:25, each = 16)
  dat$z_blocked <- rep(c(0, 1), each = 2)

  lm_bl_o <- lm_lin(Y ~ z_blocked, ~ factor(bl), data = dat, coefficient_name = "z_blocked")
  dim_bl_o <- difference_in_means(Y ~ z_blocked, blocks = bl, data = dat)

  # Not identical since row name of lm_bl_o is 2 due to the intercept
  expect_equivalent(
    tidy(lm_bl_o),
    tidy(dim_bl_o)
  )

  ## Block-clustered is equivalent to lm_lin
  ## (and indeed uses lm_robust machinery for the ests and ses, the DoF are equivalent with equal clusters
  ## by design
  dat$cl <- rep(1:200, each = 2)
  lm_blcl_o <- lm_lin(Y ~ z_blocked, ~ factor(bl), data = dat, clusters = cl, coefficient_name = "z_blocked")
  dim_blcl_o <- difference_in_means(Y ~ z_blocked, blocks = bl, clusters = cl, data = dat)

  # Not identical since row name of lm_bl_o is 2 due to the intercept
  expect_equivalent(
    tidy(lm_blcl_o),
    tidy(dim_blcl_o)
  )

})
