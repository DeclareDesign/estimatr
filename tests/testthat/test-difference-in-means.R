
context("Estimator - difference_in_means")


test_that("DIM", {
  dat <- data.frame(Y = rnorm(100), Z = sample(1:3, 100, replace = TRUE), X = rnorm(100))

  difference_in_means(Y ~ Z, condition1 = 1, condition2 = 2, data = dat)
  difference_in_means(Y ~ Z, condition1 = 2, condition2 = 1, data = dat)
  difference_in_means(Y ~ Z, condition1 = 3, condition2 = 1, data = dat)
  dimo <- difference_in_means(Y ~ Z, condition1 = 3, condition2 = 2, data = dat)

  expect_equal(
    dimo$design,
    "Standard"
  )
})

test_that("DIM arguments parsed correctly", {
  dat <- data.frame(Y = rnorm(100), Z = rbinom(100, 1, .5), X = rnorm(100))

  expect_equivalent(
    as.matrix(tidy(difference_in_means(
      Y ~ Z, data = dat, ci = FALSE
    ))[, c("p", "ci_lower", "ci_upper")]),
    matrix(NA, nrow = 1, ncol = 3)
  )

  expect_error(
    difference_in_means(Y ~ Z + X, data = dat),
    "must have only one variable on the right-hand side"
  )

  dat$bl <- rep(1:10, each = 10)
  dat$bad_cl <- rep(1:10, times = 10)
  expect_error(
    difference_in_means(Y ~ Z, blocks = bl, clusters = bad_cl, data = dat),
    "All `clusters` must be contained within `blocks`"
  )

  dat$bad_bl <- c(1, rep(2:10, length.out = 99))
  expect_error(
    difference_in_means(Y ~ Z, blocks = bad_bl, data = dat),
    "All `blocks` must have multiple units"
  )

  dat$bad_mp <- rep(1:50, each = 2)
  dat$bad_mp[dat$bad_mp == 50] <- 49
  expect_error(
    difference_in_means(Y ~ Z, blocks = bad_mp, data = dat),
    "`blocks` must either all have two units/`clusters`"
  )

  expect_error(
    difference_in_means(Y ~ Z + X, data = dat),
    "must have only one variable on the right-hand side"
  )

  # not matched pair but has some blocks with only 1 treated
  bl <- rep(1:2, each = 4)
  z <- c(1, 0, 0, 0, 1, 1, 0, 0)
  y <- rnorm(8)
  expect_error(
    difference_in_means(y ~ z, blocks = bl),
    "Must have least two treated/control"
  )
})

test_that("DIM Blocked", {
  dat <- data.frame(
    Y = rnorm(100),
    Z = rbinom(100, 1, .5),
    block = sample(c("A", "B", "C"), 100, replace = TRUE)
  )

  difference_in_means(Y ~ Z, blocks = block, data = dat)
  dim_normal <- difference_in_means(Y ~ Z, condition1 = 0, condition2 = 1, blocks = block, data = dat)
  dim_reverse <- difference_in_means(Y ~ Z, condition1 = 1, condition2 = 0, blocks = block, data = dat)

  expect_equal(
    tidy(dim_normal)[c("coefficients", "se")],
    tidy(dim_reverse)[c("coefficients", "se")] * c(-1, 1)
  )

  difference_in_means(Y ~ Z, alpha = .05, blocks = block, data = dat)
  difference_in_means(Y ~ Z, alpha = .10, blocks = block, data = dat)

  expect_equal(
    dim_normal$design,
    "Blocked"
  )
})

test_that("DIM same as t.test", {

  # test df correction
  dat <- data.frame(Y = rnorm(100), Z = rbinom(100, 1, .5), X = rnorm(100))

  expect_equal(
    unlist(
      difference_in_means(Y ~ Z, data = dat)[c("p", "ci_lower", "ci_upper", "df")],
      F,
      F
    ),
    unlist(
      with(dat, t.test(Y[Z == 1], Y[Z == 0]))[c("p.value", "conf.int", "parameter")],
      F,
      F
    )
  )
})


test_that("DIM Weighted", {
  n <- 100
  dat <- data.frame(y = rnorm(n), z = 0:1, w = 1, bl = rep(1:10, each = 10))
  dimw <- difference_in_means(y ~ z, weights = w, data = dat)
  difference_in_means(y ~ z, data = dat)

  dimbw <- difference_in_means(y ~ z, weights = w, blocks = bl, data = dat)
  difference_in_means(y ~ z, blocks = bl, data = dat)

  expect_equal(
    dimw$design,
    "Standard (weighted)"
  )

  expect_equal(
    dimbw$design,
    "Blocked (weighted)"
  )
})


test_that("DIM Clustered", {
  dat <- data.frame(
    weights = runif(100),
    weights2 = 1,
    J = rep(1:4, each = 25)
  )

  dat$Y <- rnorm(100, mean = rep(rnorm(4, sd = sqrt(0.1)), each = 25), sd = sqrt(0.9))
  dat$Z <- as.numeric(dat$J %in% 1:2)

  difference_in_means(Y ~ Z, alpha = .05, data = dat)
  dim_05 <- difference_in_means(Y ~ Z, alpha = .05, clusters = J, data = dat)
  dim_10 <- difference_in_means(Y ~ Z, alpha = .10, clusters = J, data = dat)

  expect_true(dim_05$ci_lower < dim_10$ci_lower)

  expect_equal(
    dim_10$design,
    "Clustered"
  )
})

test_that("DIM Pair Matched", {
  dat <- data.frame(
    Y = rnorm(100),
    Z = rbinom(100, 1, .5),
    weights = runif(100),
    weights2 = 1,
    block = rep(1:50, each = 2)
  )

  expect_error(
    difference_in_means(Y ~ Z, alpha = .05, blocks = block, data = dat),
    "both treatment"
  )

  dat$Z <- rep(0:1, 50)
  dim_mp <- difference_in_means(Y ~ Z, alpha = .05, blocks = block, data = dat)

  expect_equal(
    dim_mp$design,
    "Matched-pair"
  )
})


test_that("DIM Matched Pair Cluster Randomization", {
  dat <- data.frame(
    Y = rnorm(100),
    block = rep(1:25, each = 4),
    cluster = as.character(rep(1:50, each = 2)),
    Z = rep(0:1, times = 50)
  )

  expect_error(
    difference_in_means(
      Y ~ Z,
      alpha = .05,
      blocks = block,
      clusters = cluster,
      data = dat
    ),
    "same treatment condition"
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
    "both treatment conditions"
  )

  dat$Z <- rep(rep(0:1, each = 2), 25)
  dim_mpc <- difference_in_means(
    Y ~ Z,
    alpha = .05,
    blocks = block,
    clusters = cluster,
    data = dat
  )

  expect_equal(
    dim_mpc$design,
    "Matched-pair clustered"
  )
})

test_that("DIM Matched Pair Cluster Randomization = Matched Pair when cluster size is 1", {
  dat <- data.frame(
    Y = rnorm(100),
    block = rep(1:25, each = 4),
    cluster = 1:100,
    Z = rep(c(0, 0, 1, 1), times = 25)
  )

  expect_equal(
    tidy(difference_in_means(
      Y ~ Z,
      alpha = .05,
      blocks = block,
      clusters = cluster,
      data = dat
    )),
    tidy(difference_in_means(
      Y ~ Z,
      alpha = .05,
      blocks = block,
      data = dat
    ))
  )
})

test_that("DIM works with missingness", {
  dat <- data.frame(
    Y = rnorm(100),
    block = rep(1:2, each = 50),
    cluster = 1:100,
    Z = rep(c(0, 0, 1, 1), times = 25)
  )

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
    "missingness in the block"
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
    "missingness in the cluster"
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
  dat <- data.frame(
    Y = rnorm(100),
    block = rep(1:25, each = 4),
    cluster = 1:100,
    Z = rep(c(0, 0, 1, 1), times = 25)
  )

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

  expect_equivalent(
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
  dat <- data.frame(
    i = 1:10,
    Y0 = c(
      2.1, 3.5, -131.2, -1.3, -4,
      0.1, 8.1, -1.3, 1.1, 9.1
    ),
    Y1 = c(
      2.6, 3.0, -132, -0.7, -3.3,
      0.5, 24.3, -1, 1.6, 0.3
    )
  )

  # True SATE = 0.91
  trueSATE <- mean(dat$Y1) - mean(dat$Y0)

  ## Complete Randomization
  # True se(SATE_hat)
  true_seSATE <- sqrt((var(dat$Y0) + var(dat$Y1) + 2 * cov(dat$Y0, dat$Y1)) / (10 - 1))
  declaration <- randomizr::declare_ra(N = nrow(dat))
  treatment_perms <- randomizr::obtain_permutation_matrix(declaration)

  ests <- apply(
    treatment_perms,
    2,
    function(x) {
      dat$Z <- x
      dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
      dim <- difference_in_means(Y ~ Z, data = dat)
      dim$coefficients
    }
  )

  expect_equivalent(
    trueSATE,
    mean(ests)
  )

  ## cluster randomized design, 5 blocks of 2
  dat$cluster <- rep(1:5, each = 2)
  declaration <- randomizr::declare_ra(
    N = nrow(dat),
    clusters = dat$cluster
  )
  treatment_perms <- randomizr::obtain_permutation_matrix(declaration)

  ests <- apply(
    treatment_perms,
    2,
    function(x) {
      dat$Z <- x
      dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
      dim <- difference_in_means(
        Y ~ Z,
        clusters = cluster,
        data = dat
      )
      dim$coefficients
    }
  )

  expect_equivalent(
    trueSATE,
    mean(ests)
  )

  ## Matched pair design, 5 blocks of 2
  dat$blocks <- rep(1:5, each = 2)
  declaration <- randomizr::declare_ra(
    N = nrow(dat),
    blocks = dat$blocks,
    block_m = rep(1, 5)
  )
  treatment_perms <- randomizr::obtain_permutation_matrix(declaration)

  ests <- apply(
    treatment_perms,
    2,
    function(x) {
      dat$Z <- x
      dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
      dim <- difference_in_means(
        Y ~ Z,
        blocks = blocks,
        data = dat
      )
      dim$coefficients
    }
  )

  expect_equivalent(
    trueSATE,
    mean(ests)
  )

  ## block randomized design, 2 blocks of 5
  dat$blocks <- rep(1:2, each = 5)
  declaration <- randomizr::declare_ra(
    N = nrow(dat),
    blocks = dat$blocks,
    block_m = c(3, 3)
  )
  treatment_perms <- randomizr::obtain_permutation_matrix(declaration)

  ests <- apply(
    treatment_perms,
    2,
    function(x) {
      dat$Z <- x
      dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
      dim <- difference_in_means(
        Y ~ Z,
        blocks = blocks,
        data = dat
      )
      dim$coefficients
    }
  )

  expect_equivalent(
    trueSATE,
    mean(ests)
  )

  ## cluster matched pair, different sized blocks
  dat$blocks <- rep(1:3, times = c(4, 4, 2))
  dat$clusters <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 6)
  declaration <- randomizr::declare_ra(
    N = nrow(dat),
    blocks = dat$blocks,
    clusters = dat$clusters
  )
  treatment_perms <- randomizr::obtain_permutation_matrix(declaration)

  ests <- apply(
    treatment_perms,
    2,
    function(x) {
      dat$Z <- x
      dat$Y <- ifelse(dat$Z, dat$Y1, dat$Y0)
      dim <- difference_in_means(
        Y ~ Z,
        blocks = blocks,
        clusters = clusters,
        data = dat
      )
      dim$coefficients
    }
  )

  expect_equivalent(
    trueSATE,
    mean(ests)
  )
})

test_that("DIM matches lm_robust under certain conditions", {
  n <- 400
  dat <- data.frame(Y = rnorm(n))
  ## DIM and lm_robust agree without clustering except for DoF because DIM uses Satterthwaite approx
  dat$z <- c(0, 1)
  lm_o <- lm_robust(Y ~ z, data = dat)
  dim_o <- difference_in_means(Y ~ z, data = dat)
  expect_equivalent(
    tidy(lm_o)[2, 1:3],
    tidy(dim_o)[, 1:3]
  )

  ## DIM and lm_robust agree with clustering becuase DIM just uses lm_robust w/ CR2!
  dat$cl_diff_size <- sample(100, size = 400, replace = TRUE)
  dat$z_clustered <- as.numeric(dat$cl_diff_size <= 50)

  lm_cl_o <- lm_robust(Y ~ z_clustered, clusters = cl_diff_size, data = dat)
  dim_cl_o <- difference_in_means(Y ~ z_clustered, clusters = cl_diff_size, data = dat)
  expect_equivalent(
    tidy(lm_cl_o)[2, ],
    tidy(dim_cl_o)
  )

  ## Blocked design is equivalent to lm_lin
  dat$bl <- rep(1:25, each = 16)
  dat$z_blocked <- rep(c(0, 1), each = 2)

  lm_bl_o <- lm_lin(Y ~ z_blocked, ~ factor(bl), data = dat)
  dim_bl_o <- difference_in_means(Y ~ z_blocked, blocks = bl, data = dat)

  # Not identical since row name of lm_bl_o is 2 due to the intercept
  expect_equivalent(
    tidy(lm_bl_o)[2, ],
    tidy(dim_bl_o)
  )

  ## Block-clustered is equivalent to lm_lin
  ## (and indeed uses lm_robust machinery for the ests and ses, the DoF are equivalent with equal clusters
  ## by design
  dat$cl <- rep(1:200, each = 2)
  lm_blcl_o <- lm_lin(Y ~ z_blocked, ~ factor(bl), data = dat, clusters = cl)
  dim_blcl_o <- difference_in_means(Y ~ z_blocked, blocks = bl, clusters = cl, data = dat)

  # Not identical since row name of lm_bl_o is 2 due to the intercept
  expect_equivalent(
    tidy(lm_blcl_o)[2, ],
    tidy(dim_blcl_o)
  )

  # With weights now, identical to lm_robust, HC2 by force! except for matched pairs which fails
  dat$w <- runif(nrow(dat))



  # simple W
  expect_equivalent(
    tidy(lm_robust(Y ~ z, data = dat, weights = w))[2, ],
    tidy(difference_in_means(Y ~ z, data = dat, weights = w))
  )

  # blocked W
  expect_equivalent(
    tidy(lm_lin(Y ~ z_blocked, ~ factor(bl), weights = w, data = dat))[2, ],
    tidy(difference_in_means(Y ~ z_blocked, data = dat, blocks = bl, weights = w))
  )

  # blocked-clustered W (goes to CR2)
  # DF different in clustered case
  dim_bl_cl_w <- difference_in_means(
    Y ~ z_blocked,
    data = dat,
    clusters = cl,
    blocks = bl,
    weights = w
  )

  expect_equivalent(
    tidy(lm_lin(Y ~ z_blocked, ~ factor(bl), clusters = cl, weights = w, data = dat))[2, 1:3],
    tidy(dim_bl_cl_w)[, 1:3]
  )

  expect_equal(
    dim_bl_cl_w$design,
    "Block-clustered (weighted)"
  )

  # Clustered W
  dim_cl_w <- difference_in_means(
    Y ~ z_clustered,
    data = dat,
    clusters = cl_diff_size,
    weights = w
  )

  expect_equivalent(
    tidy(lm_robust(Y ~ z_clustered, clusters = cl_diff_size, weights = w, data = dat))[2, 1:3],
    tidy(dim_cl_w)[, 1:3]
  )
  expect_equal(
    dim_cl_w$design,
    "Clustered (weighted)"
  )

  # errors with matched pairs
  dat$mps <- rep(1:(n / 2), each = 2)
  dat$z_mps <- rep(0:1, times = (n / 2))
  expect_error(
    difference_in_means(Y ~ z_mps, data = dat, weights = w, blocks = mps),
    "Cannot use `weights` with matched pairs design at the moment"
  )
})
