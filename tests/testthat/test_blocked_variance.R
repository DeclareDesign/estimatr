library(estimatrZero)

# Pashley and Miratrix (2021) variance estimation for blocked designs with
# blocks of different sizes, including blocks holding a singleton treated or
# control unit. estimatr errors on those designs or falls back to the
# matched-pairs estimator for every block.

make_blocked <- function(block_sizes, m_each, seed = 1, tau = 0.5) {
  set.seed(seed)
  K <- length(block_sizes)
  bl <- rep(seq_len(K), block_sizes)
  n <- sum(block_sizes)
  y0 <- rnorm(n) + rep(rnorm(K, 0, 2), block_sizes)
  z <- randomizr::block_ra(blocks = bl, block_m = m_each)
  data.frame(y = y0 + tau * z, z = z, bl = bl)
}

# ---- the designs that used to be unusable ----

test_that("a block with a singleton arm is estimable", {
  skip_if_not_installed("randomizr")
  # 3-unit blocks with 1 treated and 2 control: estimatr errors with
  # "every block must have at least two treated and control units"
  d <- make_blocked(c(8, 8, 3, 3, 3, 3), c(4, 4, 1, 1, 1, 1))
  m <- difference_in_means(y ~ z, data = d, blocks = bl)
  expect_equal(m$design, "Hybrid blocked")
  expect_true(is.finite(m$std.error[[1]]))
  expect_true(is.finite(m$df[[1]]))
})

test_that("mixed block sizes no longer collapse to matched pairs", {
  skip_if_not_installed("randomizr")
  # estimatr warns "Using matched pairs variance estimator" and applies it to
  # every block, including the big ones
  d <- make_blocked(c(10, 10, 2, 2, 2, 2), c(5, 5, 1, 1, 1, 1))
  expect_no_warning(m <- difference_in_means(y ~ z, data = d, blocks = bl))
  expect_equal(m$design, "Hybrid blocked")
})

test_that("small blocks of varying size are estimable", {
  skip_if_not_installed("randomizr")
  d <- make_blocked(c(3, 4, 5, 6, 7, 8), rep(1, 6))
  m <- difference_in_means(y ~ z, data = d, blocks = bl)
  expect_equal(m$design, "Small blocks")
  expect_true(is.finite(m$std.error[[1]]))
})

# ---- the reference implementation ----

test_that("standard errors match blkvar, the authors' own package", {
  skip_if_not_installed("randomizr")
  skip_if_not_installed("blkvar")
  skip_if_not_installed("dplyr")
  suppressMessages(library(dplyr))    # blkvar calls dplyr::n() without importing it

  designs <- list(
    list(rep(6, 4), rep(3, 4)),
    list(c(4, 6, 8, 10), c(2, 3, 4, 5)),
    list(c(3, 4, 5, 6, 7, 8), rep(1, 6)),
    list(c(8, 8, 3, 3, 3, 3), c(4, 4, 1, 1, 1, 1)),
    list(c(9, 9, 4, 5, 6), c(4, 4, 3, 4, 5)),
    list(rep(2, 10), rep(1, 10))
  )
  for (i in seq_along(designs)) {
    d <- make_blocked(designs[[i]][[1]], designs[[i]][[2]], seed = i)
    ours <- difference_in_means(y ~ z, data = d, blocks = bl)
    theirs <- blkvar::block_estimator(Yobs = d$y, Z = d$z, B = factor(d$bl),
                                      method = "hybrid_p", throw.warnings = FALSE)
    expect_equal(ours$coefficients[[1]], theirs$ATE_hat, tolerance = 1e-10)
    expect_equal(ours$std.error[[1]], theirs$se_est, tolerance = 1e-10)
  }
})

# ---- the two special cases still reduce to what they always were ----

test_that("all-big designs keep the plain blocked estimator", {
  skip_if_not_installed("randomizr")
  d <- make_blocked(rep(6, 5), rep(3, 5))
  m <- difference_in_means(y ~ z, data = d, blocks = bl)
  expect_equal(m$design, "Blocked")

  # equation 4 by hand
  parts <- lapply(split(d, d$bl), function(x) {
    y2 <- x$y[x$z == 1]; y1 <- x$y[x$z == 0]
    c(n = nrow(x), v = var(y2) / length(y2) + var(y1) / length(y1))
  })
  parts <- do.call(rbind, parts)
  expect_equal(m$std.error[[1]],
               sqrt(sum(parts[, "n"]^2 * parts[, "v"]) / sum(parts[, "n"])^2),
               tolerance = 1e-12)
  expect_equal(m$df[[1]], nrow(d) - 2 * 5)
})

test_that("matched pairs keep the matched-pairs estimator and df", {
  skip_if_not_installed("randomizr")
  d <- make_blocked(rep(2, 12), rep(1, 12))
  m <- difference_in_means(y ~ z, data = d, blocks = bl)
  expect_equal(m$design, "Matched-pair")

  # equation 5 by hand
  tau_k <- vapply(split(d, d$bl), function(x) x$y[x$z == 1] - x$y[x$z == 0], numeric(1))
  expect_equal(m$std.error[[1]],
               sqrt(sum((tau_k - mean(tau_k))^2) / (12 * 11)), tolerance = 1e-12)
  expect_equal(m$df[[1]], 11)
})

# ---- what is refused, and why ----

test_that("a single small block is refused rather than given zero variance", {
  skip_if_not_installed("randomizr")
  # one small block leaves nothing to compare it against, and the estimator
  # would silently return a zero contribution
  d <- make_blocked(c(8, 8, 8, 3), c(4, 4, 4, 1))
  expect_error(difference_in_means(y ~ z, data = d, blocks = bl),
               "Only one block")
})

test_that("a small block holding half the small-block sample is refused", {
  skip_if_not_installed("randomizr")
  # equation 8 needs every small block under half of n_sb, or the weights go
  # negative. Equal-size small blocks avoid it by taking the matched-pairs
  # branch, so the boundary only bites when the sizes vary.
  d <- make_blocked(c(10, 3, 5), c(5, 1, 1))
  expect_error(difference_in_means(y ~ z, data = d, blocks = bl),
               "half or more")
})

test_that("blocks with no variation in treatment still error", {
  d <- data.frame(y = rnorm(8), z = c(1, 1, 0, 0, 1, 1, 1, 1),
                  bl = rep(1:2, each = 4))
  expect_error(difference_in_means(y ~ z, data = d, blocks = bl),
               "both treatment conditions")
})

# ---- the hybrid is conservative, and covers ----

test_that("hybrid variance is conservative and covers at the nominal rate", {
  skip_if_not_installed("randomizr")
  set.seed(20)
  block_sizes <- c(8, 8, 3, 3, 3, 3)
  m_each <- c(4, 4, 1, 1, 1, 1)
  K <- length(block_sizes)
  bl <- rep(seq_len(K), block_sizes)
  n <- sum(block_sizes)
  y0 <- rnorm(n) + rep(rnorm(K, 0, 1), block_sizes)
  y1 <- y0 + 0.5 + rep(rnorm(K, 0, 0.8), block_sizes)   # effects vary by block
  sate <- mean(y1 - y0)

  # potential outcomes fixed, assignment is the only randomness
  out <- replicate(1500, {
    z <- randomizr::block_ra(blocks = bl, block_m = m_each)
    fit <- difference_in_means(y ~ z, blocks = bl,
      data = data.frame(y = ifelse(z == 1, y1, y0), z = z, bl = bl))
    c(fit$coefficients[[1]], fit$std.error[[1]]^2,
      fit$conf.low[[1]], fit$conf.high[[1]])
  })

  expect_lt(abs(mean(out[1, ]) - sate), 0.05)                   # unbiased
  expect_gt(mean(out[2, ]) / var(out[1, ]), 1)                  # conservative
  expect_gt(mean(out[3, ] <= sate & out[4, ] >= sate), 0.94)    # covers
})

# ---- cluster-randomized blocks (estimatr #336 is unit-level only) ----

# Pashley and Miratrix cover treatment assigned WITHIN blocks, not blocks of
# clusters ("we are focusing on treatment assigned within cluster", section 2),
# so the hybrid above does not extend here. What is enforced instead is that a
# block must be able to estimate its own variance.

make_clustered <- function(cl_per_block, m_cl, cl_size = 4, seed = 1) {
  set.seed(seed)
  K <- length(cl_per_block)
  bl_of_cl <- rep(seq_len(K), cl_per_block)
  ncl <- length(bl_of_cl)
  cl <- rep(seq_len(ncl), each = cl_size)
  bl <- rep(bl_of_cl, each = cl_size)
  z <- rep(randomizr::block_ra(blocks = bl_of_cl, block_m = m_cl), each = cl_size)
  data.frame(y = rnorm(length(z)) + rep(rnorm(ncl, 0, 1.5), each = cl_size) + 0.5 * z,
             z = z, cl = cl, bl = bl)
}

test_that("a block with one treated cluster is refused, not silently wrong", {
  skip_if_not_installed("randomizr")
  d <- make_clustered(c(6, 6, 6, 4), c(3, 3, 3, 1))
  expect_error(
    difference_in_means(y ~ z, data = d, blocks = bl, clusters = cl),
    "only one treated or one control `cluster`"
  )
  # the message names the offending block
  expect_error(
    difference_in_means(y ~ z, data = d, blocks = bl, clusters = cl), "4"
  )
})

test_that("matched-pair clustered designs are exempt from that check", {
  skip_if_not_installed("randomizr")
  # every block here has one treated and one control cluster by construction,
  # and the variance comes from across blocks, so it is estimable
  d <- make_clustered(c(2, 2, 2, 2), rep(1, 4))
  m <- difference_in_means(y ~ z, data = d, blocks = bl, clusters = cl)
  expect_equal(m$design, "Matched-pair clustered")
  expect_true(is.finite(m$std.error[[1]]))
})

test_that("clustered blocks with two or more clusters per arm still work", {
  skip_if_not_installed("randomizr")
  d <- make_clustered(c(6, 6), c(3, 3))
  m <- difference_in_means(y ~ z, data = d, blocks = bl, clusters = cl)
  expect_equal(m$design, "Block-clustered")
  expect_true(is.finite(m$std.error[[1]]))
})

test_that("within-block cluster variance is exact when both arms have 2+ clusters", {
  skip_if_not_installed("randomizr")
  # Enumerate every assignment of 6 clusters with 3 treated: with the
  # cluster-level potential outcomes fixed, the mean estimated variance must
  # equal the true variance of the estimator.
  set.seed(77)
  cl_size <- 4; ncl <- 6
  y0 <- rnorm(ncl * cl_size) + rep(rnorm(ncl, 0, 1.5), each = cl_size)
  y1 <- y0 + 0.5
  cl <- rep(seq_len(ncl), each = cl_size)
  combos <- utils::combn(ncl, 3)

  res <- apply(combos, 2, function(tr) {
    z <- as.integer(cl %in% tr)
    f <- difference_in_means(y ~ z, clusters = cl,
      data = data.frame(y = ifelse(z == 1, y1, y0), z = z, cl = cl))
    c(f$coefficients[[1]], f$std.error[[1]]^2)
  })
  true_var <- var(res[1, ]) * (ncol(combos) - 1) / ncol(combos)
  expect_equal(mean(res[2, ]), true_var, tolerance = 1e-8)
})
