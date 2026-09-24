library(estimatr)

# difference_in_means's own behaviour: which design it recognises, how it reads
# condition1 and condition2, and what ci = FALSE withholds.
#
# The estimate and variance of each design against the regression it
# describes are in test_equivalence.R; blocked designs with small blocks, and
# blocks of clusters, in test_blocked_variance.R; refusals in test_errors.R.

balanced_in_blocks <- function(blocks) {
  as.integer(unlist(lapply(split(seq_along(blocks), blocks), function(idx) {
    z <- integer(length(idx))
    z[sample(seq_along(z), length(idx) %/% 2L)] <- 1L
    z
  })))
}

test_that("each design is recognised and named", {
  set.seed(42)
  n <- 200
  d <- data.frame(y = rnorm(n), z = rbinom(n, 1, 0.5), w = runif(n),
                  cl = rep(1:20, 10), block = rep(1:20, each = 10))
  d$z_block <- rep(rep(c(0L, 1L), 5L), 20L)
  d$z_cl <- as.integer(d$cl %% 2 == 0)
  d$pair <- rep(1:100, each = 2)
  d$z_pair <- rep(0:1, 100)

  expect_equal(difference_in_means(y ~ z, data = d)$design, "Standard")
  blocked <- difference_in_means(y ~ z_block, blocks = block, data = d)
  expect_equal(blocked$design, "Blocked")
  expect_equal(unname(blocked$df), n - 2 * 20)
  expect_equal(difference_in_means(y ~ z_cl, clusters = cl, data = d)$design, "Clustered")
  expect_equal(difference_in_means(y ~ z_pair, blocks = pair, data = d)$design, "Matched-pair")
  expect_equal(difference_in_means(y ~ z, weights = w, data = d)$design, "Standard (weighted)")
  expect_equal(difference_in_means(y ~ z_block, weights = w, blocks = block, data = d)$design,
               "Blocked (weighted)")

  # Two clusters in each block, one of them treated.
  set.seed(42)
  bc <- data.frame(Y = rnorm(80), cl = rep(1:8, 10), bl = rep(rep(1:4, each = 2), 10))
  bc$Z <- as.integer(bc$cl %in% c(1, 3, 5, 7))
  pair_clustered <- difference_in_means(Y ~ Z, clusters = cl, blocks = bl, data = bc)
  expect_equal(pair_clustered$design, "Matched-pair clustered")
  expect_true(is.finite(pair_clustered$std.error[["Z"]]))
})

test_that("with more than two conditions, condition1 and condition2 pick the pair", {
  set.seed(3)
  d <- data.frame(Y = rnorm(100), Z = sample(1:3, 100, replace = TRUE))
  m12 <- difference_in_means(Y ~ Z, condition1 = 1, condition2 = 2, data = d)
  welch <- t.test(d$Y[d$Z == 2], d$Y[d$Z == 1])
  expect_equal(unname(m12$coefficients), mean(d$Y[d$Z == 2]) - mean(d$Y[d$Z == 1]),
               tolerance = 1e-12)
  expect_equal(unname(m12$std.error), welch$stderr, tolerance = 1e-12)
  expect_equal(unname(m12$df), unname(welch$parameter), tolerance = 1e-12)
  expect_equal(m12$design, "Standard")
})

test_that("reversing the conditions negates the estimate and leaves the standard error", {
  set.seed(42)
  d <- data.frame(Y = rnorm(80), Z = rbinom(80, 1, 0.5))
  forward <- difference_in_means(Y ~ Z, condition1 = 0, condition2 = 1, data = d)
  reverse <- difference_in_means(Y ~ Z, condition1 = 1, condition2 = 0, data = d)
  expect_equal(forward$coefficients[[1]], -reverse$coefficients[[1]])
  expect_equal(forward$std.error[[1]], reverse$std.error[[1]])

  set.seed(42)
  b <- data.frame(Y = rnorm(60), Z = rbinom(60, 1, 0.5),
                  bl = sample(c("A", "B", "C"), 60, replace = TRUE))
  forward <- difference_in_means(Y ~ Z, blocks = bl, condition1 = 0, condition2 = 1, data = b)
  reverse <- difference_in_means(Y ~ Z, blocks = bl, condition1 = 1, condition2 = 0, data = b)
  expect_equal(forward$coefficients[[1]], -reverse$coefficients[[1]])
  expect_equal(forward$std.error[[1]], reverse$std.error[[1]])
  expect_equal(forward$design, "Blocked")
})

test_that("ci = FALSE withholds the p-value and the interval", {
  set.seed(1)
  d <- data.frame(Y = rnorm(40), Z = rbinom(40, 1, 0.5))
  full <- difference_in_means(Y ~ Z, data = d)
  m <- difference_in_means(Y ~ Z, data = d, ci = FALSE)
  expect_equal(m$std.error, full$std.error)
  expect_true(is.na(m$p.value[[1]]))
  expect_true(is.na(m$conf.low[[1]]))
  expect_true(is.na(m$conf.high[[1]]))
})

test_that("a cbind() outcome is refused in the right words", {
  # The multivariate formula otherwise reached the blocked path and was refused
  # for want of "both treatment conditions within each block", which describes
  # a different design problem and sends the reader to look at their blocks.
  set.seed(343)
  d <- data.frame(y = rnorm(60), y2 = rnorm(60), z = rbinom(60, 1, 0.5))

  expect_error(
    difference_in_means(cbind(y, y2) ~ z, data = d),
    "does not support multiple outcomes"
  )
})

test_that("update() refits a difference_in_means", {
  # The fit stored no call, so update() stopped with "need an object with call
  # component" where horvitz_thompson and the regression fits all refit.
  set.seed(343)
  d <- data.frame(y = rnorm(60), y2 = rnorm(60), z = rbinom(60, 1, 0.5))

  m <- difference_in_means(y ~ z, data = d)
  expect_false(is.null(m$call))

  same <- update(m, . ~ .)
  expect_equal(same$coefficients, m$coefficients)
  expect_equal(same$std.error, m$std.error)

  other <- update(m, y2 ~ .)
  expect_equal(other$coefficients,
               difference_in_means(y2 ~ z, data = d)$coefficients)
})
