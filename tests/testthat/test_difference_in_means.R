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

test_that("C13: se_type = \"none\" on a blocked design withholds the variance only", {
  # The Pashley-Miratrix branch has its own `se_type == "none"` arm, separate
  # from the one every other design takes, and it had never been run: the
  # estimate must be the one the default fit gives, with the variance and the
  # df withheld rather than a different point estimate.
  set.seed(343)
  d <- data.frame(y = rnorm(200),
                  z = rep(rep(0:1, 5), 20),
                  bl = rep(1:20, each = 10))

  none <- difference_in_means(y ~ z, blocks = bl, data = d, se_type = "none")
  full <- difference_in_means(y ~ z, blocks = bl, data = d)

  expect_equal(none$design, "Blocked")
  expect_equal(none$coefficients, full$coefficients)
  expect_true(is.na(none$std.error))
  expect_true(is.na(none$df))
  expect_true(is.na(none$p.value))
  expect_false(is.na(full$std.error))
})

test_that("C13: blocks of mixed cluster count say so and use the matched-pair estimator", {
  # Two clusters in one block and four in another: the design is neither
  # matched pairs nor a block design with estimable within-block variance, and
  # the matched-pair estimator is used across blocks with a warning that says
  # which assumption was made.
  set.seed(343)
  cl <- rep(1:6, each = 3)
  bl <- c(rep(1, 6), rep(2, 12))
  d <- data.frame(y = rnorm(18), cl = cl, bl = bl,
                  z = as.integer(cl %in% c(2, 4, 6)))

  expect_warning(
    m <- difference_in_means(y ~ z, blocks = bl, clusters = cl, data = d),
    "two units/`clusters` while other blocks have more"
  )
  expect_equal(m$design, "Matched-pair clustered")

  # Two blocks, so the across-block estimator has one degree of freedom.
  expect_equal(unname(m$df), 1)
})

test_that("C13: a block with a single cluster is refused on the clustered path", {
  # check_clusters_blocks() counts clusters rather than units here, so a block
  # of six units in one cluster is a block of one.
  set.seed(343)
  cl <- c(rep(1L, 6), rep(3:6, each = 3))
  bl <- c(rep(1, 6), rep(2, 12))
  d <- data.frame(y = rnorm(18), cl = cl, bl = bl,
                  z = as.integer(cl %in% c(1, 4, 6)))

  expect_error(
    difference_in_means(y ~ z, blocks = bl, clusters = cl, data = d),
    "All `blocks` must have multiple units"
  )
})

test_that("C13: a block with only one treatment condition is refused in its own words", {
  # difference_in_means_internal() is called once per block as well as on the
  # whole sample, and this is the per-block refusal: the block has clusters in
  # both arms by count but every unit in it is treated.
  set.seed(343)
  cl <- rep(1:6, each = 3)
  bl <- c(rep(1, 6), rep(2, 12))
  d <- data.frame(y = rnorm(18), cl = cl, bl = bl,
                  z = as.integer(cl %in% c(2, 4, 6)))
  d$z[d$bl == 1] <- 1L

  expect_error(
    suppressWarnings(
      difference_in_means(y ~ z, blocks = bl, clusters = cl, data = d)
    ),
    "Must have units with both treatment conditions"
  )
})

test_that("C13: a factor treatment takes its conditions from the levels", {
  # parse_conditions() reads levels(droplevels()) rather than sort(unique()) for
  # a factor, so an unused level cannot become a condition and the contrast
  # follows the level order rather than the alphabet.
  set.seed(343)
  d <- data.frame(y = rnorm(60),
                  z = factor(rep(c("trt", "ctl"), 30), levels = c("trt", "ctl")))
  d$y <- d$y + (d$z == "trt")

  m <- difference_in_means(y ~ z, data = d)

  # levels() puts trt first, so trt is condition1 and the contrast is ctl - trt.
  expect_equal(m$term, "zctl")
  expect_lt(m$coefficients[[1]], 0)

  # An unused level is dropped rather than demanded as a third condition.
  d$z <- factor(as.character(d$z), levels = c("trt", "ctl", "never"))
  expect_equal(difference_in_means(y ~ z, data = d)$coefficients, m$coefficients)
})

test_that("C13: naming one condition infers the other from the two present", {
  # Two separate branches in parse_conditions, neither run: given condition2
  # alone the other value becomes condition1, and given condition1 alone the
  # other becomes condition2. Naming either one alone is the same contrast as
  # naming neither, up to the sign the choice implies.
  set.seed(343)
  d <- data.frame(y = rnorm(60), z = factor(rep(c("ctl", "trt"), 30)))

  both <- difference_in_means(y ~ z, data = d)
  c2_only <- difference_in_means(y ~ z, data = d, condition2 = "ctl")
  c1_only <- difference_in_means(y ~ z, data = d, condition1 = "trt")

  # condition2 = "ctl" leaves "trt" as condition1, so the contrast reverses,
  # and the coefficient is named for the condition2 it now carries.
  expect_equal(c2_only$term, "zctl")
  expect_equal(unname(c2_only$coefficients), unname(-both$coefficients))
  # condition1 = "trt" leaves "ctl" as condition2, the same reversal.
  expect_equal(c1_only$term, "zctl")
  expect_equal(unname(c1_only$coefficients), unname(-both$coefficients))
  expect_equal(unname(c1_only$std.error), unname(both$std.error))
  expect_equal(unname(c2_only$std.error), unname(both$std.error))
})

test_that("C13: weights with matched pairs are refused, which is what makes one branch unreachable", {
  # `difference_in_means()` has a "Matched-pair" arm at the pair-matched
  # blocked path (the unit-randomized across-pair variance, Gerber & Green
  # 2012 p. 77 eq. 3.16). Reaching it needs blocks, no clusters and weights,
  # because the unweighted unclustered case is taken by the Pashley-Miratrix
  # branch above it; and weights with matched pairs are refused here, one
  # call earlier. The arm is therefore dead while this refusal stands, and
  # this test is what will report it if the refusal is ever lifted.
  set.seed(343)
  np <- 40
  d <- data.frame(y = rnorm(np * 2), z = rep(0:1, np),
                  pr = rep(seq_len(np), each = 2), w = runif(np * 2, 0.5, 2))

  expect_error(
    difference_in_means(y ~ z, blocks = pr, weights = w, data = d),
    "Cannot use `weights` with matched pairs design"
  )
})
