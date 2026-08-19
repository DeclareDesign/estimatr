# Data used by the tests that compare this package against estimatr 1.0.6.
#
# These live in one place because two things consume them: the test files, and
# `data-raw/make_estimatr_reference.R`, which runs the same fits under estimatr
# 1.0.6 and records the answers. If the two ever built their data separately
# the recorded reference would drift from what the tests compare against, and
# nothing would report it. Each builder seeds its own RNG, so it does not
# matter what order they are called in or what has drawn before them.

ref_data_vs <- function() {
  set.seed(42)
  n <- 200
  dat <- data.frame(
    y  = rnorm(n),
    y2 = rnorm(n),
    x  = rnorm(n),
    z  = rbinom(n, 1, 0.5),
    cl = rep(1:20, 10),
    bl = rep(1:20, each = 10),
    w  = runif(n, 0.5, 2.0),
    iv = rnorm(n) + rnorm(n, 0, 0.3)
  )
  dat$z_block <- rep(rep(c(0L, 1L), 5L), 20L)
  dat$cl_z <- as.integer(dat$cl %% 2 == 0)
  dat
}

ref_data_pairs <- function() {
  set.seed(99)
  n_pairs <- 50
  data.frame(
    y  = rnorm(n_pairs * 2),
    z  = rep(c(0L, 1L), n_pairs),
    pr = rep(seq_len(n_pairs), each = 2)
  )
}

ref_data_fe <- function() {
  set.seed(42)
  n <- 200
  data.frame(
    y  = rnorm(n),
    y2 = rnorm(n),
    x  = rnorm(n),
    z  = rbinom(n, 1, 0.5),
    bl = rep(1:20, each = 10),
    cl = rep(1:10, 20),
    iv = rnorm(n) + rnorm(n, 0, 0.3)
  )
}

ref_data_fe_weighted <- function() {
  dat <- ref_data_fe()
  set.seed(99)
  dat$w <- runif(nrow(dat), 0.5, 2)
  dat
}

# The three leverage configurations are drawn in sequence from one seed, so
# they must be built together rather than one at a time.
ref_data_leverage <- function() {
  set.seed(7)
  mk <- function(n, G, unbalanced = FALSE) {
    g <- if (unbalanced) {
      factor(sample(seq_len(G), n, replace = TRUE))
    } else {
      factor(rep(seq_len(G), length.out = n))
    }
    data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n), g = g,
               wts = runif(n, 0.5, 3))
  }
  list(
    balanced   = mk(200, 20),
    unbalanced = mk(200, 17, unbalanced = TRUE),
    many_small = mk(300, 60)
  )
}

ref_data_ht <- function() {
  set.seed(42)
  data.frame(y = rnorm(40))
}

ref_data_post <- function() {
  set.seed(1)
  data.frame(y = rnorm(40), z = rep(0:1, 20), x = rnorm(40))
}

ref_data_lin <- function() {
  set.seed(1)
  data.frame(
    y  = rnorm(90),
    zf = factor(rep(c("a", "b", "c"), length.out = 90)),
    zb = rep(0:1, length.out = 90),
    zn = rep(c(0, 2, 5), length.out = 90),
    x  = rnorm(90),
    w  = rnorm(90)
  )
}

# The lm_lin predict cases, shared between the tests and the recorder.
ref_lin_cases <- function() {
  list(
    list(lbl = "factor treatment",        z = "zf", cov = ~ x),
    list(lbl = "factor, two covariates",  z = "zf", cov = ~ x + w),
    list(lbl = "binary treatment",        z = "zb", cov = ~ x),
    list(lbl = "binary, two covariates",  z = "zb", cov = ~ x + w),
    list(lbl = "numeric multi-valued",    z = "zn", cov = ~ x),
    list(lbl = "numeric multi, two covs", z = "zn", cov = ~ x + w)
  )
}
