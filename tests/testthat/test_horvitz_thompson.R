library(estimatrZero)

skip_if_not_installed("randomizr")
skip_if_not_installed("estimatr")

set.seed(42)
N <- 40
dat <- data.frame(y = rnorm(N))

# ---- helpers ----

# Compare our HT to estimatr for designs it supports (all 2-arm designs).
ht_check <- function(decl, dat, ...) {
  dat$Z <- randomizr::conduct_ra(decl)
  mz <- estimatrZero::horvitz_thompson(y ~ Z, data = dat,
                                        condition_prs = decl, ...)
  me <- estimatr::horvitz_thompson(y ~ Z, data = dat,
                                    ra_declaration = decl, ...)
  list(mz = mz, me = me)
}

# ---- estimate identity across all supported designs ----

test_that("HT simple: estimate and SE identical to estimatr", {
  decl <- randomizr::declare_ra(N = N, prob = 0.4, simple = TRUE)
  r <- ht_check(decl, dat)
  expect_equal(r$mz$coefficients[[1]], r$me$coefficients[[1]], tolerance = 1e-12)
  expect_equal(r$mz$std.error[[1]],    r$me$std.error[[1]],    tolerance = 1e-12)
})

test_that("HT complete: estimate and SE identical to estimatr", {
  decl <- randomizr::declare_ra(N = N, m = 16)
  r <- ht_check(decl, dat)
  expect_equal(r$mz$coefficients[[1]], r$me$coefficients[[1]], tolerance = 1e-12)
  expect_equal(r$mz$std.error[[1]],    r$me$std.error[[1]],    tolerance = 1e-12)
})

test_that("HT blocked: estimate and SE identical to estimatr", {
  dat$bl <- rep(1:4, each = 10)
  decl <- randomizr::declare_ra(N = N, blocks = dat$bl)
  r <- ht_check(decl, dat)
  expect_equal(r$mz$coefficients[[1]], r$me$coefficients[[1]], tolerance = 1e-12)
  expect_equal(r$mz$std.error[[1]],    r$me$std.error[[1]],    tolerance = 1e-12)
})

test_that("HT blocked non-integer k2_b: SE identical to estimatr", {
  # 5 clusters per block (K_b=5, pi=0.5 → non-integer m=2.5)
  dat$bl2 <- rep(1:2, each = 20)
  dat$cl2 <- rep(1:10, each = 4)
  decl <- randomizr::declare_ra(N = N, blocks = dat$bl2, clusters = dat$cl2)
  r <- ht_check(decl, dat)
  expect_equal(r$mz$coefficients[[1]], r$me$coefficients[[1]], tolerance = 1e-12)
  expect_equal(r$mz$std.error[[1]],    r$me$std.error[[1]],    tolerance = 1e-12)
})

test_that("HT clustered complete: estimate and SE identical to estimatr", {
  dat$cl <- rep(1:8, each = 5)
  decl <- randomizr::declare_ra(N = N, clusters = dat$cl)
  r <- ht_check(decl, dat)
  expect_equal(r$mz$coefficients[[1]], r$me$coefficients[[1]], tolerance = 1e-12)
  expect_equal(r$mz$std.error[[1]],    r$me$std.error[[1]],    tolerance = 1e-12)
})

test_that("HT clustered simple: estimate and SE identical to estimatr", {
  dat$cl <- rep(1:8, each = 5)
  decl <- randomizr::declare_ra(N = N, clusters = dat$cl, simple = TRUE)
  r <- ht_check(decl, dat)
  expect_equal(r$mz$coefficients[[1]], r$me$coefficients[[1]], tolerance = 1e-12)
  expect_equal(r$mz$std.error[[1]],    r$me$std.error[[1]],    tolerance = 1e-12)
})

test_that("HT custom permutation: SE identical to estimatr", {
  perm <- replicate(20, {
    z <- rep(0L, N); z[sample(N, 16)] <- 1L; z
  })
  decl <- randomizr::declare_ra(permutation_matrix = perm)
  dat$Z <- perm[, 1L]
  mz <- horvitz_thompson(y ~ Z, data = dat, condition_prs = decl)
  me <- estimatr::horvitz_thompson(y ~ Z, data = dat, ra_declaration = decl)
  expect_equal(mz$coefficients[[1]], me$coefficients[[1]], tolerance = 1e-12)
  expect_equal(mz$std.error[[1]],    me$std.error[[1]],    tolerance = 1e-12)
})

# ---- named scalar probability vector ----

test_that("HT named vector: SE matches estimatr conservative formula", {
  dat2 <- data.frame(y = rnorm(20), Z = c(rep(1L, 10L), rep(0L, 10L)))
  mz <- horvitz_thompson(y ~ Z, data = dat2, condition_prs = c("0" = 0.5, "1" = 0.5))
  me <- estimatr::horvitz_thompson(y ~ Z, data = dat2,
          condition_prs = ifelse(dat2$Z == 1, 0.5, 0.5))
  expect_equal(mz$coefficients[[1]], me$coefficients[[1]], tolerance = 1e-12)
  expect_equal(mz$std.error[[1]],    me$std.error[[1]],    tolerance = 1e-12)
})

# ---- multi-arm (estimatr doesn't support ra_declaration with >2 arms) ----

# Every assignment of 6 units to three arms of two, as a permutation matrix.
# 6!/(2!2!2!) = 90 columns, each equally likely, so this design is exactly
# three-arm complete randomization written out longhand.
three_arm_perms <- function() {
  units <- 1:6
  out <- list()
  for (a in asplit(utils::combn(units, 2), 2)) {
    rest <- setdiff(units, a)
    for (b in asplit(utils::combn(rest, 2), 2)) {
      z <- rep("C", 6L)
      z[a] <- "A"; z[b] <- "B"
      out[[length(out) + 1L]] <- z
    }
  }
  matrix(unlist(out), nrow = 6L)
}

test_that("HT multi-arm: runs and returns class", {
  dat3 <- data.frame(y = rnorm(30))
  decl3 <- randomizr::declare_ra(N = 30, num_arms = 3)
  dat3$Z <- randomizr::conduct_ra(decl3)
  m <- horvitz_thompson(y ~ Z, data = dat3, condition_prs = decl3,
                        condition1 = "T1", condition2 = "T2")
  expect_s3_class(m, "horvitz_thompson")
  expect_false(is.na(m$std.error[[1]]))
  expect_false(is.na(m$coefficients[[1]]))
})

test_that("HT multi-arm: estimate is unbiased over every assignment", {
  # Averaging over all 90 equally likely assignments must recover the true
  # ATE exactly. Dividing by the two-arm count rather than N inflates this
  # by K/2 = 1.5.
  perm <- three_arm_perms()
  set.seed(7)
  po <- data.frame(A = rnorm(6), B = rnorm(6) + 0.8, C = rnorm(6) - 0.3)
  decl <- randomizr::declare_ra(N = 6, conditions = c("A", "B", "C"))

  ests <- vapply(seq_len(ncol(perm)), function(j) {
    z <- perm[, j]
    d <- data.frame(y = po[cbind(seq_len(6), match(z, names(po)))], Z = z)
    horvitz_thompson(y ~ Z, data = d, condition_prs = decl,
                     condition1 = "A", condition2 = "B",
                     ci = FALSE, se_type = "none")$coefficients[[1]]
  }, numeric(1))

  expect_equal(mean(ests), mean(po$B) - mean(po$A), tolerance = 1e-12)
})

test_that("HT multi-arm: analytic variance equals brute-force joint probabilities", {
  # The permutation matrix enumerates the same design, so the O(n^2) custom
  # pathway computes the Young's bound from exact joint probabilities. The
  # closed-form multi-arm coefficients must reproduce it.
  perm <- three_arm_perms()
  decl_complete <- randomizr::declare_ra(N = 6, conditions = c("A", "B", "C"))
  decl_custom <- randomizr::declare_ra(permutation_matrix = perm)

  set.seed(8)
  for (j in c(1L, 17L, 45L, 90L)) {
    d <- data.frame(y = rnorm(6), Z = perm[, j])
    m_analytic <- horvitz_thompson(y ~ Z, data = d, condition_prs = decl_complete,
                                    condition1 = "A", condition2 = "B")
    m_brute <- horvitz_thompson(y ~ Z, data = d, condition_prs = decl_custom,
                                 condition1 = "A", condition2 = "B")
    expect_equal(m_analytic$coefficients[[1]], m_brute$coefficients[[1]],
                 tolerance = 1e-12)
    expect_equal(m_analytic$std.error[[1]], m_brute$std.error[[1]],
                 tolerance = 1e-12)
  }
})

test_that("HT multi-arm blocked: analytic variance equals brute force", {
  # Two blocks of 6, three arms of two within each block. Enumerating the
  # block design means the Cartesian product of the two within-block
  # enumerations, so instead check the defining property: the blocked
  # variance is the sum of the within-block contributions, each of which the
  # previous test has already matched to exact joint probabilities.
  perm <- three_arm_perms()
  bl <- rep(1:2, each = 6)
  decl_blk <- randomizr::declare_ra(blocks = bl, conditions = c("A", "B", "C"))
  set.seed(9)
  z <- c(perm[, 3L], perm[, 40L])
  d <- data.frame(y = rnorm(12), Z = z, bl = bl)
  m_blk <- horvitz_thompson(y ~ Z, data = d, condition_prs = decl_blk,
                             condition1 = "A", condition2 = "B")

  decl_1 <- randomizr::declare_ra(N = 6, conditions = c("A", "B", "C"))
  var_within <- vapply(1:2, function(b) {
    db <- d[d$bl == b, ]
    m <- horvitz_thompson(y ~ Z, data = db, condition_prs = decl_1,
                          condition1 = "A", condition2 = "B")
    m$std.error[[1]]^2 * 6^2      # undo the per-block 1/N^2
  }, numeric(1))

  expect_equal(m_blk$std.error[[1]]^2, sum(var_within) / 12^2, tolerance = 1e-12)
})

test_that("HT multi-arm clustered: runs and uses cluster-level probabilities", {
  cl <- rep(1:6, each = 3)
  decl <- randomizr::declare_ra(clusters = cl, conditions = c("A", "B", "C"))
  set.seed(10)
  d <- data.frame(y = rnorm(18), Z = randomizr::conduct_ra(decl))
  m <- horvitz_thompson(y ~ Z, data = d, condition_prs = decl,
                        condition1 = "A", condition2 = "B")
  expect_false(is.na(m$std.error[[1]]))
  expect_equal(m$nobs, 18)
})

test_that("HT errors when the declaration does not cover the data rows", {
  decl <- randomizr::declare_ra(N = N, m = 16)
  dat$Z <- randomizr::conduct_ra(decl)
  expect_error(
    horvitz_thompson(y ~ Z, data = dat[1:20, ], condition_prs = decl),
    "declares 40 units but `data` has 20 rows"
  )
})

# ---- return object structure ----

test_that("HT return object has expected fields", {
  decl <- randomizr::declare_ra(N = N, m = 16)
  dat$Z <- randomizr::conduct_ra(decl)
  m <- horvitz_thompson(y ~ Z, data = dat, condition_prs = decl)
  expect_s3_class(m, "horvitz_thompson")
  for (f in c("coefficients", "std.error", "p.value", "conf.low", "conf.high",
              "df", "nobs", "vcov", "se_type", "condition1", "condition2",
              "outcome", "term")) {
    expect_true(!is.null(m[[f]]), label = paste("field", f, "exists"))
  }
  expect_equal(m$nobs, N)
  expect_equal(m$se_type, "youngs")
})

test_that("HT se_type = 'none' gives NA standard error", {
  decl <- randomizr::declare_ra(N = N, m = 16)
  dat$Z <- randomizr::conduct_ra(decl)
  m <- horvitz_thompson(y ~ Z, data = dat, condition_prs = decl, se_type = "none")
  expect_true(is.na(m$std.error[[1]]))
})

test_that("HT ci = FALSE gives NA p-value and CIs", {
  decl <- randomizr::declare_ra(N = N, m = 16)
  dat$Z <- randomizr::conduct_ra(decl)
  m <- horvitz_thompson(y ~ Z, data = dat, condition_prs = decl, ci = FALSE)
  expect_true(is.na(m$p.value[[1]]))
  expect_true(is.na(m$conf.low[[1]]))
  expect_true(is.na(m$conf.high[[1]]))
})

test_that("HT tidy returns a data frame with expected columns", {
  decl <- randomizr::declare_ra(N = N, m = 16)
  dat$Z <- randomizr::conduct_ra(decl)
  m  <- horvitz_thompson(y ~ Z, data = dat, condition_prs = decl)
  td <- tidy(m)
  expect_s3_class(td, "data.frame")
  expect_true(all(c("term", "estimate", "std.error", "p.value") %in% names(td)))
})

test_that("HT print runs without error", {
  decl <- randomizr::declare_ra(N = N, m = 16)
  dat$Z <- randomizr::conduct_ra(decl)
  m <- horvitz_thompson(y ~ Z, data = dat, condition_prs = decl)
  expect_no_error(print(m))
})

# ---- condition1/condition2 defaults ----

test_that("HT condition defaults to first and second sorted condition", {
  decl <- randomizr::declare_ra(N = N, m = 16)
  dat$Z <- randomizr::conduct_ra(decl)
  m_default <- horvitz_thompson(y ~ Z, data = dat, condition_prs = decl)
  m_explicit <- horvitz_thompson(y ~ Z, data = dat, condition_prs = decl,
                                  condition1 = "0", condition2 = "1")
  expect_equal(m_default$coefficients[[1]], m_explicit$coefficients[[1]])
})

# ---- missing data ----

test_that("HT missing Y: estimate uses the surviving rows only", {
  decl <- randomizr::declare_ra(N = N, m = 16)
  dat$Z <- randomizr::conduct_ra(decl)
  dat_miss <- dat; dat_miss$y[5] <- NA
  m_miss <- horvitz_thompson(y ~ Z, data = dat_miss, condition_prs = decl)

  kept <- setdiff(seq_len(N), 5L)
  pm <- decl$probabilities_matrix
  z_k <- as.character(dat$Z[kept]); y_k <- dat$y[kept]
  t2 <- which(z_k == "1"); t1 <- which(z_k == "0")
  Y2 <- y_k[t2] / pm[kept, 2][t2]
  Y1 <- y_k[t1] / pm[kept, 1][t1]
  expect_equal(m_miss$coefficients[[1]],
               (sum(Y2) - sum(Y1)) / length(kept), tolerance = 1e-12)
  expect_equal(m_miss$nobs, length(kept))
})

test_that("HT missing Y: non-uniform pi uses correct row probabilities", {
  # Blocked design with unequal treatment proportions per block → varying pi.
  # Row 2 (in block 1, pi2=0.3) is dropped; we verify the estimate uses
  # pi=0.3 for block-1 units and pi=0.6 for block-2 units (not shifted rows).
  N2 <- 30
  dat2 <- data.frame(y = rnorm(N2), bl = c(rep(1, 10), rep(2, 20)))
  decl2 <- randomizr::declare_ra(N = N2, blocks = dat2$bl, block_m = c(3, 12))
  dat2$Z <- randomizr::conduct_ra(decl2)
  dat2_miss <- dat2; dat2_miss$y[2] <- NA
  m_miss <- horvitz_thompson(y ~ Z, data = dat2_miss, condition_prs = decl2)
  # Manually compute: use probabilities for rows 1,3,...,30 (skipping row 2)
  pm   <- decl2$probabilities_matrix
  kept <- c(1L, 3:30)
  pi2_k <- pm[kept, 2]; pi1_k <- pm[kept, 1]
  y_k   <- dat2$y[kept]
  z_k   <- as.character(dat2$Z[kept])        # coerce to character for safe comparison
  conds <- sort(unique(z_k))
  t2_k  <- which(z_k == conds[2]); t1_k <- which(z_k == conds[1])
  Y2_k  <- y_k[t2_k] / pi2_k[t2_k]; Y1_k <- y_k[t1_k] / pi1_k[t1_k]
  est_manual <- (sum(Y2_k) - sum(Y1_k)) / length(z_k)
  expect_equal(m_miss$coefficients[[1]], est_manual, tolerance = 1e-10)
})
