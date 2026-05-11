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

test_that("HT missing Y drops the row", {
  decl <- randomizr::declare_ra(N = N, m = 16)
  dat$Z <- randomizr::conduct_ra(decl)
  dat_miss <- dat; dat_miss$y[5] <- NA
  m_miss <- horvitz_thompson(y ~ Z, data = dat_miss, condition_prs = decl)
  m_drop <- horvitz_thompson(y ~ Z, data = dat[-5, ], condition_prs = decl)
  expect_equal(m_miss$coefficients[[1]], m_drop$coefficients[[1]], tolerance = 1e-10)
})
