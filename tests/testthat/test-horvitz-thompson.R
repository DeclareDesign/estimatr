context("Estimator - horvitz_thompson")


test_that("Horvitz-Thompson matches d-i-m under certain conditions", {
  n <- 4
  dat <- data.frame(
    y0 = rnorm(n),
    z = rep(0:1, each = n / 2),
    ps = rep(0.5, n)
  )

  dat$y1 <- dat$y0 + 0.43
  dat$y <- ifelse(dat$z, dat$y1, dat$y0)

  expect_equal(
    horvitz_thompson(
      y ~ z,
      condition_prs = ps,
      data = dat
    )$coefficients,
    difference_in_means(
      y ~ z,
      data = dat
    )$coefficients
  )
})

test_that("Horvitz-Thompson works in simple case", {
  n <- 40
  dat <- data.frame(
    y = rnorm(n)
  )
  simp_decl <- randomizr::declare_ra(N = n, prob = 0.4, simple = T)
  dat$z <- randomizr::conduct_ra(simp_decl)

  ht_simp <- horvitz_thompson(
    y ~ z,
    data = dat,
    declaration = simp_decl,
    return_condition_pr_mat = TRUE
  )

  # Works with constant effects assumption
  ht_const <- horvitz_thompson(
    y ~ z,
    data = dat,
    declaration = simp_decl,
    se_type = "constant"
  )

  # picks out right declaration
  ht_rev <- horvitz_thompson(
    y ~ z,
    data = dat,
    condition1 = 1,
    condition2 = 0,
    declaration = simp_decl,
    return_condition_pr_mat = TRUE
  )

  # Fails properly if condition in treatment but not in declaration
  dat$z[1] <- 2
  expect_error(
    horvitz_thompson(
      y ~ z,
      data = dat,
      condition1 = 0,
      condition2 = 2,
      declaration = simp_decl
    )
  )

  expect_equal(
    tidy(ht_simp)[, c("coefficients", "se")],
    tidy(ht_rev)[, c("coefficients", "se")] * c(-1, 1)
  )

  # Simple designs needn't use condition matrix as joint prs are product of marginals
  expect_equal(
    ht_simp$condition_pr_mat,
    NULL
  )

  # complete randomization works as well
  comp_decl <- randomizr::declare_ra(N = n, prob = 0.4, simple = FALSE)
  dat$z_comp <- randomizr::conduct_ra(comp_decl)
  dat$pr_comp <- 0.4
  # We can learn it! or you can tell us
  expect_equal(
    horvitz_thompson(y ~ z_comp, data = dat, simple = FALSE),
    horvitz_thompson(y ~ z_comp, data = dat, declaration = comp_decl)
  )
  expect_equal(
    horvitz_thompson(y ~ z_comp, data = dat, simple = FALSE),
    horvitz_thompson(y ~ z_comp, data = dat, simple = FALSE, condition_prs = pr_comp)
  )

  # error if you pass wrong prs
  dat$pr_wrong <- dat$pr_comp
  dat$pr_wrong[1] <- 0.5
  expect_error(
    horvitz_thompson(y ~ z_comp, data = dat, simple = FALSE, condition_prs = pr_wrong),
    "Treatment probabilities must be fixed for complete randomized designs"
  )
})

test_that("Horvitz-Thompson works with clustered data", {
  n <- 8
  dat <- data.frame(
    y = rnorm(n),
    cl = rep(1:4, each = 2)
  )

  ## Complete random sample, clustered
  clust_crs_decl <- randomizr::declare_ra(N = nrow(dat), clusters = dat$cl, prob = 0.5)
  dat$z <- randomizr::conduct_ra(clust_crs_decl)

  # Regular SE using Young's inequality
  ht_crs_decl <- horvitz_thompson(y ~ z, data = dat, declaration = clust_crs_decl)

  expect_equal(
    ht_crs_decl$df,
    NA
  )

  # Can infer probabilities as well
  expect_equal(
    ht_crs_decl,
    horvitz_thompson(y ~ z, data = dat, clusters = cl, simple = FALSE)
  )

  # And constant effects error for non-simple designs
  expect_error(
    horvitz_thompson(y ~ z, data = dat, declaration = clust_crs_decl, se_type = "constant"),
    "`se_type` = 'constant' only supported for simple random"
  )

  ## Simple random sample, clustered
  clust_srs_decl <- randomizr::declare_ra(
    N = nrow(dat),
    clusters = dat$cl,
    prob = 0.4,
    simple = T
  )

  # With declaration
  # Regular SE using Young's inequality
  ht_srs_decl <- horvitz_thompson(y ~ z, data = dat, declaration = clust_srs_decl)

  # Not the same because second doesn't know it's clustered!
  # Just passing mat
  clust_srs_mat <- declaration_to_condition_pr_mat(clust_srs_decl)
  expect_is(
    all.equal(
      ht_srs_decl,
      horvitz_thompson(y ~ z, data = dat, condition_pr_mat = clust_srs_mat)
    ),
    "character"
  )

  # works if I also pass cluster
  expect_identical(
    ht_srs_decl,
    horvitz_thompson(y ~ z, data = dat, clusters = cl, condition_pr_mat = clust_srs_mat)
  )

  # Can infer from number of treated clusters per block the treatment pr
  clbl_dat <- data.frame(
    cl_new = cl_new <- c(1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8),
    bl = rep(1:3, each = 4),
    y = rnorm(12)
  )
  # pr = 0.25 in first, 0.5 in second
  blcl_ra <- randomizr::declare_ra(blocks = clbl_dat$bl, clusters = clbl_dat$cl_new, m = c(1, 2, 1))
  clbl_dat$z_clbl <- blcl_ra$ra_function()
  expect_equivalent(
    horvitz_thompson(y ~ z_clbl, data = clbl_dat, declaration = blcl_ra),
    horvitz_thompson(y ~ z_clbl, data = clbl_dat, blocks = bl, clusters = cl_new)
  )

  # should work with just a column if SRS!
  dat$ps <- 0.4
  expect_identical(
    ht_srs_decl,
    horvitz_thompson(y ~ z, data = dat, clusters = cl, condition_prs = ps)
  )


  # And constant effects
  # Only work for simple for now
  expect_error(
    horvitz_thompson(y ~ z, data = dat, declaration = clust_srs_decl, se_type = "constant"),
    "`se_type` = 'constant' only supported for simple random designs at the moment"
  )

  # Fails with condition_pr varying within cluster
  dat$p_wrong <- dat$ps
  dat$p_wrong[1] <- 0.545

  expect_error(
    horvitz_thompson(y ~ z, data = dat, clusters = cl, condition_prs = p_wrong),
    "`condition_prs` must be constant within `cluster`"
  )

  # Or pr outside of 0 1
  dat$p_wrong[1] <- 1.5
  expect_error(
    horvitz_thompson(y ~ z, data = dat, clusters = cl, condition_prs = p_wrong),
    "`condition_prs` must be a vector of positive values no greater than 1"
  )

  # or treatment varying within a cluster
  dat$z_wrong <- dat$z
  dat$z_wrong[1:2] <- c(0, 1)
  table(dat$z_wrong, dat$cl)
  expect_error(
    horvitz_thompson(y ~ z_wrong, data = dat, clusters = cl, condition_prs = ps),
    "Treatment condition must be constant within `clusters`"
  )
})

# TODO test missingness works as expected
test_that("Horvitz-Thompson works with missingness", {
  n <- 40
  dat <- data.frame(
    y = rnorm(n),
    bl = rep(1:10, each = 4),
    ps = 0.35
  )

  decl <- randomizr::declare_ra(n, prob = 0.35)
  dat$z <- randomizr::conduct_ra(decl)
  missing_dat <- dat
  missing_dat$y[1] <- NA

  expect_error(
    horvitz_thompson(y ~ z, data = missing_dat, declaration = decl),
    NA
  )

  missing_dat$ps[2] <- NA
  dat$drop_these <- c(1, 1, rep(0, times = n - 2))

  expect_warning(
    ht_miss <- horvitz_thompson(y ~ z, data = missing_dat, condition_prs = ps),
    "missingness in the condition_pr"
  )

  expect_equal(
    horvitz_thompson(y ~ z, data = dat, condition_prs = ps, subset = drop_these == 0),
    ht_miss
  )
})

# test blocks in the data
test_that("Estimating Horvitz-Thompson can be done two ways with blocks", {
  n <- 40
  dat <- data.frame(
    y = rnorm(n),
    bl = rep(1:10, each = 4)
  )

  bl_ra <- randomizr::declare_ra(blocks = dat$bl)
  dat$z <- randomizr::conduct_ra(bl_ra)
  bl_pr_mat <- declaration_to_condition_pr_mat(bl_ra)

  # This creates estimates within blocks and then joins them together using the common
  # formula
  ht_declare_bl <- horvitz_thompson(y ~ z, data = dat, declaration = bl_ra)
  # This estimates the treatment effect at once using only condition_pr_mat
  ht_condmat_bl <- horvitz_thompson(y ~ z, data = dat, condition_pr_mat = bl_pr_mat)

  expect_equivalent(
    tidy(ht_declare_bl),
    tidy(ht_condmat_bl)
  )

  dat$mps <- rep(1:20, each = 2)
  mp_ra <- randomizr::declare_ra(blocks = dat$mps)
  dat$z <- randomizr::conduct_ra(mp_ra)
  mp_pr_mat <- declaration_to_condition_pr_mat(mp_ra)

  ht_declare_mp <- horvitz_thompson(y ~ z, data = dat, declaration = mp_ra)
  # This estimates the treatment effect at once using only condition_pr_mat
  ht_condmat_mp <- horvitz_thompson(y ~ z, data = dat, condition_pr_mat = mp_pr_mat)

  expect_equivalent(
    tidy(ht_declare_mp),
    tidy(ht_condmat_mp)
  )
})

# errors when arguments are passed that shouldn't be together
test_that("Horvitz-Thompson properly checks arguments and data", {
  n <- 8
  dat <- data.frame(
    y = rnorm(n),
    ps = 0.4,
    z = sample(rep(0:1, each = n / 2)),
    x = runif(n),
    cl = rep(1:4, each = 2),
    bl = rep(1:2, each = 4)
  )
  decl <- randomizr::declare_ra(N = n, prob = 0.4, simple = FALSE)

  # default is mean(ps)
  expect_identical(
    horvitz_thompson(y ~ z, data = dat),
    horvitz_thompson(y ~ z, data = dat, condition_prs = rep(0.5, times = nrow(dat)))
  )

  expect_error(
    horvitz_thompson(y ~ z, data = dat, condition_prs = ps, declaration = decl),
    "Cannot use `declaration` with any of"
  )

  expect_error(
    horvitz_thompson(y ~ z, data = dat, condition_pr_mat = declaration_to_condition_pr_mat(decl), declaration = decl),
    "Cannot use `declaration` with any of"
  )

  expect_error(
    horvitz_thompson(y ~ z + x, data = dat, declaration = decl),
    "must have only one variable on the right-hand side"
  )

  expect_error(
    horvitz_thompson(y ~ z, data = dat, declaration = randomizr::declare_ra(N = n + 1, prob = 0.4)),
    "N|declaration"
  )

  ht_o <- horvitz_thompson(y ~ z, data = dat, ci = FALSE)
  expect_equivalent(
    as.matrix(tidy(horvitz_thompson(y ~ z, data = dat, ci = FALSE))[, c("p", "ci_lower", "ci_upper")]),
    matrix(NA, nrow = 1, ncol = 3)
  )

  # Reserved variable names
  dat[[".clusters_ddinternal"]] <- 1
  expect_error(
    horvitz_thompson(
      y ~ z,
      data = dat,
      declaration = randomizr::declare_ra(clusters = dat$cl)
    ),
    ".clusters_ddinternal"
  )

  dat[[".blocks_ddinternal"]] <- 1
  expect_error(
    horvitz_thompson(
      y ~ z,
      data = dat,
      declaration = randomizr::declare_ra(blocks = dat$bl)
    ),
    ".blocks_ddinternal"
  )

  dat[[".treatment_prob_ddinternal"]] <- 1
  expect_error(
    horvitz_thompson(
      y ~ z,
      data = dat,
      declaration = randomizr::declare_ra(N = n)
    ),
    ".treatment_prob_ddinternal"
  )


  # condition pr mat is the wrong size
  expect_error(
    horvitz_thompson(
      y ~ z,
      data = dat,
      condition_pr_mat = matrix(rnorm(4), 2, 2)
    ),
    "cleaning the data"
  )

  # subset and condition_pr_mat checked
})

test_that("Works without variation in treatment", {
  set.seed(1)
  dat <- data.frame(
    y = rnorm(20),
    bl = 1:5,
    ps = 0.4
  )

  # Simple case
  dat$z_const <- 1

  ht_const_1 <- horvitz_thompson(
    y ~ z_const,
    data = dat
  )

  ht_const_cond1 <- horvitz_thompson(
    y ~ z_const,
    data = dat,
    condition2 = 1
  )

  expect_equivalent(
    ht_const_1,
    ht_const_cond1
  )


  expect_equivalent(ht_const_1$coefficients, mean(dat$y))
  expect_equal(ht_const_1$se, 1 / (nrow(dat)) * sqrt(sum(dat$y ^ 2)))


  expect_equal(
    ht_const_1$df,
    NA
  )

  ht_const <- horvitz_thompson(
    y ~ z_const,
    data = dat,
    condition_prs = ps
  )

  expect_equivalent(ht_const$coefficients, mean(dat$y / dat$ps))
  expect_equal(ht_const$se, 1 / (nrow(dat)) * sqrt(sum((dat$y / dat$ps) ^ 2)))

  ## Blocks and all are treated
  ht_block <- horvitz_thompson(
    y ~ z_const,
    data = dat,
    blocks = bl,
    condition_prs = ps,
    return_condition_pr_mat = TRUE
  )

  # with blocks SE is different because not simple any more
  expect_equivalent(ht_block$coefficients, mean(dat$y / dat$ps))
  # expect_equal(ht_block$se, 1/(nrow(dat)) * sqrt(sum((dat$y / dat$ps)^2)))

  ## Blocks and some are treated!
  dat$z_diff <- as.numeric(dat$bl <= 2)
  ht_block <- horvitz_thompson(
    y ~ z_diff,
    data = dat,
    blocks = bl,
    condition_prs = rep(0.4, nrow(dat))
  )
  ht_block

  # With only one treatment, but value is 0, still put it as treatment!!
  # But note we leave a hint in the coefficient name
  dat$z <- 0
  ht_zero <- horvitz_thompson(
    y ~ z,
    data = dat,
    blocks = bl,
    condition_prs = rep(0.5, nrow(dat))
  )

  expect_identical(ht_zero$coefficient_name, "z0")

  # Drop name if they specify the only treatment as condition1
  ht_rev <- horvitz_thompson(
    y ~ z,
    data = dat,
    blocks = bl,
    condition1 = 0,
    condition_prs = rep(0.5, nrow(dat))
  )

  expect_identical(ht_rev$coefficient_name, "z")

  # This is only true because condition prs are 0.5
  expect_identical(
    tidy(ht_zero)[c("coefficients", "se")],
    tidy(ht_rev)[c("coefficients", "se")] * c(-1, 1)
  )

  # Some weird specifications that hit unusual parts of the variance
  cpm <- diag(0.5, nrow = 4, ncol = 4)
  y <- rnorm(2)
  t <- c(0, 1)
  expect_error(
    horvitz_thompson(y ~ t, condition_pr_mat = cpm),
    NA
  )
  t <- c(1, 1)
  expect_error(
    horvitz_thompson(y ~ t, condition_pr_mat = cpm),
    NA
  )

})

test_that("multi-valued treatments not allowed in declaration", {
  dat <- data.frame(
    y = rnorm(20),
    ps = 0.4
  )

  decl_multi <- randomizr::declare_ra(N = 20, prob_each = c(0.4, 0.4, 0.2))
  dat$z <- randomizr::conduct_ra(decl_multi)

  expect_error(
    horvitz_thompson(y ~ z, data = dat, declaration = decl_multi),
    "Cannot use horvitz_thompson\\(\\) with a `declaration` with"
  )

  # will work but you have to get the PRs right!
  ht_condition <- horvitz_thompson(
    y ~ z,
    data = dat,
    condition_prs = ps,
    condition1 = "T1",
    condition2 = "T2"
  )

  subdat <- dat[dat$z != "T3", ]
  ht_subdat <- horvitz_thompson(
    y ~ z,
    data = subdat,
    condition_prs = ps
  )

  ht_subset <- horvitz_thompson(
    y ~ z,
    data = dat,
    subset = z != "T3",
    condition_prs = ps
  )

  expect_equal(
    ht_condition,
    ht_subdat
  )
  expect_equal(
    ht_condition,
    ht_subset
  )
})
