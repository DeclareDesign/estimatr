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
    coef(horvitz_thompson(
      y ~ z,
      condition_prs = ps,
      data = dat
    )),
    coef(difference_in_means(
      y ~ z,
      data = dat
    ))
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
    ra_declaration = simp_decl,
    return_condition_pr_mat = TRUE
  )

  # Also with no SEs
  ht_simp_no <- horvitz_thompson(
    y ~ z,
    data = dat,
    ra_declaration = simp_decl,
    return_condition_pr_mat = TRUE,
    se_type = "none"
  )

  expect_equal(
    ht_simp$coefficients,
    ht_simp_no$coefficients
  )

  expect_equivalent(
    tidy(ht_simp_no)[c("std.error", "p.value", "ci.lower", "ci.upper")],
    rep(NA_real_, 4)
  )

  # Works with constant effects assumption
  ht_const <- horvitz_thompson(
    y ~ z,
    data = dat,
    ra_declaration = simp_decl,
    se_type = "constant"
  )

  # picks out right declaration
  ht_rev <- horvitz_thompson(
    y ~ z,
    data = dat,
    condition1 = 1,
    condition2 = 0,
    ra_declaration = simp_decl,
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
      ra_declaration = simp_decl
    )
  )

  expect_equal(
    tidy(ht_simp)[, c("estimate", "std.error")],
    tidy(ht_rev)[, c("estimate", "std.error")] * c(-1, 1)
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
    ht_comp <- horvitz_thompson(y ~ z_comp, data = dat, simple = FALSE),
    horvitz_thompson(y ~ z_comp, data = dat, ra_declaration = comp_decl)
  )
  expect_equal(
    ht_comp,
    horvitz_thompson(y ~ z_comp, data = dat, simple = FALSE, condition_prs = pr_comp)
  )

  # Also with no SEs
  ht_comp_no <- horvitz_thompson(y ~ z_comp, data = dat, simple = FALSE, se_type = "none")

  expect_equal(
    ht_comp$coefficients,
    ht_comp_no$coefficients
  )

  expect_equivalent(
    tidy(ht_comp_no)[c("std.error", "p.value", "ci.lower", "ci.upper")],
    rep(NA_real_, 4)
  )

  # error if you pass wrong prs
  dat$pr_wrong <- dat$pr_comp
  dat$pr_wrong[1] <- 0.5
  expect_error(
    horvitz_thompson(y ~ z_comp, data = dat, simple = FALSE, condition_prs = pr_wrong),
    "Treatment probabilities must be fixed for complete randomized designs"
  )

  # Works without data frame!
  ht_with <- with(
    dat,
    horvitz_thompson(y ~ z_comp, simple = FALSE, condition_prs = pr_comp)
  )

  pr_comp <- dat$pr_comp
  y <- dat$y
  z_comp <- dat$z_comp
  ht_glob <- horvitz_thompson(y ~ z_comp, simple = FALSE, condition_prs = pr_comp)

  expect_equal(
    ht_with,
    ht_glob
  )

  # with declaration
  ht_nod <- horvitz_thompson(y ~ z_comp, ra_declaration = comp_decl)
  ht_d <- horvitz_thompson(y ~ z_comp, data = dat, ra_declaration = comp_decl)
  expect_equal(
    tidy(ht_nod),
    tidy(ht_d)
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
  ht_crs_decl <- horvitz_thompson(y ~ z, data = dat, ra_declaration = clust_crs_decl)

  expect_true(
    !is.na(ht_crs_decl$coefficients)
  )

  expect_equivalent(
    ht_crs_decl$df,
    NA
  )

  # Also with no SEs
  ht_crs_decl_no <- horvitz_thompson(y ~ z, data = dat, ra_declaration = clust_crs_decl, se_type = "none")

  expect_equal(
    ht_crs_decl$coefficients,
    ht_crs_decl_no$coefficients
  )

  expect_equivalent(
    tidy(ht_crs_decl_no)[c("std.error", "p.value", "ci.lower", "ci.upper")],
    rep(NA_real_, 4)
  )

  # Can infer probabilities as well
  expect_equal(
    ht_crs_decl,
    horvitz_thompson(y ~ z, data = dat, clusters = cl, simple = FALSE)
  )

  # And constant effects error for non-simple designs
  expect_error(
    horvitz_thompson(y ~ z, data = dat, ra_declaration = clust_crs_decl, se_type = "constant"),
    "`se_type` = 'constant' only supported for simple random"
  )

  ## Simple random sample, clustered
  clust_srs_decl <- randomizr::declare_ra(
    N = nrow(dat),
    clusters = dat$cl,
    prob = 0.4,
    simple = TRUE
  )

  # With declaration
  # Regular SE using Young's inequality
  ht_srs_decl <- horvitz_thompson(y ~ z, data = dat, ra_declaration = clust_srs_decl)

  # Also with no SEs
  ht_srs_decl_no <- horvitz_thompson(y ~ z, data = dat, ra_declaration = clust_srs_decl, se_type = "none")

  expect_equal(
    ht_srs_decl$coefficients,
    ht_srs_decl_no$coefficients
  )

  expect_equivalent(
    tidy(ht_srs_decl_no)[c("std.error", "p.value", "ci.lower", "ci.upper")],
    rep(NA_real_, 4)
  )

  # Not the same because second doesn't know it's clustered!
  # Just passing mat
  clust_srs_mat <- declaration_to_condition_pr_mat(clust_srs_decl)
  expect_is(
    all.equal(
      ht_srs_decl,
      ht_srs_nodecl <- horvitz_thompson(y ~ z, data = dat, condition_pr_mat = clust_srs_mat)
    ),
    "character"
  )

  # Also with no SEs
  ht_srs_nodecl_no <-  horvitz_thompson(y ~ z, data = dat, condition_pr_mat = clust_srs_mat, se_type = "none")

  expect_equal(
    ht_srs_nodecl$coefficients,
    ht_srs_nodecl_no$coefficients
  )

  # works if I also pass cluster
  expect_identical(
    ht_srs_decl,
    ht_srs_cl <- horvitz_thompson(y ~ z, data = dat, clusters = cl, condition_pr_mat = clust_srs_mat)
  )

  # Also with no SEs
  ht_srs_cl_no <- horvitz_thompson(y ~ z, data = dat, clusters = cl, condition_pr_mat = clust_srs_mat, se_type = "none")

  expect_equal(
    ht_srs_cl$coefficients,
    ht_srs_cl_no$coefficients
  )

  # Can infer from number of treated clusters per block the treatment pr
  clbl_dat <- data.frame(
    cl_new = cl_new <- c(1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8),
    bl = rep(1:3, each = 4),
    y = rnorm(12)
  )
  # pr = 0.25 in first, 0.5 in second
  blcl_ra <- randomizr::declare_ra(blocks = clbl_dat$bl, clusters = clbl_dat$cl_new, block_m = c(1, 2, 1))
  clbl_dat$z_clbl <- randomizr::conduct_ra(blcl_ra)
  expect_equivalent(
    horvitz_thompson(y ~ z_clbl, data = clbl_dat, ra_declaration = blcl_ra),
    horvitz_thompson(y ~ z_clbl, data = clbl_dat, blocks = bl, clusters = cl_new)
  )

  # should work with just a column if SRS!
  dat$ps <- 0.4
  expect_identical(
    ht_srs_decl,
    ht_srs_prs <- horvitz_thompson(y ~ z, data = dat, clusters = cl, condition_prs = ps)
  )

  # Also with no SEs
  ht_srs_prs_no <- horvitz_thompson(y ~ z, data = dat, clusters = cl, condition_prs = ps, se_type = "none")

  expect_equal(
    ht_srs_prs$coefficients,
    ht_srs_prs_no$coefficients
  )

  # And constant effects
  # Only work for simple for now
  expect_error(
    horvitz_thompson(y ~ z, data = dat, ra_declaration = clust_srs_decl, se_type = "constant"),
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
  decl$probabilities_matrix
  nrow(missing_dat)
  expect_error(
    horvitz_thompson(y ~ z, data = missing_dat, ra_declaration = decl),
    NA
  )
  # Test that we didn't edit the declaration in the users env
  # Should work a second time
  expect_error(
    horvitz_thompson(y ~ z, data = missing_dat, ra_declaration = decl),
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
  ht_declare_bl <- horvitz_thompson(y ~ z, data = dat, ra_declaration = bl_ra)
  # This estimates the treatment effect at once using only condition_pr_mat
  ht_condmat_bl <- horvitz_thompson(y ~ z, data = dat, condition_pr_mat = bl_pr_mat)

  expect_equivalent(
    tidy(ht_declare_bl),
    tidy(ht_condmat_bl)
  )

  # Also with no SEs
  ht_declare_bl_no <- horvitz_thompson(y ~ z, data = dat, ra_declaration = bl_ra, se_type = "none")
  ht_condmat_bl_no <- horvitz_thompson(y ~ z, data = dat, condition_pr_mat = bl_pr_mat, se_type = "none")

  expect_equal(
    ht_declare_bl$coefficients,
    ht_declare_bl_no$coefficients
  )
  expect_equal(
    ht_condmat_bl$coefficients,
    ht_condmat_bl_no$coefficients
  )

  dat$mps <- rep(1:20, each = 2)
  mp_ra <- randomizr::declare_ra(blocks = dat$mps)
  dat$z <- randomizr::conduct_ra(mp_ra)
  mp_pr_mat <- declaration_to_condition_pr_mat(mp_ra)

  ht_declare_mp <- horvitz_thompson(y ~ z, data = dat, ra_declaration = mp_ra)
  # This estimates the treatment effect at once using only condition_pr_mat
  ht_condmat_mp <- horvitz_thompson(y ~ z, data = dat, condition_pr_mat = mp_pr_mat)

  expect_equivalent(
    tidy(ht_declare_mp),
    tidy(ht_condmat_mp)
  )

  # block messages when passing with simple = TRUE flag, not otherwise
  dat$p <- tapply(dat$z, dat$bl, mean)[dat$bl]
  expect_message(
    ht_declare_mp <- horvitz_thompson(y ~ z, data = dat, blocks = bl, condition_prs = p, simple = TRUE),
    "Assuming complete random assignment of clusters within blocks."
  )

  expect_message(
    ht_declare_mp <- horvitz_thompson(y ~ z, data = dat, blocks = bl, condition_prs = p, simple = FALSE),
    NA
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
    horvitz_thompson(y ~ z, data = dat, condition_prs = ps, ra_declaration = decl),
    "Cannot use `ra_declaration` with any of"
  )

  expect_error(
    horvitz_thompson(y ~ z, data = dat, condition_pr_mat = declaration_to_condition_pr_mat(decl), ra_declaration = decl),
    "Cannot use `ra_declaration` with any of"
  )

  expect_error(
    horvitz_thompson(y ~ z + x, data = dat, ra_declaration = decl),
    "must have only one variable on the right-hand side"
  )

  expect_error(
    horvitz_thompson(y ~ z, data = dat, ra_declaration = randomizr::declare_ra(N = n + 1, prob = 0.4)),
    "variable lengths differ"
  )

  ht_o <- horvitz_thompson(y ~ z, data = dat, ci = FALSE)
  expect_equivalent(
    as.matrix(tidy(horvitz_thompson(y ~ z, data = dat, ci = FALSE))[, c("p.value", "ci.lower", "ci.upper")]),
    matrix(NA, nrow = 1, ncol = 3)
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


  expect_equivalent(coef(ht_const_1), mean(dat$y))
  expect_equivalent(ht_const_1$std.error, 1 / (nrow(dat)) * sqrt(sum(dat$y ^ 2)))


  expect_equivalent(
    ht_const_1$df,
    NA
  )

  ht_const <- horvitz_thompson(
    y ~ z_const,
    data = dat,
    condition_prs = ps
  )

  expect_equivalent(coef(ht_const), mean(dat$y / dat$ps))
  expect_equivalent(ht_const$std.error, 1 / (nrow(dat)) * sqrt(sum((dat$y / dat$ps) ^ 2)))

  ## Blocks and all are treated
  ht_block <- horvitz_thompson(
    y ~ z_const,
    data = dat,
    blocks = bl,
    condition_prs = ps,
    return_condition_pr_mat = TRUE
  )

  # with blocks SE is different because not simple any more
  expect_equivalent(coef(ht_block), mean(dat$y / dat$ps))
  # expect_equal(ht_block$std.error, 1/(nrow(dat)) * sqrt(sum((dat$y / dat$ps)^2)))

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

  expect_identical(ht_zero$term, "z0")

  # Drop name if they specify the only treatment as condition1
  ht_rev <- horvitz_thompson(
    y ~ z,
    data = dat,
    blocks = bl,
    condition1 = 0,
    condition_prs = rep(0.5, nrow(dat))
  )

  expect_identical(ht_rev$term, "z")

  # This is only true because condition prs are 0.5
  expect_identical(
    tidy(ht_zero)[c("estimate", "std.error")],
    tidy(ht_rev)[c("estimate", "std.error")] * c(-1, 1)
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

test_that("multi-valued treatments not allowed in ra_declaration", {
  dat <- data.frame(
    y = rnorm(20),
    ps = 0.4
  )

  decl_multi <- randomizr::declare_ra(N = 20, prob_each = c(0.4, 0.4, 0.2))
  dat$z <- randomizr::conduct_ra(decl_multi)

  expect_error(
    horvitz_thompson(y ~ z, data = dat, ra_declaration = decl_multi),
    "Cannot use horvitz_thompson\\(\\) with a `ra_declaration` with"
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
