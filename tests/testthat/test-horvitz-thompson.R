context("Horvitz Thompson")


test_that("Horvitz-Thompson matches d-i-m under certain conditions", {

  n <- 4
  dat <- data.frame(y0 = rnorm(n),
                    z = rep(0:1, each = n/2),
                    ps = rep(0.5, n))

  dat$y1 <- dat$y0 + 0.43
  dat$y <- ifelse(dat$z, dat$y1, dat$y0)

  expect_equal(
    horvitz_thompson(y ~ z,
                     condition_prs = ps,
                     data = dat)$est,
    difference_in_means(y ~ z,
                        data = dat)$est
  )

})

test_that("Horvitz-Thompson works with clustered data", {

  n <- 8
  dat <- data.frame(y = rnorm(n),
                    cl = rep(1:4, each = 2))

  ## Complete random sample, clustered
  clust_crs_decl <- randomizr::declare_ra(N = nrow(dat), clust_var = dat$cl, prob = 0.4)
  dat$z <- randomizr::conduct_ra(clust_crs_decl)

  # Regular SE using Young's inequality
  ht_crs_decl <- horvitz_thompson(y ~ z, data = dat, declaration = clust_crs_decl)
  # Also can just pass probability matrix
  clust_crs_mat <- declaration_to_condition_pr_mat(clust_crs_decl)
  ht_crs_mat <- horvitz_thompson(y ~ z, data = dat, condition_pr_mat = clust_crs_mat)

  expect_identical(
    ht_crs_decl,
    ht_crs_mat
  )

  # Also should be same with collapsed totals
  expect_identical(
    horvitz_thompson(y ~ z, data = dat, declaration = clust_crs_decl, collapsed = T),
    horvitz_thompson(y ~ z, data = dat, clusters = cl, condition_pr_mat = clust_crs_mat, collapsed = T)
  )

  # Have to specify clusters for collapsed estimator if you don't pass declaration
  expect_error(
    horvitz_thompson(y ~ z, data = dat, condition_pr_mat = clust_crs_mat, collapsed = T),
    "collapsed estimator only works"
  )

  # And constant effects
  expect_identical(
    horvitz_thompson(y ~ z, data = dat, declaration = clust_crs_decl, se_type = 'constant'),
    horvitz_thompson(y ~ z, data = dat, condition_pr_mat = clust_crs_mat, se_type = 'constant')
  )

  ## Simple random sample, clustered
  clust_srs_decl <- randomizr::declare_ra(N = nrow(dat),
                                          clust_var = dat$cl,
                                          prob = 0.4,
                                          simple = T)

  # With declaration
  # TODO Update once new randomizr is on CRAN
  # For now, the old version of declare_ra stores no information about whether
  # a cluster randomized design is simple or complete. Have to
  # comment this out for now.
  if (FALSE) {
    # Regular SE using Young's inequality
    ht_srs_decl <- horvitz_thompson(y ~ z, data = dat, declaration = clust_srs_decl)
    # Also can just pass probability matrix
    clust_srs_mat <- declaration_to_condition_pr_mat(clust_srs_decl)

    expect_identical(
      ht_srs_decl,
      horvitz_thompson(y ~ z, data = dat, condition_pr_mat = clust_srs_mat)
    )

    # should work with just a column if SRS!
    dat$ps <- 0.4
    expect_identical(
      ht_srs_decl,
      horvitz_thompson(y ~ z, data = dat, clusters = cl, condition_prs = ps)
    )

    # Also should be same with collapsed totals
    # matrix approach
    ht_cl_srs_collapsed <- horvitz_thompson(y ~ z, data = dat, declaration = clust_srs_decl, collapsed = T)
    expect_identical(
      ht_cl_srs_collapsed,
      horvitz_thompson(y ~ z, data = dat, clusters = cl, condition_pr_mat = clust_srs_mat, collapsed = T)
    )

    # condition var name approach
    expect_identical(
      ht_cl_srs_collapsed,
      horvitz_thompson(y ~ z, data = dat, clusters = cl, condition_prs = ps, collapsed = T)
    )

    # And constant effects
    # matrix approach
    ht_cl_srs_const <- horvitz_thompson(y ~ z, data = dat, declaration = clust_srs_decl, se_type = 'constant')
    expect_identical(
      ht_cl_srs_const,
      horvitz_thompson(y ~ z, data = dat, condition_pr_mat = clust_srs_mat, se_type = 'constant')
    )

    # condition var name approach
    expect_identical(
      ht_cl_srs_const,
      horvitz_thompson(y ~ z, data = dat, clusters = cl, condition_prs = ps, se_type = 'constant')
    )
  }


})

# test missingness works as expected
# test blocks in the data

# test auxiliary funtions (get condition pr)
test_that("gen_pr_matrix_complete works as expected", {
  # TODO test other methods
  n <- 5
  prs <- rep(0.4, times = n)
  pr_mat <- gen_pr_matrix_complete(prs)

  # FALSE until randomizr on CRAN
  if (FALSE) {
    perms <- randomizr::obtain_permutation_matrix(randomizr::declare_ra(N = n, prob = prs[1]))

    # From Chris Kennedy (https://github.com/ck37)
    # for the htestimate package (https://github.com/ck37/htestimate)
    stacked_inds <- matrix(nrow = 2 * n, ncol = ncol(perms))

    for (assign in 1:2) {
      indicator_matrix <- as.numeric(perms == (assign - 1))
      stacked_inds[(n*(assign-1)+1):(n*assign), ] <- indicator_matrix
    }

    # Use the stacked indicator matrices to calculate the probability matrix.
    result <- stacked_inds %*% t(stacked_inds) / ncol(perms)
    expect_equivalent(pr_mat, result)

  }
})

# errors when arguments are passed that shouldn't be together
test_that("Horvitz-Thompson properly checks arguments", {

  n <- 8
  dat <- data.frame(y = rnorm(n),
                    ps = 0.4,
                    z = sample(rep(0:1, each = n/2)),
                    x = runif(n))
  decl <- randomizr::declare_ra(N = n, prob = 0.4, simple = F)

  expect_error(
    horvitz_thompson(y ~ z, data = dat, condition_prs = ps, declaration = decl),
    "Cannot use declaration with any of"
  )

  expect_error(
    horvitz_thompson(y ~ z, data = dat, condition_pr_mat = declaration_to_condition_pr_mat(decl), declaration = decl),
    "Cannot use declaration with any of"
  )

  expect_error(
    horvitz_thompson(y ~ z + x, data = dat, declaration = decl),
    "formula"
  )

  expect_error(
    horvitz_thompson(y ~ z, data = dat, declaration = randomizr::declare_ra(N = n+1, prob = 0.4)),
    "N|declaration"
  )

})

