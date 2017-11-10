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
  clust_srs_mat
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

})

# test missingness works as expected
# test blocks in the data
# test auxiliary funtions (get condition pr)

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

