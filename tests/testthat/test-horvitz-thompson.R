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
                     condition_pr_variable_name = ps,
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
  ht_crs_mat <- horvitz_thompson(y ~ z, data = dat, condition_pr_matrix = clust_crs_mat)

  # Also should be same with collapsed totals
  expect_identical(
    horvitz_thompson(y ~ z, data = dat, declaration = clust_crs_decl, collapsed = T),
    horvitz_thompson(y ~ z, data = dat, condition_pr_matrix = clust_crs_mat, collapsed = T)
  )

  # And constant effects
  expect_identical(
    horvitz_thompson(y ~ z, data = dat, declaration = clust_crs_decl, se_type = 'constant'),
    horvitz_thompson(y ~ z, data = dat, condition_pr_matrix = clust_crs_mat, se_type = 'constant')
  )

  ## Simple random sample, clustered
  clust_srs_decl <- randomizr::declare_ra(N = nrow(dat),
                                          clust_var = dat$cl,
                                          prob = 0.4,
                                          simple = F)

  # With declaration
  # Regular SE using Young's inequality
  ht_srs_decl <- horvitz_thompson(y ~ z, data = dat, declaration = clust_srs_decl)
  # Also can just pass probability matrix
  clust_srs_mat <- declaration_to_condition_pr_mat(clust_srs_decl)

  expect_identical(
    ht_srs_decl,
    horvitz_thompson(y ~ z, data = dat, condition_pr_matrix = clust_srs_mat)
  )

  # should work with just a column if SRS!
  dat$ps <- 0.4
  clust_srs_mat
  expect_identical(
    ht_srs_decl,
    horvitz_thompson(y ~ z, data = dat, cluster_variable_name = cl, condition_pr_variable_name = ps)
  )

  # Also should be same with collapsed totals
  # matrix approach
  ht_cl_srs_collapsed <- horvitz_thompson(y ~ z, data = dat, declaration = clust_srs_decl, collapsed = T)
  expect_identical(
    ht_cl_srs_collapsed,
    horvitz_thompson(y ~ z, data = dat, condition_pr_matrix = clust_srs_mat, collapsed = T)
  )

  # condition var name approach
  expect_identical(
    ht_cl_srs_collapsed,
    horvitz_thompson(y ~ z, data = dat, cluster_variable_name = cl, condition_pr_variable_name = ps, collapsed = T)
  )

  # And constant effects
  # matrix approach
  ht_cl_srs_const <- horvitz_thompson(y ~ z, data = dat, declaration = clust_srs_decl, se_type = 'constant')
  expect_identical(
    ht_cl_srs_const,
    horvitz_thompson(y ~ z, data = dat, condition_pr_matrix = clust_srs_mat, se_type = 'constant')
  )

  # condition var name approach
  expect_identical(
    ht_cl_srs_const,
    horvitz_thompson(y ~ z, data = dat, cluster_variable_name = cl, condition_pr_variable_name = ps, se_type = 'constant')
  )

})

# test missingness works as expected
# errors when arguments are passed that shouldn't be together
