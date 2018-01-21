context("Helper - HT condition_pr_matrix")

test_that("condition_pr_matrix behaves as expected", {

  # Complete randomization
  n <- 5
  prs <- rep(0.4, times = n)
  comp_ra <- randomizr::declare_ra(N = n, prob = prs[1])
  perms <- randomizr::obtain_permutation_matrix(comp_ra)

  expect_equal(
    declaration_to_condition_pr_mat(comp_ra),
    permutations_to_condition_pr_mat(perms)
  )

  # Simple randomization
  prs <- rep(0.4, times = n)
  simp_ra <- randomizr::declare_ra(N = n, prob = prs[1], simple = TRUE)

  # perms <- randomizr::obtain_permutation_matrix(simp_ra)
  # Won't work because some permutations are more likely than others
  # So instead we just resample and set the tolerance
  perms <- replicate(15000, simp_ra$ra_function())
  # Won't be equal because some permutations are more likely than others in
  # this case
  expect_equal(
    declaration_to_condition_pr_mat(simp_ra),
    permutations_to_condition_pr_mat(perms),
    tolerance = 0.02
  )

  # Blocked case
  dat <- data.frame(
    bl = c("A", "B", "A", "B", "B", "B"),
    pr = c(0.5, 0.25, 0.5, 0.25, 0.25, 0.25)
  )

  bl_ra <- randomizr::declare_ra(blocks = dat$bl, block_m = c(1, 1))
  bl_perms <- randomizr::obtain_permutation_matrix(bl_ra)

  expect_equal(
    declaration_to_condition_pr_mat(bl_ra),
    permutations_to_condition_pr_mat(bl_perms)
  )

  # Cluster complete case
  dat <- data.frame(
    cl = c("A", "B", "A", "C", "A", "B")
  )

  cl_ra <- randomizr::declare_ra(clusters = dat$cl, m = 1)
  cl_perms <- randomizr::obtain_permutation_matrix(cl_ra)

  expect_equal(
    declaration_to_condition_pr_mat(cl_ra),
    permutations_to_condition_pr_mat(cl_perms)
  )

  # Blocked and clustered
  dat <- data.frame(
    bl = c("A", "B", "B", "B", "A", "A", "B", "B"),
    cl = c(1, 2, 3, 3, 4, 4, 5, 5)
  )

  bl_cl_ra <- randomizr::declare_ra(clusters = dat$cl, blocks = dat$bl, block_m  = c(1, 2))
  bl_cl_perms <- randomizr::obtain_permutation_matrix(bl_cl_ra)

  expect_equal(
    declaration_to_condition_pr_mat(bl_cl_ra),
    permutations_to_condition_pr_mat(bl_cl_perms)
  )

  # Custom case
  cust_perms <- cbind(c(1, 0, 1, 0), c(1, 1, 0, 0))
  cust_ra <- randomizr::declare_ra(permutation_matrix = cust_perms)

  expect_equal(
    declaration_to_condition_pr_mat(cust_ra),
    permutations_to_condition_pr_mat(cust_perms)
  )

})
