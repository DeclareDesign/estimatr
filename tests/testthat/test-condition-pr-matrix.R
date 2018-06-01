context("Helper - HT condition_pr_matrix")
n <- 5


test_that("Checks class", {
    # Errors appropriately
    expect_error(
      declaration_to_condition_pr_mat(rbinom(5, 1, 0.5)),
      "`ra_declaration` must be an object of class 'ra_declaration'"
    )
})

test_that("Complete randomization", {

  prs <- rep(0.4, times = n)
  comp_ra <- randomizr::declare_ra(N = n, prob = prs[1])
  perms <- randomizr::obtain_permutation_matrix(comp_ra)
  expect_equal(
    declaration_to_condition_pr_mat(comp_ra),
    permutations_to_condition_pr_mat(perms)
  )

})

test_that("declaration to condition_pr_mat errors", {

  expect_error(
    declaration_to_condition_pr_mat(randomizr::declare_ra(N = n), 1, NULL),
    "Cannot have `condition2 == NULL`"
  )
  expect_error(
    declaration_to_condition_pr_mat(randomizr::declare_ra(N = n), NULL, 1),
    "Cannot have `condition1 == NULL`"
  )
  expect_error(
    declaration_to_condition_pr_mat(rbinom(5, 1, 0.5)),
    "`ra_declaration` must be an object of class 'ra_declaration'"
  )
})

test_that("condition args work properly", {

  # Condition args work properly
  mat01 <- declaration_to_condition_pr_mat(
    randomizr::declare_ra(N = n, prob = 0.4),
    0,
    1
  )
  mat10 <- declaration_to_condition_pr_mat(
    randomizr::declare_ra(N = n, prob = 0.4),
    1,
    0
  )

  # Diagonals are just flipped, check the names!
  # colnames(mat01)
  # colnames(mat10)
  expect_equal(mat01, mat10[rownames(mat01), colnames(mat01)])
})

test_that("Complete randomization with number of treated units not fixed", {

  #
  comp_odd_ra <- randomizr::declare_ra(N = 3, prob = 0.5)
  perms <- randomizr::obtain_permutation_matrix(comp_odd_ra)

  decl_cond_pr_mat <- declaration_to_condition_pr_mat(comp_odd_ra)

  # following passes so just use perms instead of get_perms
  # get_perms <- replicate(40000, conduct_ra(comp_odd_ra))
  # expect_true(
  #    max(permutations_to_condition_pr_mat(perms) -
  #         round(permutations_to_condition_pr_mat(get_perms), 3)) < 0.01
  # )

  expect_equal(
    decl_cond_pr_mat,
    permutations_to_condition_pr_mat(perms)
  )

})

test_that("Complete randomization with non 0.5 as remainder", {
  comp_odd_ra <- randomizr::declare_ra(N = 3, prob = 0.4)
  decl_cond_pr_mat <- declaration_to_condition_pr_mat(comp_odd_ra)

  set.seed(40)
  get_perms <- replicate(5000, randomizr::conduct_ra(comp_odd_ra))
  expect_equal(
    decl_cond_pr_mat,
    permutations_to_condition_pr_mat(get_perms),
    tolerance = 0.01
  )
})
test_that("Simple ra", {

  # Simple randomization
  prs <- rep(0.4, times = n)
  simp_ra <- randomizr::declare_ra(N = n, prob = prs[1], simple = TRUE)

  # perms <- randomizr::obtain_permutation_matrix(simp_ra)
  # Won't work because some permutations are more likely than others
  # So instead we just resample and set the tolerance
  perms <- replicate(10000, randomizr::conduct_ra(simp_ra))
  # Won't be equal because some permutations are more likely than others in
  # this case
  expect_equal(
    declaration_to_condition_pr_mat(simp_ra),
    permutations_to_condition_pr_mat(perms),
    tolerance = 0.02
  )
})
test_that("Blocked complete ra", {

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
})
test_that("Blocked complete ra with remainder", {
  dat <- data.frame(
    bl = c("A", "B", "A", "B", "B", "B"),
    pr = c(0.5, 0.25, 0.5, 0.25, 0.25, 0.25)
  )

  # with remainder
  bl <- c("A", "B", "A", "A", "B", "B") # Is this used anywhere?

  bl_ra <- randomizr::declare_ra(blocks = dat$bl, prob = 0.4)
  bl_perms <- replicate(5000, randomizr::conduct_ra(bl_ra))

  expect_equal(
    declaration_to_condition_pr_mat(bl_ra),
    permutations_to_condition_pr_mat(bl_perms),
    tolerance = 0.02
  )
})
test_that("Clustered complete ra", {

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

  # with remainder
  cl_ra <- randomizr::declare_ra(clusters = dat$cl, prob = 0.5)
  cl_perms <- randomizr::obtain_permutation_matrix(cl_ra)

  # lapply(1:ncol(cl_perms), function(x) table(dat$cl, cl_perms[, x]))
  expect_equal(
    declaration_to_condition_pr_mat(cl_ra),
    permutations_to_condition_pr_mat(cl_perms)
  )
})

test_that("Clustered ra", {

  # Cluster simple ? Should this be simple or no? --NJF
  dat <- data.frame(
    cl = c("A", "B", "A", "C", "A", "B")
  )

  dat$prs <- 0.3
  cl_simp_ra <- randomizr::declare_ra(clusters = dat$cl, prob = dat$prs[1])
  cl_simp_perms <- randomizr::obtain_permutation_matrix(cl_simp_ra)

  cl_simp_cpm <- declaration_to_condition_pr_mat(cl_simp_ra)

  expect_is(
    all.equal(
      cl_simp_cpm,
      permutations_to_condition_pr_mat(cl_simp_perms),
      check.attributes = FALSE
    ),
    "character"
  )

  cl_simp_sim_perms <- replicate(5000, randomizr::conduct_ra(cl_simp_ra))

  expect_equal(
    cl_simp_cpm,
    permutations_to_condition_pr_mat(cl_simp_sim_perms),
    tolerance = 0.01
  )

})

test_that("Blocked and Clustered ra", {

  # Blocked and clustered
  dat <- data.frame(
    bl = c("A", "B", "B", "B", "A", "A", "B", "B"),
    cl = c(1, 2, 3, 3, 4, 4, 5, 5)
  )

  bl_cl_ra <- randomizr::declare_ra(clusters = dat$cl, blocks = dat$bl, block_m = c(1, 2))
  bl_cl_perms <- randomizr::obtain_permutation_matrix(bl_cl_ra)

  expect_equal(
    declaration_to_condition_pr_mat(bl_cl_ra),
    permutations_to_condition_pr_mat(bl_cl_perms)
  )
})

test_that("Blocked and clusted ra with remainder", {

  # with remainder
  dat <- data.frame(
    bl = c("A", "B", "B", "B", "A", "A", "B", "B"),
    cl = c(1, 2, 3, 3, 4, 4, 5, 5)
  )

  bl_cl_ra <- randomizr::declare_ra(clusters = dat$cl, blocks = dat$bl, prob = 0.5)
  bl_cl_perms <- randomizr::obtain_permutation_matrix(bl_cl_ra)

  expect_equal(
    declaration_to_condition_pr_mat(bl_cl_ra),
    permutations_to_condition_pr_mat(bl_cl_perms)
  )
})

test_that("Custom ra", {
  cust_perms <- cbind(c(1, 0, 1, 0), c(1, 1, 0, 0))
  cust_ra <- randomizr::declare_ra(permutation_matrix = cust_perms)

  expect_equal(
    declaration_to_condition_pr_mat(cust_ra),
    permutations_to_condition_pr_mat(cust_perms)
  )
})

test_that("Errors for things that we can't support", {

  #
  # multiple armed experiments
  mult_ra <- randomizr::declare_ra(N = 10, prob_each = c(0.2, 0.2, 0.6))
  expect_error(
    declaration_to_condition_pr_mat(mult_ra),
    "`ra_declaration` must have only two arms when passed directly"
  )

  # Permutation error
  expect_error(
    permutations_to_condition_pr_mat(matrix(c(1, 2, 2, 1), nrow = 2)),
    "Matrix of `permutations` must be comprised of only 0s and 1s"
  )

  # Not unique treatment prob for all clusters when complete randomized
  expect_error(
    gen_pr_matrix_cluster(
      c(1, 1, 2, 2),
      treat_probs = runif(4),
      simple = FALSE
    ),
    "Treatment probabilities cannot vary within blocks"
  )
})

test_that("probability not fixed within blocks", {
  bl_small <- randomizr::declare_ra(
    blocks = c(1, 1, 2, 2),
    prob = 0.4
  )
  assign(
    "probabilities_matrix",
    matrix(
      c(0.4, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3),
      ncol = 2,
      dimnames = list(NULL, c("prob_0", "prob_1"))
    ),
    bl_small
  )
  expect_error(
    declaration_to_condition_pr_mat(bl_small),
    "Treatment probabilities must be fixed within blocks for block randomized"
  )
})

test_that("N=2, m=1", {

  comp <- randomizr::declare_ra(N = 2, m = 1)
  assign(
    "probabilities_matrix",
    matrix(
      c(0.4, 0.5, 0.6, 0.5),
      ncol = 2,
      dimnames = list(NULL, c("prob_0", "prob_1"))
    ),
    comp
  )
  expect_error(
    declaration_to_condition_pr_mat(comp),
    "Treatment probabilities must be fixed for complete randomized designs"
  )

  # error in internal function
  expect_error(
    estimatr:::gen_pr_matrix_block(c(1, 2), c(1, 2)),
    "Must specify one of `t`, `p2`, or `p1`"
  )


})
