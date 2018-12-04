obtain <- function(ra_declaration, condition) {
  if (requireNamespace("randomizr", quietly = TRUE)) {
    randomizr::obtain_condition_probabilities(ra_declaration, condition)
  } else {
    ra_declaration$probability_matrix[, paste0("prob_", condition)]
  }
}

#' Builds condition probability matrices for Horvitz-Thompson estimation from
#' \pkg{randomizr} declaration
#'
#' @param ra_declaration An object of class \code{"ra_declaration"}, generated
#' by the \code{\link[randomizr]{declare_ra}} function in \pkg{randomizr}. This
#' object contains the experimental design that will be represented in a
#' condition probability matrix
#' @param condition1 The name of the first condition, often the control group. If \code{NULL},
#' defaults to first condition in randomizr declaration. Either both \code{condition1}
#' and \code{condition2} have to be specified or both left as \code{NULL}.
#' @param condition2 The name of the second condition, often the treatment group. If \code{NULL},
#' defaults to second condition in randomizr declaration. Either both \code{condition1}
#' and \code{condition2} have to be specified or both left as \code{NULL}.
#' @param prob_matrix An optional probability matrix to override the one in
#' \code{ra_declaration}
#'
#' @details This function takes a \code{"ra_declaration"}, generated
#' by the \code{\link[randomizr]{declare_ra}} function in \pkg{randomizr} and
#' returns a 2n*2n matrix that can be used to fully specify the design for
#' \code{\link{horvitz_thompson}} estimation. This is done by passing this
#' matrix to the \code{condition_pr_mat} argument of
#' \code{\link{horvitz_thompson}}.
#'
#' Currently, this function can learn the condition probability matrix for a
#' wide variety of randomizations: simple, complete, simple clustered, complete
#' clustered, blocked, block-clustered.
#'
#' A condition probability matrix is made up of four submatrices, each of which
#' corresponds to the
#' joint and marginal probability that each observation is in one of the two
#' treatment conditions.
#'
#' The upper-left quadrant is an n*n matrix. On the diagonal is the marginal
#' probability of being in condition 1, often control, for every unit
#' (Pr(Z_i = Condition1) where Z represents the vector of treatment conditions).
#' The off-diagonal elements are the joint probabilities of each unit being in
#' condition 1 with each other unit, Pr(Z_i = Condition1, Z_j = Condition1)
#' where i indexes the rows and j indexes the columns.
#'
#' The upper-right quadrant is also an n*n matrix. On the diagonal is the joint
#' probability of a unit being in condition 1 and condition 2, often the
#' treatment, and thus is always 0. The off-diagonal elements are the joint
#' probability of unit i being in condition 1 and unit j being in condition 2,
#' Pr(Z_i = Condition1, Z_j = Condition2).
#'
#' The lower-left quadrant is also an n*n matrix. On the diagonal is the joint
#' probability of a unit being in condition 1 and condition 2, and thus is
#' always 0. The off-diagonal elements are the joint probability of unit i
#' being in condition 2 and unit j being in condition 1,
#' Pr(Z_i = Condition2, Z_j = Condition1).
#'
#' The lower-right quadrant is an n*n matrix. On the diagonal is the marginal
#' probability of being in condition 2, often treatment, for every unit
#' (Pr(Z_i = Condition2)). The off-diagonal elements are the joint probability
#' of each unit being in condition 2 together,
#' Pr(Z_i = Condition2, Z_j = Condition2).
#'
#' @return a numeric 2n*2n matrix of marginal and joint condition treatment
#' probabilities to be passed to the \code{condition_pr_mat} argument of
#' \code{\link{horvitz_thompson}}. See details.
#'
#' @seealso \code{\link{permutations_to_condition_pr_mat}}
#'
#' @examples
#'
#' # Learn condition probability matrix from complete blocked design
#' library(randomizr)
#' n <- 100
#' dat <- data.frame(
#'   blocks = sample(letters[1:10], size = n, replace = TRUE),
#'   y = rnorm(n)
#' )
#'
#' # Declare complete blocked randomization
#' bl_declaration <- declare_ra(blocks = dat$blocks, prob = 0.4, simple = FALSE)
#' # Get probabilities
#' block_pr_mat <- declaration_to_condition_pr_mat(bl_declaration, 0, 1)
#' # Do randomiztion
#' dat$z <- conduct_ra(bl_declaration)
#'
#' horvitz_thompson(y ~ z, data = dat, condition_pr_mat = block_pr_mat)
#'
#' # When you pass a declaration to horvitz_thompson, this function is called
#'
#' # Equivalent to above call
#' horvitz_thompson(y ~ z, data = dat, ra_declaration = bl_declaration)
#'
#' @export
declaration_to_condition_pr_mat <- function(ra_declaration,
                                            condition1 = NULL,
                                            condition2 = NULL,
                                            prob_matrix = NULL) {
  if (!(inherits(ra_declaration, "ra_declaration"))) {
    stop("`ra_declaration` must be an object of class 'ra_declaration'")
  }

  if (!is.numeric(prob_matrix)) {
    prob_matrix <- ra_declaration$probabilities_matrix
  }

  if (ncol(prob_matrix) > 2) {
    stop(
      "`ra_declaration` must have only two arms when passed directly to ",
      "declaration_to_condition_pr_mat()"
    )
  }

  if (is.null(condition1) && is.null(condition2)) {
    condition1 <- ra_declaration$conditions[1]
    condition2 <- ra_declaration$conditions[2]
  } else if (is.null(condition1) && !is.null(condition2)) {
    stop(
      "Cannot have `condition1 == NULL` and `condition2 != NULL`"
    )
  } else if (!is.null(condition1) && is.null(condition2)) {
    stop(
      "Cannot have `condition2 == NULL` and `condition1 != NULL`"
    )
  }

  p1 <- obtain(
    ra_declaration,
    condition1
  )
  p2 <- obtain(
    ra_declaration,
    condition2
  )

  n <- nrow(prob_matrix)

  if (inherits(ra_declaration, "ra_simple")) {
    v <- c(p1, p2)
    condition_pr_matrix <- tcrossprod(v)
    diag(condition_pr_matrix) <- v
    condition_pr_matrix[cbind(n + 1:n, 1:n)] <- 0
    condition_pr_matrix[cbind(1:n, n + 1:n)] <- 0
  } else if (inherits(ra_declaration, "ra_complete")) {
    if (length(unique(p2)) > 1) {
      stop(
        "Treatment probabilities must be fixed for complete randomized designs"
      )
    }

    condition_pr_matrix <-
      gen_pr_matrix_complete(
        pr = p2[1],
        n_total = n
      )
  } else if (inherits(ra_declaration, "ra_clustered")) {
    condition_pr_matrix <- gen_pr_matrix_cluster(
      clusters = ra_declaration$clusters,
      treat_probs = p2,
      simple = ra_declaration$simple
    )
  } else if (inherits(ra_declaration, "ra_blocked")) {
    condition_pr_matrix <- gen_pr_matrix_block(
      blocks = ra_declaration$blocks,
      clusters = NULL,
      p1 = p1,
      p2 = p2
    )
  } else if (inherits(ra_declaration, "ra_blocked_and_clustered")) {
    condition_pr_matrix <- gen_pr_matrix_block(
      blocks = ra_declaration$blocks,
      clusters = ra_declaration$clusters,
      p1 = p1,
      p2 = p2
    )
  } else if (inherits(ra_declaration, "ra_custom")) {
    # Use permutation matrix
    return(permutations_to_condition_pr_mat(ra_declaration$permutation_matrix))
  }

  # Add names
  colnames(condition_pr_matrix) <- rownames(condition_pr_matrix) <-
    c(paste0(condition1, "_", 1:n), paste0(condition2, "_", 1:n))

  return(condition_pr_matrix)
}

#' Builds condition probability matrices for Horvitz-Thompson estimation from
#' permutation matrix
#'
#' @param permutations A matrix where the rows are units and the columns are
#' different treatment permutations; treated units must be represented with a
#' 1 and control units with a 0
#'
#' @details This function takes a matrix of permutations, for example from
#' the \code{\link[randomizr]{obtain_permutation_matrix}} function in
#' \pkg{randomizr} or through simulation and returns a 2n*2n matrix that can
#' be used to fully specify the design for \code{\link{horvitz_thompson}}
#' estimation. You can read more about these matrices in the documentation for
#' the \code{\link{declaration_to_condition_pr_mat}} function.
#'
#' This is done by passing this matrix to the \code{condition_pr_mat} argument
#' of
#'
#' @seealso  \code{\link[randomizr]{declare_ra}},
#' \code{\link{declaration_to_condition_pr_mat}}
#'
#' @return a numeric 2n*2n matrix of marginal and joint condition treatment
#' probabilities to be passed to the \code{condition_pr_mat} argument of
#' \code{\link{horvitz_thompson}}.
#'
#' @examples
#'
#' # Complete randomization
#' perms <- replicate(1000, sample(rep(0:1, each = 50)))
#' comp_pr_mat <- permutations_to_condition_pr_mat(perms)
#'
#' # Arbitrary randomization
#' possible_treats <- cbind(
#'   c(1, 1, 0, 1, 0, 0, 0, 1, 1, 0),
#'   c(0, 1, 1, 0, 1, 1, 0, 1, 0, 1),
#'   c(1, 0, 1, 1, 1, 1, 1, 0, 0, 0)
#' )
#' arb_pr_mat <- permutations_to_condition_pr_mat(possible_treats)
#' # Simulating a column to be realized treatment
#' z <- possible_treats[, sample(ncol(possible_treats), size = 1)]
#' y <- rnorm(nrow(possible_treats))
#' horvitz_thompson(y ~ z, condition_pr_mat = arb_pr_mat)
#'
#' @export
permutations_to_condition_pr_mat <- function(permutations) {
  N <- nrow(permutations)

  if (!all(permutations %in% c(0, 1))) {
    stop("Matrix of `permutations` must be comprised of only 0s and 1s")
  }

  condition_pr_matrix <- tcrossprod(rbind(1 - permutations, permutations)) / ncol(permutations)

  colnames(condition_pr_matrix) <- rownames(condition_pr_matrix) <-
    c(paste0("0_", 1:N), paste0("1_", 1:N))

  return(condition_pr_matrix)
}


#' Generate condition probability matrix given clusters and probabilities
#'
#' @param clusters A vector of clusters
#' @param treat_probs A vector of treatment (condition 2) probabilities
#' @param simple A boolean for whether the assignment is a random sample
#' assignment (TRUE, default) or complete random assignment (FALSE)
#'
#' @return a numeric 2n*2n matrix of marginal and joint condition treatment
#' probabilities to be passed to the \code{condition_pr_mat} argument of
#' \code{\link{horvitz_thompson}}.
#'
#' @seealso \code{\link{declaration_to_condition_pr_mat}}
#'
#' @export
gen_pr_matrix_cluster <- function(clusters, treat_probs, simple) {
  n <- length(clusters)
  cluster_lists <- split(1:n, clusters, drop = TRUE)
  n_clust <- length(cluster_lists)

  unique_first_in_cl <- !duplicated(clusters)

  cluster_marginal_probs <-
    treat_probs[unique_first_in_cl]

  # Container mats
  # Get cluster condition_pr_matrices
  # Complete random sampling
  if (is.null(simple) || !simple) {
    if (length(unique(cluster_marginal_probs)) > 1) {
      stop(
        "Treatment probabilities cannot vary within blocks for ",
        "block-clustered randomized designs and cannot vary within the whole ",
        "sample for complete cluster randomized designs"
      )
    }

    prs <- gen_joint_pr_complete(cluster_marginal_probs[1], n_clust)

    # This definitely could be optimized
    mat_00 <- matrix(prs[["00"]], n, n)
    mat_10 <- matrix(prs[["10"]], n, n)
    mat_11 <- matrix(prs[["11"]], n, n)

    for (i in 1:n_clust) {
      mat_11[cluster_lists[[i]], cluster_lists[[i]]] <-
        cluster_marginal_probs[i]

      mat_00[cluster_lists[[i]], cluster_lists[[i]]] <-
        1 - cluster_marginal_probs[i]

      mat_10[cluster_lists[[i]], cluster_lists[[i]]] <-
        0
    }

    condition_pr_matrix <-
      rbind(
        cbind(mat_00, mat_10),
        cbind(mat_10, mat_11)
      )
  } else if (simple) { # cluster, simple randomized

    # container mats
    mat_00 <- mat_01 <- mat_10 <- mat_11 <-
      matrix(NA, nrow = n, ncol = n)
    for (i in seq_along(cluster_lists)) {
      for (j in seq_along(cluster_lists)) {
        if (i == j) {
          mat_11[cluster_lists[[i]], cluster_lists[[j]]] <-
            cluster_marginal_probs[i]

          mat_00[cluster_lists[[i]], cluster_lists[[j]]] <-
            1 - cluster_marginal_probs[i]

          mat_01[cluster_lists[[i]], cluster_lists[[j]]] <-
            0

          mat_10[cluster_lists[[i]], cluster_lists[[j]]] <-
            0
        } else {
          mat_11[cluster_lists[[i]], cluster_lists[[j]]] <-
            cluster_marginal_probs[i] *
              cluster_marginal_probs[j]

          mat_00[cluster_lists[[i]], cluster_lists[[j]]] <-
            (1 - cluster_marginal_probs[i]) *
              (1 - cluster_marginal_probs[j])

          mat_01[cluster_lists[[i]], cluster_lists[[j]]] <-
            (1 - cluster_marginal_probs[i]) *
              cluster_marginal_probs[j]

          mat_10[cluster_lists[[i]], cluster_lists[[j]]] <-
            cluster_marginal_probs[i] *
              (1 - cluster_marginal_probs[j])
        }
      }
    }

    condition_pr_matrix <-
      rbind(
        cbind(mat_00, mat_01),
        cbind(mat_10, mat_11)
      )
  }

  return(condition_pr_matrix)
}


gen_pr_matrix_block <- function(blocks,
                                clusters,
                                p2 = NULL,
                                p1 = NULL,
                                t = NULL,
                                condition2 = NULL) {
  n <- length(blocks)
  # Assume complete randomization
  condition_pr_matrix <- matrix(NA, nrow = 2 * n, ncol = 2 * n)

  # Split by block and get complete randomized values within each block
  id_dat <- data.frame(ids = 1:n, stringsAsFactors = FALSE)
  if (!is.null(p2)) {
    id_dat$p2 <- p2
  }
  if (!is.null(p1)) {
    id_dat$p1 <- p1
  }
  if (!is.null(t)) {
    id_dat$t <- t
  }

  if (is.null(t) && is.null(p2) && is.null(p1)) {
    stop("Must specify one of `t`, `p2`, or `p1`")
  }

  clustered <- !is.null(clusters)
  if (clustered) {
    id_dat$clusters <- clusters
  }

  block_dat <- split(
    id_dat,
    blocks,
    drop = TRUE
  )

  n_blocks <- length(block_dat)

  for (i in seq_along(block_dat)) {
    ids <- c(block_dat[[i]]$ids, n + block_dat[[i]]$ids)

    if (clustered) {
      if (is.null(block_dat[[i]]$p2)) {
        # learn prs
        cluster_treats <- get_cluster_treats(block_dat[[i]], condition2)
        block_dat[[i]]$p2 <- mean(cluster_treats$treat_clust)
      }

      if (is.null(block_dat[[i]]$p1)) {
        block_dat[[i]]$p1 <- 1 - block_dat[[i]]$p2
      }

      # Has to be complete randomization of clusters
      condition_pr_matrix[ids, ids] <-
        gen_pr_matrix_cluster(
          clusters = block_dat[[i]]$clusters,
          treat_probs = block_dat[[i]]$p2,
          simple = FALSE
        )
    } else {
      if (length(unique(block_dat[[i]]$p2)) > 1) {
        stop(
          "Treatment probabilities must be fixed within blocks for block ",
          "randomized designs"
        )
      }

      if (is.null(block_dat[[i]]$p2)) {
        # learn prs
        block_dat[[i]]$p2 <- mean(block_dat[[i]]$t)
      }

      if (is.null(block_dat[[i]]$p1)) {
        block_dat[[i]]$p1 <- 1 - block_dat[[i]]$p2
      }

      condition_pr_matrix[ids, ids] <-
        gen_pr_matrix_complete(
          pr = block_dat[[i]]$p2[1],
          n_total = length(block_dat[[i]]$p2)
        )
    }
  }

  for (i in seq_along(block_dat)) {
    ids <- c(block_dat[[i]]$ids, n + block_dat[[i]]$ids)

    for (j in seq_along(block_dat)) {
      if (i != j) {
        condition_pr_matrix[
          ids,
          c(block_dat[[j]]$ids, n + block_dat[[j]]$ids)
        ] <- tcrossprod(
          c(block_dat[[i]]$p1, block_dat[[i]]$p2),
          c(block_dat[[j]]$p1, block_dat[[j]]$p2)
        )
      }
    }
  }


  return(condition_pr_matrix)
}

gen_pr_matrix_complete <- function(pr, n_total) {
  prs <- gen_joint_pr_complete(pr, n_total)

  pr00_mat <- matrix(prs[["00"]], nrow = n_total, ncol = n_total)
  diag(pr00_mat) <- 1 - pr
  pr10_mat <- matrix(prs[["10"]], nrow = n_total, ncol = n_total)
  diag(pr10_mat) <- 0
  pr11_mat <- matrix(prs[["11"]], nrow = n_total, ncol = n_total)
  diag(pr11_mat) <- pr

  pr_mat <- cbind(
    rbind(pr00_mat, pr10_mat),
    rbind(pr10_mat, pr11_mat)
  )

  return(pr_mat)
}

gen_joint_pr_complete <- function(pr, n_total) {
  n_treated <- pr * n_total
  remainder <- n_treated %% 1

  n_treated_floor <- floor(n_treated)
  n_control <- n_total - n_treated_floor

  prs <- list()

  prs[["11"]] <-
    remainder * # pr(M)
    ((n_treated_floor + 1) / n_total) * # pr(j = 1 | M)
    (n_treated_floor / (n_total - 1)) + # pr(i = 1 | j = 1, M)
    (1 - remainder) * # pr(M')
      (n_treated_floor / n_total) * # pr(j = 1 | M')
      ((n_treated_floor - 1) / (n_total - 1)) # pr(i = 1 | j = 1, M')


  prs[["10"]] <-
    remainder * # pr(M)
    ((n_control - 1) / n_total) * # pr(j = 0 | M)
    ((n_treated_floor + 1) / (n_total - 1)) + # pr(i = 1 | j = 0, M)
    (1 - remainder) * # pr(M')
      (n_control / n_total) * # pr(j = 0 | M')
      (n_treated_floor / (n_total - 1)) # pr(i = 1 | j = 0, M')

  prs[["00"]] <-
    remainder * # pr(M)
    ((n_control - 1) / n_total) * # pr(j = 0 | M)
    ((n_control - 2) / (n_total - 1)) + # pr(i = 0 | j = 0, M)
    (1 - remainder) * # pr(M')
      (n_control / n_total) * # pr(j = 0 | M')
      ((n_control - 1) / (n_total - 1)) # pr(i = 0 | j = 0, M')

  return(prs)
}


get_cluster_treats <- function(data, condition2) {
  cluster_dat <- split(
    data$t,
    data$clusters,
    drop = TRUE
  )

  n_clust <- length(cluster_dat)
  treat_clust <- numeric(n_clust)

  for (i in seq_along(cluster_dat)) {
    if (length(unique(cluster_dat[[i]])) > 1) {
      stop("Treatment condition must be constant within `clusters`")
    }

    treat_clust[i] <- as.numeric(cluster_dat[[i]][1] == condition2)
  }

  return(list(n_clust = n_clust, treat_clust = treat_clust))
}
