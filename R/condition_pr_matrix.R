#' Build condition probaability matrix for Horvitz-Thompson estimation from randomization declaration
#'
#' @param declaration An object of class 'ra_declaration' that contains the experimental design
#'
#' @export
declaration_to_condition_pr_mat <- function(declaration) {

  if (class(declaration) != 'ra_declaration') {
    stop("declaration must be an object of class 'ra_declaration'")
  }

  if (declaration$cleaned_arguments$num_arms > 2) {
    stop("The 'declaration' argument can only be used with a binary treatment variable.")
  }

  declaration_call <- as.list(declaration$original_call)
  simple <- eval(declaration_call$simple)

  if (declaration$ra_type == 'simple') {

    mat_00 <- tcrossprod(declaration$probabilities_matrix[, 1])
    diag(mat_00) <- declaration$probabilities_matrix[, 1]
    mat_11 <- tcrossprod(declaration$probabilities_matrix[, 2])
    diag(mat_11) <- declaration$probabilities_matrix[, 2]

    joint_probs <- tcrossprod(declaration$probabilities_matrix[, 1],
                              declaration$probabilities_matrix[, 2])

    mat_01 <- joint_probs
    mat_01 <- copy_upper_to_lower_triangle(mat_01)
    mat_10 <- joint_probs
    mat_10 <- copy_lower_to_upper_triangle(mat_10)
    diag(mat_01) <- diag(mat_10) <- 0

    condition_pr_matrix <-
      rbind(cbind(mat_00, mat_01),
            cbind(mat_10, mat_11))

  } else if (declaration$ra_type == 'complete') {

    condition_pr_matrix <-
      gen_pr_matrix_complete(declaration$probabilities_matrix[, 2])

  } else if (declaration$ra_type == 'clustered') {

    #print(simple)
    condition_pr_matrix <-
      gen_pr_matrix_cluster(
        clusters = declaration$clust_var,
        treat_probs = declaration$probabilities_matrix[, 2],
        simple = simple
      )

  } else if (declaration$ra_type == 'blocked') {
    stop('blocked designs cannot be read from declare_ra for now.')
  }

  return(condition_pr_matrix)

}

#' @export
gen_pr_matrix_cluster <- function(clusters, treat_probs, simple) {

  n <- length(clusters)
  cluster_lists <- split(seq_len(n), clusters)
  n_clust <- length(cluster_lists)

  unique_first_in_cl <- !duplicated(clusters)
  cluster_marginal_probs <-
    treat_probs[unique_first_in_cl]

  # Complete random sampling
  if (is.null(simple) || !simple) {

    ## todo: make work with odd number of clusters
    ## Conditional probabilities
    # p(j==0|i==0)
    pr_j0_given_i0 <-
      (
        (n_clust * (1-cluster_marginal_probs)) # total 0s
        - 1 # remove i == 0
      ) /
      (n_clust - 1) # remaining units

    # p(j==0|i==1)
    pr_j0_given_i1 <-
      (
        (n_clust * (1-cluster_marginal_probs)) # total 0s
      ) /
      (n_clust - 1) # remaining units

    # p(j==1|i==0)
    pr_j1_given_i0 <-
      (
        n_clust * cluster_marginal_probs # total 1s
      ) /
      (n_clust - 1) # remaining units


    # p(j==1|i==1)
    pr_j1_given_i1 <-
      (
        (n_clust * cluster_marginal_probs) # total 1s
        - 1 # remove i == 1
      ) /
      (n_clust - 1) # remaining units

  } else if (simple) { # cluster, simple randomized
    pr_j0_given_i0 <- pr_j0_given_i1 <-
      1 - cluster_marginal_probs

    pr_j1_given_i0 <- pr_j1_given_i1 <-
      cluster_marginal_probs
  }


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
          pr_j1_given_i1[j]

        mat_00[cluster_lists[[i]], cluster_lists[[j]]] <-
          (1 - cluster_marginal_probs[i]) *
          pr_j0_given_i0[j]

        mat_01[cluster_lists[[i]], cluster_lists[[j]]] <-
          (1 - cluster_marginal_probs[i]) *
          pr_j1_given_i0[j]

        mat_10[cluster_lists[[i]], cluster_lists[[j]]] <-
          cluster_marginal_probs[i] *
          pr_j0_given_i1[j]

      }
    }
  }

  condition_pr_matrix <-
    rbind(cbind(mat_00, mat_01),
          cbind(mat_10, mat_11))

  return(condition_pr_matrix)
}

#' Build condition probaability matrix for Horvitz-Thompson estimation from treatment permutations
#'
#' @param permutations A matrix where the rows are units and the columns are different treatment permutations; treated units must be represented with a 1 and control units with a 0
#'
#' @export
permutations_to_condition_pr_mat <- function(permutations) {

  N <- nrow(permutations)

  if (!all(permutations %in% c(0, 1))) {
    stop("Permutations matrix must only have 0s and 1s in it.")
  }

  condition_pr_matrix <- tcrossprod(rbind(permutations, 1 - permutations)) / ncol(permutations)


  colnames(condition_pr_matrix) <- rownames(condition_pr_matrix) <-
    c(paste0("0_", 1:N), paste0("1_", 1:N))

  return(condition_pr_matrix)

}

## Helper functions
copy_upper_to_lower_triangle <- function(mat) {
  mat[lower.tri(mat, diag = F)] <- t(mat)[lower.tri(mat)]
  return(mat)
}

copy_lower_to_upper_triangle <- function(mat) {
  mat[upper.tri(mat, diag = F)] <- t(mat)[upper.tri(mat)]
  return(mat)
}
