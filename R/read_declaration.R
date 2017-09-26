copy_upper_to_lower_triangle <- function(mat) {
  mat[lower.tri(mat, diag = F)] <- t(mat)[lower.tri(mat)]
  return(mat)
}

## todo: figure out which assignment is the "treatment" in Z, or the binary variable
#' @export
read_declaration <-
  function(
    declaration,
    estimator
  ) {

    if (declaration$cleaned_arguments$num_arms > 2) {
      stop("The 'declaration' argument can only be used with a binary treatment variable.")
    }

    if (is.null(declaration$block_var)) {
      blocks <- NULL
    } else {
      blocks <- declaration$block_var
    }

    if (is.null(declaration$clust_var)) {
      clusters <- NULL
    } else {
      clusters <- declaration$clust_var
    }

    ## take IPW as weights
    if (length(unique(declaration$probabilities_matrix[, 1])) == 1) {
      weights <- NULL
    } else {
      weights <- 1 / declaration$probabilities_matrix[, 2] ## todo: complement for control units
    }

    if (estimator == 'horvitz_thompson') {
      condition_probabilities <- declaration$probabilities_matrix[, 2]
      if (declaration$ra_type == 'simple') {

        condition_pr_matrix <- NULL

      } else if (declaration$ra_type == 'complete') {
        condition_pr_matrix <-
          gen_pr_mat_complete(declaration$probabilities_matrix[, 2])

       } else if (declaration$ra_type == 'clustered') {

         n <- nrow(declaration$probabilities_matrix)
         cluster_lists <- split(seq_len(n), declaration$clust_var)
         n_clust <- length(cluster_lists)

         unique_first_in_cl <- !duplicated(declaration$clust_var)
         cluster_marginal_probs <- declaration$probabilities_matrix[unique_first_in_cl, 2]

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

         ## container mats

         mat_00 <- mat_01 <- mat_10 <- mat_11 <-
           matrix(NA, nrow = n, ncol = n)

         for(i in seq_along(cluster_lists)) {
           for(j in seq_along(cluster_lists)) {
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

         condition_pr_matrix <- rbind(cbind(mat_00, mat_01),
                                      cbind(mat_10, mat_11))
      }
    }
    return(
      list(
        blocks = blocks,
        clusters = clusters,
        weights = weights,
        condition_probabilities = condition_probabilities,
        condition_pr_matrix = condition_pr_matrix
      )
    )
  }

#' @export
get_pr_matrix_complete <-
  function(condition_probabilities) {
    n <- length(condition_probabilities)
    ## todo: only get one of the off-diagonals if symmetric, which it should be
    mat_11 <- mat_00 <- mat_01 <- mat_10 <- matrix(nrow = n, ncol = n)

    ## todo: check if all prs are the same, as they should be
    ## todo: speed by putting in C++ or using the tirangles, or find a matrix way to do it
    for(i in 1:n) {
      for(j in 1:n) {
        mat_11[i,j] <- condition_probabilities[i] * ((n*condition_probabilities[j])-1)/(n-1)
        mat_00[i,j] <- (1-condition_probabilities[i]) * ((n*(1-condition_probabilities[j]))-1)/(n-1)
        mat_01[i,j] <- (1-condition_probabilities[i]) * ((n*condition_probabilities[j]))/(n-1)
        mat_10[i,j] <- condition_probabilities[i] * ((n*(1-condition_probabilities[j])))/(n-1)
      }
    }

    diag(mat_11) <- condition_probabilities
    diag(mat_00) <- (1-condition_probabilities)
    diag(mat_01) <- 0
    diag(mat_10) <- 0

    condition_pr_matrix <- cbind(rbind(mat_00, mat_01),
                                 rbind(mat_10, mat_11))

    return(condition_pr_matrix)
  }
