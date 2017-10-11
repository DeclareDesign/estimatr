copy_upper_to_lower_triangle <- function(mat) {
  mat[lower.tri(mat, diag = F)] <- t(mat)[lower.tri(mat)]
  return(mat)
}

## todo: figure out which assignment is the "treatment" in Z, or the binary variable
#' Learns experimental design from randomization declaration
#'
#' @param declaration An object of class 'ra_declaration' that contains the experimental design
#' @param estimator The estimator that needs to learn the experimental design
#'
#' @export
read_declaration <-
  function(
    declaration,
    estimator
  ) {

    if (class(declaration) != 'ra_declaration') {
      stop("declaration must be an object of class 'ra_declaration'")
    }

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
          gen_pr_matrix_complete(declaration$probabilities_matrix[, 2])

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
