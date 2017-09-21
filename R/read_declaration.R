copy_upper_to_lower_triangle <- function(mat) {
  mat[lower.tri(mat, diag = F)] <- t(mat)[lower.tri(mat)]
  return(mat)
}

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

    if (is.null(declaration$cluster_var)) {
      cluster <- NULL
    } else {
      cluster <- declaration$cluster_var
    }

    ## take IPW as weights
    if (length(unique(declaration$probabilities_matrix[, 1])) == 1) {
      weights <- NULL
    } else {
      weights <- 1 / declaration$probabilities_matrix[, 2] ## todo: complement for control units
    }

    if (estimator == 'horvitz_thompson') {
      if (is.null(blocks) & is.null(cluster)) {
        if (declaration$ra_type == 'simple') {
          condition_probabilities <- declaration$probabilities_matrix[, 2]
          condition_probability_matrix <- NULL
        } else {
          condition_probabilities <- declaration$probabilities_matrix[, 2]

          n <- length(condition_probabilities)
          ## todo: only get one of the off-diagonals if symmetric, which it should be
          mat_11 <- mat_00 <- mat_01 <- mat_10 <- matrix(nrow = n, ncol = n)

          ## todo: check if all prs are the same, as they should be
          ## todo: speed by putting in C++ or using the tirangles, or find a matrix way to do it
          for(i in 1:n) {
            for(j in 1:n) {
              mat_11[i,j] <- probs[i] * ((n*probs[j])-1)/(n-1)
              mat_00[i,j] <- (1-probs[i]) * ((n*(1-probs[j]))-1)/(n-1)
              mat_01[i,j] <- (1-probs[i]) * ((n*probs[j]))/(n-1)
              mat_10[i,j] <- probs[i] * ((n*(1-probs[j])))/(n-1)
            }
          }

          diag(mat_11) <- probs
          diag(mat_00) <- (1-probs)
          diag(mat_01) <- 0
          diag(mat_10) <- 0
          condition_probability_matrix <- cbind(rbind(mat_00, mat_01),
                                                rbind(mat_10, mat_11))
        }
      }
    }

    print(condition_probabilities)
    return(
      list(
        blocks = blocks,
        cluster = cluster,
        weights = weights,
        condition_probabilities = condition_probabilities,
        condition_probability_matrix = condition_probability_matrix
      )
    )
  }
