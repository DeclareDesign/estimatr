# This function parses condition names for HT and DiM estimators
parse_conditions <- function(treatment, condition1, condition2, estimator) {
  if (is.factor(treatment)) {
    condition_names <- levels(droplevels(treatment))
  } else {
    condition_names <- sort(unique(treatment))
  }

  if (any(!(c(condition1, condition2) %in% condition_names))) {
    stop("`condition1` and `condition2` must be values found in the treatment")
  }

  n_conditions <- length(condition_names)

  conditions <- list(NULL, NULL)

  if (n_conditions > 2) {
    if (is.null(condition1) || is.null(condition2)) {
      stop(
        "Treatment has > 2 values; must specify both 'condition1' and ",
        "'condition2' or use a treatment with only 2 values."
      )
    } else {
      conditions[1:2] <- c(condition1, condition2)
    }
  } else if (n_conditions == 2) {
    if (is.null(condition1) && is.null(condition2)) {
      conditions[1:2] <- condition_names
    } else if (!is.null(condition2)) {
      conditions[1:2] <- c(setdiff(condition_names, condition2), condition2)
    } else if (!is.null(condition1)) {
      conditions[1:2] <- c(condition1, setdiff(condition_names, condition1))
    } else {
      conditions[1:2] <- c(condition1, condition2)
    }
  } else if (n_conditions == 1) {
    # Allowable for HT estimator
    if (estimator != "horvitz_thompson") {
      stop(
        "Must have more than one value in treatment unless using Horvitz-",
        "Thompson estimator."
      )
    }

    if (is.null(condition1) && is.null(condition2)) {
      conditions[2] <- condition_names
    } else if (!is.null(condition2)) {
      conditions[2] <- condition2
    } else if (!is.null(condition1)) {
      conditions[1] <- condition1
    }
  }

  return(conditions)
}

# This function ensures that blocks and  clusters have been specified correctly
check_clusters_blocks <- function(data) {
  if (!is.null(data$cluster)) {
    one_block_per_clust <-
      tapply(data$block, data$cluster, function(x) all(x == x[1]))

    # Check that clusters nest within blocks
    if (any(!one_block_per_clust)) {
      stop("All clusters must be contained within blocks")
    }

    # get number of clusters per block
    clust_per_block <- tapply(
      data$cluster,
      data$block,
      function(x) length(unique(x))
    )
  } else {
    clust_per_block <- tabulate(as.factor(data$block))
  }

  return(clust_per_block)
}

## todo: figure out which assignment is the "treatment" in Z, or the binary variable
## unused for now!
##' Learns experimental design from randomization declaration
##'
##' @param declaration An object of class 'ra_declaration' that contains the experimental design
##' @param estimator The estimator that needs to learn the experimental design
##'
# read_declaration <- function(declaration,
#                              estimator) {
#
#   if (class(declaration) != 'ra_declaration') {
#     stop("declaration must be an object of class 'ra_declaration'")
#   }
#
#   if (declaration$cleaned_arguments$num_arms > 2) {
#     stop("The 'declaration' argument can only be used with a binary treatment variable.")
#   }
#
#   declaration_call <- as.list(declaration$original_call)
#   simple <- eval(declaration_call$simple)
#
#   ret <- list()
#
#   ## take IPW as weights
#   # if (length(unique(declaration$probabilities_matrix[, 2])) == 1) {
#   #   ret[["weights"]] <- 1 / declaration$probabilities_matrix[, 2] ## todo: complement for control units
#   # }
#
#   if (estimator == 'horvitz_thompson') {
#     ret[["condition_probabilities"]] <-
#       declaration$probabilities_matrix[, 2]
#
#     if (declaration$ra_type == 'complete') {
#       ret[["condition_pr_matrix"]] <-
#         gen_pr_matrix_complete(declaration$probabilities_matrix[, 2])
#
#
#     } else if (declaration$ra_type == 'clustered') {
#
#       n <- nrow(declaration$probabilities_matrix)
#       cluster_lists <- split(seq_len(n), declaration$clust_var)
#       n_clust <- length(cluster_lists)
#
#       unique_first_in_cl <- !duplicated(declaration$clust_var)
#       cluster_marginal_probs <-
#         declaration$probabilities_matrix[unique_first_in_cl, 2]
#
#       # Complete random sampling
#       if (is.null(simple) || !simple) {
#
#         ## todo: make work with odd number of clusters
#         ## Conditional probabilities
#         # p(j==0|i==0)
#         pr_j0_given_i0 <-
#           (
#             (n_clust * (1-cluster_marginal_probs)) # total 0s
#             - 1 # remove i == 0
#           ) /
#           (n_clust - 1) # remaining units
#
#         # p(j==0|i==1)
#         pr_j0_given_i1 <-
#            (
#              (n_clust * (1-cluster_marginal_probs)) # total 0s
#            ) /
#            (n_clust - 1) # remaining units
#
#         # p(j==1|i==0)
#         pr_j1_given_i0 <-
#            (
#              n_clust * cluster_marginal_probs # total 1s
#            ) /
#            (n_clust - 1) # remaining units
#
#
#         # p(j==1|i==1)
#         pr_j1_given_i1 <-
#            (
#              (n_clust * cluster_marginal_probs) # total 1s
#              - 1 # remove i == 1
#            ) /
#            (n_clust - 1) # remaining units
#
#       } else if (simple) { # cluster, simple randomized
#         pr_j0_given_i0 <- pr_j0_given_i1 <-
#           1 - cluster_marginal_probs
#
#         pr_j1_given_i0 <- pr_j1_given_i1 <-
#           cluster_marginal_probs
#       }
#
#
#       # container mats
#       mat_00 <- mat_01 <- mat_10 <- mat_11 <-
#         matrix(NA, nrow = n, ncol = n)
#
#       for (i in seq_along(cluster_lists)) {
#         for (j in seq_along(cluster_lists)) {
#           if (i == j) {
#             mat_11[cluster_lists[[i]], cluster_lists[[j]]] <-
#               cluster_marginal_probs[i]
#
#             mat_00[cluster_lists[[i]], cluster_lists[[j]]] <-
#               1 - cluster_marginal_probs[i]
#
#             mat_01[cluster_lists[[i]], cluster_lists[[j]]] <-
#               0
#
#             mat_10[cluster_lists[[i]], cluster_lists[[j]]] <-
#               0
#
#           } else {
#             mat_11[cluster_lists[[i]], cluster_lists[[j]]] <-
#               cluster_marginal_probs[i] *
#               pr_j1_given_i1[j]
#
#             mat_00[cluster_lists[[i]], cluster_lists[[j]]] <-
#               (1 - cluster_marginal_probs[i]) *
#               pr_j0_given_i0[j]
#
#             mat_01[cluster_lists[[i]], cluster_lists[[j]]] <-
#               (1 - cluster_marginal_probs[i]) *
#               pr_j1_given_i0[j]
#
#             mat_10[cluster_lists[[i]], cluster_lists[[j]]] <-
#               cluster_marginal_probs[i] *
#               pr_j0_given_i1[j]
#
#           }
#         }
#       }
#
#       ret[["condition_pr_matrix"]] <- rbind(cbind(mat_00, mat_01),
#                                             cbind(mat_10, mat_11))
#
#     } else if (declaration$ra_type == 'blocked') {
#       stop('blocked designs cannot be read from declare_ra right now')
#     }
#
#   }
#
#   return(ret)
# }
