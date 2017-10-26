#' Horvitz-Thompson estimator of treatment effects
#'
#' @param formula An object of class "formula", such as Y ~ Z
#' @param condition_pr_variable_name A bare (unquoted) name of the variable with the condition 2 (treatment) probabilities.
#' @param condition_pr_matrix A 2n * 2n matrix of marginal and joint probabilities of all units in condition1 and condition2; see details.
#' @param block_variable_name An optional bare (unquoted) name of the block variable. Use for blocked designs only.
#' @param cluster_variable_name An optional bare (unquoted) name of the variable that corresponds to the clusters in the data; used for cluster randomized designs. For blocked designs, clusters must be within blocks.
#' @param declaration An object of class "ra_declaration", from the randomizr package. Can substitute for manually specifying clusters, blocks, and condition probabilities.
#' @param data A data.frame.
#' @param subset An optional bare (unquoted) expression specifying a subset of observations to be used.
#' @param alpha The significance level, 0.05 by default.
#' @param condition1 names of the conditions to be compared. Effects are estimated with condition1 as control and condition2 as treatment. If unspecified, condition1 is the "first" condition and condition2 is the "second" according to r defaults.
#' @param condition2 names of the conditions to be compared. Effects are estimated with condition1 as control and condition2 as treatment. If unspecified, condition1 is the "first" condition and condition2 is the "second" according to r defaults.
#' @param constant_effects should constant effects be assumed to estimate the covariance between potential outcomes, default = FALSE.
#'
#' @details This function implements the Horvitz-Thompson estimator for treatment effects.
#'
#' @export
#'
horvitz_thompson <-
  function(formula,
           condition_pr_variable_name = NULL,
           condition_pr_matrix = NULL,
           block_variable_name = NULL,
           cluster_variable_name = NULL,
           declaration = NULL,
           estimator = 'ht',
           constant_effects = FALSE,
           data,
           subset = NULL,
           alpha = .05,
           condition1 = NULL,
           condition2 = NULL) {

    #-----
    # Check arguments
    #-----
    if (length(all.vars(formula[[3]])) > 1) {
      stop(
        "The formula should only include one variable on the right-hand side: the treatment variable."
      )
    }

    #-----
    # Parse arguments, clean data
    #-----
    ## todo: deal with declarations and missingness
    ## todo: warn about colliding arguments
    if (!is.null(declaration)) {
      declare_out <- read_declaration(declaration,
                                      'horvitz_thompson')
      clusters <- declare_out$clusters
      blocks <- declare_out$blocks

      if (!is.null(blocks)) {
        if (length(blocks) != nrow(data)) {
          stop("the block variable in your declaration is of different length than your data")
        }
      }

      if (!is.null(clusters)) {
        if (length(clusters) != nrow(data)) {
          stop("the cluster variable in your declaration is of different length than your data")
        }
      }

      condition_probabilities <- declare_out$condition_probabilities
      condition_pr_matrix <- declare_out$condition_pr_matrix
    }

    ## Get subset of data
    condition_call <- substitute(subset)

    if (!is.null(condition_call)){
      r <- eval(condition_call, data)
      data <- data[r, ]
    }

    missing_data <- !complete.cases(data[, all.vars(formula)])

    if (!is.null(declaration)) {

      if (!is.null(condition_call)) {
        condition_probabilities <- condition_probabilities[r]
        condition_pr_matrix <- condition_pr_matrix[rep(r, 2), rep(r, 2)]
        clusters <- clusters[r]
        blocks <- blocks[r]
      }

      any_missing <- missing_data
      if (!is.null(blocks)) {
        blocks_missing <- find_warn_missing(blocks, "blocking")
      }

      if (!is.null(clusters)) {
        clusters_missing <- find_warn_missing(clusters, "cluster")
      }

      if (!is.null(clusters) & is.null(blocks)) {
        any_missing <- any_missing | clusters_missing
        clusters <- clusters[!any_missing]
      } else if (is.null(clusters) & !is.null(blocks)) {
        any_missing <- any_missing | blocks_missing
        blocks <- blocks[!any_missing]
      } else if (!is.null(clusters) & !is.null(clusters)) {
        any_missing <- any_missing | blocks_missing | clusters_missing
        blocks <- blocks[!any_missing]
        clusters <- clusters[!any_missing]
      }

      data <- data[!any_missing, ]

    } else {

      ## Drop rows with missingness
      data <- data[!missing_data, ]

      ## Get condition probabilities
      if (!is.null(substitute(condition_pr_variable_name))) {
        condition_probabilities <- eval(substitute(condition_pr_variable_name), data)

        if (any(condition_probabilities < 0 | condition_probabilities > 1)) {
          stop(
            "Some condition probabilities are outside of [0, 1]."
          )
        }

        cond_prob_missing <- find_warn_missing(condition_probabilities,
                                               "condition probabilities")

        condition_probabilities_name <- deparse(substitute(condition_pr_variable_name))
      } else {
        condition_probabilities <- NULL
      }

      ## Get block variable
      if (!is.null(substitute(block_variable_name))) {
        blocks <- eval(substitute(block_variable_name), data)
        if (is.factor(blocks)) {
          blocks <- droplevels(blocks)
        }

        blocks_missing <- find_warn_missing(blocks, "blocking")

      } else {
        blocks <- NULL
      }

      ## Get cluster variable
      if (!is.null(substitute(cluster_variable_name))) {
        clusters <- droplevels(as.factor(eval(substitute(cluster_variable_name), data)))

        clusters_missing <- find_warn_missing(clusters, "cluster")

        cluster_variable_name <- deparse(substitute(cluster_variable_name))

      } else {
        clusters <- NULL
      }

      ## remove missingness
      if (is.null(clusters) & is.null(blocks)) {
        any_missing <- cond_prob_missing
      } else if (!is.null(clusters) & is.null(blocks)) {
        any_missing <- cond_prob_missing | clusters_missing
        clusters <- clusters[!any_missing]
      } else if (is.null(clusters) & !is.null(blocks)) {
        any_missing <- cond_prob_missing | blocks_missing
        blocks <- blocks[!any_missing]
      } else {
        any_missing <- cond_prob_missing | blocks_missing | clusters_missing
        blocks <- blocks[!any_missing]
        clusters <- clusters[!any_missing]
      }

      condition_probabilities <- condition_probabilities[!any_missing]
      data <- data[!any_missing, ]

      if (!is.null(clusters) & is.null(declaration)) {
        # check condition ps constant within cluster
        if(any(!tapply(condition_probabilities, clusters, function(x) all(x == x[1])))) {
          stop("condition probabilities must be constant within cluster.")
        }
      }
    }

    #-----
    # Estimation
    #-----

    if (is.null(blocks)){
      return_list <-
        horvitz_thompson_internal(
          formula,
          condition_probabilities = condition_probabilities,
          condition_pr_matrix = condition_pr_matrix,
          condition1 = condition1,
          condition2 = condition2,
          data = data,
          clusters = clusters,
          estimator = estimator,
          constant_effects = constant_effects,
          alpha = alpha
        )

      return_list$df <- with(return_list,
                              N - 2)
      return_list$p <- with(return_list,
                             2 * pt(abs(est / se), df = df, lower.tail = FALSE))
      return_list$ci_lower <- with(return_list,
                                    est - qt(1 - alpha / 2, df = df) * se)
      return_list$ci_upper <- with(return_list,
                                    est + qt(1 - alpha / 2, df = df) * se)

    } else {

      if (!is.null(condition_pr_matrix)) {
        stop("Blocks currently only supported for simple random assignment of units or clusters within blocks.")
      }

      if (!is.null(clusters)) {

        ## Check that clusters nest within blocks
        if (!all(tapply(blocks, clusters, function(x)
          all(x == x[1])))) {
          stop("All units within a cluster must be in the same block.")
        }

        ## get number of clusters per block
        clust_per_block <- tapply(clusters,
                                  blocks,
                                  function(x) length(unique(x)))
      } else {
        clust_per_block <- tabulate(as.factor(blocks))
      }

      ## Check if design is pair matched
      if (any(clust_per_block == 1)) {
        stop(
          "Some blocks have only one unit or cluster. Blocks must have multiple units or clusters."
        )
      } else if (all(clust_per_block == 2)) {
        stop(
          "Cannot compute variance for Horvitz-Thompson estimator with one treated and control unit per block."
        )
      } else if (any(clust_per_block == 2) & any(clust_per_block > 2)) {
        stop(
          "Some blocks have two units or clusters. All blocks must be of size 3 or greater."
        )
      }

      block_dfs <- split(data, blocks)

      block_estimates <- lapply(block_dfs, function(x) {
        horvitz_thompson_internal(
          formula,
          data = x,
          condition1 = condition1,
          condition2 = condition2,
          condition_probabilities_name = condition_probabilities_name,
          condition_pr_matrix = condition_pr_matrix,
          cluster_variable_name = cluster_variable_name,
          alpha = alpha
        )
      })

      block_estimates <- do.call(rbind, block_estimates)

      N_overall <- with(block_estimates, sum(N))

      n_blocks <- nrow(block_estimates)

      diff <- with(block_estimates, sum(est * N/N_overall))

      se <- with(block_estimates, sqrt(sum(se^2 * (N/N_overall)^2)))

      ## we don't know if this is correct!
      df <- n_blocks - 2
      p <- 2 * pt(abs(diff / se), df = df, lower.tail = FALSE)
      ci_lower <- diff - qt(1 - alpha / 2, df = df) * se
      ci_upper <- diff + qt(1 - alpha / 2, df = df) * se

      return_list <-
        list(
          est = diff,
          se = se,
          p = p,
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          df = df
        )

    }

    #-----
    # Build and return output
    #-----

    return_list$alpha <- alpha

    return_list$coefficient_name <- all.vars(formula[[3]])

    attr(return_list, "class") <- "horvitz_thompson"

    return(return_list)
  }

var_ht_total_no_cov <-
  function(y, ps) {
    sum((1 - ps) * ps * y^2)
  }


horvitz_thompson_internal <-
  function(formula,
           condition_probabilities = NULL,
           condition_probabilities_name = NULL,
           condition_pr_matrix = NULL,
           condition1 = NULL,
           condition2 = NULL,
           data,
           clusters = NULL,
           cluster_variable_name = NULL,
           pair_matched = FALSE,
           estimator = 'ht',
           constant_effects = FALSE,
           alpha = .05) {

    Y <- data[, all.vars(formula[[2]])]
    t <- data[, all.vars(formula[[3]])]
    if (is.factor(t)) {
      condition_names <- levels(t)
    } else{
      condition_names <- sort(unique(t))
    }

    if (length(condition_names) == 1) {
      stop("Must have units with both treatment conditions within each block.")
    }

    if (is.null(condition1) & is.null(condition2)) {
      condition1 <- condition_names[1]
      condition2 <- condition_names[2]
    }

    if (!is.null(cluster_variable_name)) {
      clusters <- droplevels(as.factor(data[, cluster_variable_name]))
    }

    if (!is.null(condition_probabilities_name)) {
      condition_probabilities <- data[, condition_probabilities_name]
    }

    # Check that treatment status is uniform within cluster, checked here
    # so that the treatment vector t doesn't have to be built anywhere else
    if (!is.null(clusters)) {

      if (!all(tapply(t, clusters, function(x) all(x == x[1])))) {
        stop(
          "All units within a cluster must have the same treatment condition."
        )
      }
    }

    print(table(condition_pr_matrix))

    N <- length(Y)

    ps2 <- condition_probabilities[t == condition2]
    ps1 <- 1 - condition_probabilities[t == condition1]

    Y2 <- Y[t == condition2] / ps2
    Y1 <- Y[t == condition1] / ps1

    # Trying out estimator from Middleton & Aronow 2015 page 51
    if (!is.null(clusters) & estimator == 'ma') {

      k <- length(unique(clusters))
      c_2 <- clusters[t==condition2]
      c_1 <- clusters[t==condition1]
      k_2 <- length(unique(c_2))
      k_1 <- length(unique(c_1))

      totals_2 <- tapply(Y[t == condition2], c_2, sum)
      totals_1 <- tapply(Y[t == condition1], c_1, sum)

      diff <- (mean(totals_2) - mean(totals_1)) * k / N

      se <-
        sqrt(
          k^2 / N^2 * (
            mean((totals_2 - mean(totals_2))^2) / (k_2 - 1) +
            mean((totals_1 - mean(totals_1))^2) / (k_1 - 1)
          )
        )

    } else {

    diff <- (sum(Y2) - sum(Y1)) / N

    if (is.null(condition_pr_matrix)) {
      # Simple random assignment
      # joint inclusion probabilities = product of marginals
      if (constant_effects) {
        # TODO i
        # rescaled potential outcomes
        y0 <- ifelse(t==condition1,
                     Y / (1-condition_probabilities),
                     (Y - diff) / (1-condition_probabilities))
        y1 <- ifelse(t==condition2,
                     Y / condition_probabilities,
                     (Y + diff) / condition_probabilities)

        # print(var_ht_total_no_cov(y1, condition_probabilities) +
        #       var_ht_total_no_cov(y0, 1 - condition_probabilities))
        #
        # print(-2*sum(c(Y[t == condition2], Y[t == condition1] + diff) * c(Y[t == condition2] - diff, Y[t == condition1])))

        se <-
          sqrt(
            (
              var_ht_total_no_cov(y1, condition_probabilities) +
              var_ht_total_no_cov(y0, 1 - condition_probabilities) +
              2 * sum(c(Y[t == condition2], Y[t == condition1] + diff) * c(Y[t == condition2] - diff, Y[t == condition1]))

            )
          ) / N
      } else {
        se <-
          sqrt(
            var_ht_total_no_cov(Y2, ps2) / length(Y2)^2 +
            var_ht_total_no_cov(Y1, ps1) / length(Y1)^2
          )
      }
    } else {
      # Complete random assignment
      if (constant_effects) {
        # rescaled potential outcomes
        y0 <- ifelse(t==condition1,
                     Y / (1-condition_probabilities),
                     (Y - diff) / (1-condition_probabilities))
        y1 <- ifelse(t==condition2,
                     Y / condition_probabilities,
                     (Y + diff) / condition_probabilities)
        se <-
          sqrt(
            ht_var_total(
              c(-y0, y1),
              condition_pr_matrix
            )
          ) / N
      } else {
        # No constant effects assumption
        se <-
          sqrt(
            ht_var_partial(
              Y2,
              condition_pr_matrix[(N + which(t == condition2)), (N + which(t == condition2))]
            ) +
            ht_var_partial(
              Y1,
              condition_pr_matrix[which(t == condition1), which(t == condition1)]
            ) -
            2 * ht_covar_partial(Y2,
                                 Y1,
                                 condition_pr_matrix[(N + which(t == condition2)), which(t == condition1)],
                                 ps2,
                                 ps1) +
            sum(Y[t == condition2]^2/ps2) +
            sum(Y[t == condition1]^2/ps1)
          ) / N
      }
      }
    }

    return_list <-
      list(
        est = diff,
        se = se,
        N = N
      )

    return(return_list)

  }
