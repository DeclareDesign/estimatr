#' Horvitz-Thompson estimator of treatment effects
#'
#' @param formula An object of class "formula", such as Y ~ Z
#' @param declaration An object of class "ra_declaration", from the randomizr package
#' @param condition_pr_variable_name A bare (unquoted) name of the variable with the condition 2 (treatment) probabilities
#' @param condition_pr_matrix A 2n * 2n matrix of marginal and joint probabilities of all units in condition1 and condition2; see details
#' @param block_variable_name An optional bare (unquoted) name of the block variable. Use for blocked designs only.
#' @param cluster_variable_name An optional bare (unquoted) name of the variable that corresponds to the clusters in the data; used for cluster randomized designs. For blocked designs, clusters must be within blocks.
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
           declaration = NULL,
           block_variable_name = NULL,
           cluster_variable_name = NULL,
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
      declare_out <- read_declaration(declaration, 'horvitz_thompson')
      condition_probabilities <- declare_out$condition_probabilities
      condition_pr_matrix <- declare_out$condition_pr_matrix
      clusters <- declare_out$clusters
      blocks <- declare_out$blocks

      if (length(blocks) != nrow(data)) {
        stop("the variables in your declaration are of different length than your data")
      }

    }

    ## Get subset of data
    condition_call <- substitute(subset)

    if (!is.null(condition_call)){
      r <- eval(condition_call, data)
      data <- data[r, ]
    }

    missing_data <- !complete.cases(data[, all.vars(formula)])

    if (!is.null(declaration)) {
      ## todo: if this subsetting is slow, could check to see if any(subset) is needed first
      condition_probabilities <- condition_probabilities[r]
      condition_pr_matrix <- condition_pr_matrix[rep(r, 2), rep(r, 2)]
      clusters <- clusters[r]
      blocks <- blocks[r]

      blocks_missing <- find_warn_missing(blocks[!missing_data], "blocking")
      clusters_missing <- find_warn_missing(clusters[!missing_data], "cluster")

      blocks <- blocks[!missing_data][!blocks_missing & !clusters_missing]
      clusters <- clusters[!missing_data][!blocks_missing & !clusters_missing]

      data <- data[!missing_data, , drop = F][!blocks_missing & !clusters_missing, drop = F]

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
        cluster <- droplevels(as.factor(eval(substitute(cluster_variable_name), data)))

        cluster_missing <- find_warn_missing(cluster, "cluster")

        cluster_variable_name <- deparse(substitute(cluster_variable_name))

      } else {
        cluster <- NULL
      }

      ## remove missingness
      if (is.null(cluster) & is.null(blocks)) {
        any_missing <- cond_prob_missing
      } else if (!is.null(cluster) & is.null(blocks)) {
        any_missing <- cond_prob_missing | cluster_missing
        cluster <- cluster[!any_missing]
      } else if (is.null(cluster) & !is.null(blocks)) {
        any_missing <- cond_prob_missing | blocks_missing
        blocks <- blocks[!any_missing]
      } else {
        any_missing <- cond_prob_missing | blocks_missing | cluster_missing
        blocks <- blocks[!any_missing]
        cluster <- cluster[!any_missing]
      }

      condition_probabilities <- condition_probabilities[!any_missing]
      data <- data[!any_missing, ]

      if (!is.null(cluster)) {
        # check condition ps constant within cluster
        if(any(!tapply(condition_probabilities, cluster, function(x) all(x == x[1])))) {
          stop("condition probabilities must be constant within cluster.")
        }
      }

    }

    #-----
    # Estimation
    #-----

    if (is.null(blocks)){
      return_frame <- horvitz_thompson_internal(
        formula,
        condition_probabilities = condition_probabilities,
        condition_pr_matrix = condition_pr_matrix,
        condition1 = condition1,
        condition2 = condition2,
        data = data,
        weights = weights,
        cluster = cluster,
        constant_effects = constant_effects,
        alpha = alpha
      )

      return_frame$df <- with(return_frame,
                              N - 2)
      return_frame$p <- with(return_frame,
                             2 * pt(abs(est / se), df = df, lower.tail = FALSE))
      return_frame$ci_lower <- with(return_frame,
                                    est - qt(1 - alpha / 2, df = df) * se)
      return_frame$ci_upper <- with(return_frame,
                                    est + qt(1 - alpha / 2, df = df) * se)

    } else {

      if (!is.null(cluster)) {

        ## Check that clusters nest within blocks
        if (!all(tapply(blocks, cluster, function(x)
          all(x == x[1])))) {
          stop("All units within a cluster must be in the same block.")
        }

        ## get number of clusters per block
        clust_per_block <- tapply(cluster,
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
          cluster_variable_name = cluster_variable_name,
          weights = weights,
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

      return_frame <- data.frame(
        est = diff,
        se = se,
        p = p,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        df = df,
        stringsAsFactors = FALSE
      )

    }

    #-----
    # Build and return output
    #-----

    return_frame$coefficient_name <- all.vars(formula[[3]])

    rownames(return_frame) <- NULL

    return(
      return_frame[
        ,
        c("coefficient_name", "est", "se", "p", "ci_lower", "ci_upper", "df")
      ]
    )

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
           weights = NULL,
           cluster = NULL,
           cluster_variable_name = NULL,
           pair_matched = FALSE,
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
      cluster <- droplevels(as.factor(data[, cluster_variable_name]))
    }

    if (!is.null(condition_probabilities_name)) {
      condition_probabilities <- data[, condition_probabilities_name]
    }

    # Check that treatment status is uniform within cluster, checked here
    # so that the treatment vector t doesn't have to be built anywhere else
    if (!is.null(cluster)) {

      if (!all(tapply(t, cluster, function(x) all(x == x[1])))) {
        stop(
          "All units within a cluster must have the same treatment condition."
        )
      }
    }

    N <- length(Y)

    ps2 <- condition_probabilities[t == condition2]
    ps1 <- 1 - condition_probabilities[t == condition1]

    Y2 <- Y[t == condition2] / ps2
    Y1 <- Y[t == condition1] / ps1

    diff <- (sum(Y2) - sum(Y1)) / N
    if (is.null(cluster)) {

      if (is.null(condition_pr_matrix)) {
        # SRS
        if (constant_effects) {
          se <-
            sqrt(
              (var_ht_total_no_cov(c((Y[t==condition1] + diff)/(1-ps1), Y2), c(1 - ps1, ps2)) +
                 var_ht_total_no_cov(c(Y1, (Y[t==condition2] - diff)/(1-ps2)), c(ps1, 1 - ps2)))
              / N^2
            )
        } else {
          se <-
            sqrt(
              (var_ht_total_no_cov(Y2, ps2) / length(Y2)^2 +
                 var_ht_total_no_cov(Y1, ps1) / length(Y1)^2)
            )
        }
      } else {
        # CRS
        if (constant_effects) {
          # rescaled potential outcomes
          y0 <- ifelse(t==condition1,
                           Y / (1-condition_probabilities),
                           (Y - diff) / (1 - condition_probabilities))
          y1 <- ifelse(t==condition2,
                         Y / condition_probabilities,
                         (Y + diff) / condition_probabilities)
          se <-
            sqrt(
              (
                ht_var_total2(
                  y1,
                  condition_pr_matrix[(N+1):(2*N), (N+1):(2*N)]
                ) +
                  ht_var_total2(
                    y0,
                    condition_pr_matrix[1:N, 1:N]
                  ) -
                  2 * ht_covar_total(
                    y0 = y0,
                    y1 = y1,
                    pj = condition_pr_matrix[1:N, (N+1):(2*N)],
                    p00 = condition_pr_matrix[1:N, 1:N],
                    p11 = condition_pr_matrix[(N+1):(2*N), (N+1):(2*N)]
                  )
              ) / N^2
            )
        } else {
          se <-
            sqrt(
              ht_var_total2(
                Y2,
                condition_pr_matrix[(N+1):(2*N), (N+1):(2*N)]
               ) / length(Y2)^2 +
              ht_var_total2(
                Y1,
                condition_pr_matrix[1:N, 1:N]
               ) / length(Y1)^2
            )

        }
      }

    } else {

      if(constant_effects)
        stop("constant effects not currently supported")
      k <- length(unique(cluster))

      se <-
        sqrt(
          ht_var_total_clusters(Y2, ps2, cluster[t == condition2]) / (length(Y2)^2) +
            ht_var_total_clusters(Y1, ps1, cluster[t == condition1]) / (length(Y1)^2)
        )

    }

    return_frame <- data.frame(
      est = diff,
      se = se,
      N = N
    )

    return(return_frame)

  }
