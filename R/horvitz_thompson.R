#' Horvitz-Thompson estimator of treatment effects
#'
#' @param formula An object of class "formula", such as Y ~ Z
#' @param inclusion_pr_variable_name A bare (unquoted) name of the variable with the inclusion probabilities (or probabilities of being in condition 2)
#' @param block_variable_name An optional bare (unquoted) name of the block variable. Use for blocked designs only.
#' @param cluster_variable_name An optional bare (unquoted) name of the variable that corresponds to the clusters in the data; used for cluster randomized designs. For blocked designs, clusters must be within blocks.
#' @param data A data.frame.
#' @param weights An optional vector of weights (not yet implemented).
#' @param subset An optional bare (unquoted) expression specifying a subset of observations to be used.
#' @param alpha The significance level, 0.05 by default.
#' @param condition1 names of the conditions to be compared. Effects are estimated with condition1 as control and condition2 as treatment. If unspecified, condition1 is the "first" condition and condition2 is the "second" according to r defaults.
#' @param condition2 names of the conditions to be compared. Effects are estimated with condition1 as control and condition2 as treatment. If unspecified, condition1 is the "first" condition and condition2 is the "second" according to r defaults.
#'
#'
#' @details This function implements the Horvitz-Thompson estimator for treatment effects.
#'
#' @export
#'
horvitz_thompson <-
  function(formula,
           inclusion_pr_variable_name,
           block_variable_name = NULL,
           cluster_variable_name = NULL,
           condition1 = NULL,
           condition2 = NULL,
           constant_effects = FALSE,
           data,
           weights = NULL,
           subset = NULL,
           alpha = .05) {

    if (length(all.vars(formula[[3]])) > 1) {
      stop(
        "The formula should only include one variable on the right-hand side: the treatment variable."
      )
    }

    ## Get subset of data
    condition_call <- substitute(subset)

    if (!is.null(condition_call)){
      r <- eval(condition_call, data)
      data <- data[r, ]
    }

    ## Drop rows with missingness
    missing_data <- !complete.cases(data[, all.vars(formula)])
    data <- data[!missing_data, ]

    ## Get inclusion probabilities
    inclusion_probabilities <- eval(substitute(inclusion_pr_variable_name), data)

    if (any(inclusion_probabilities < 0 | inclusion_probabilities > 1)) {
      stop(
        "Some inclusion probabilities are outside of [0, 1]."
      )
    }

    incl_prob_missing <- find_warn_missing(inclusion_probabilities,
                                           "inclusion probabilities")


    inclusion_probabilities <- inclusion_probabilities[!incl_prob_missing]
    data <- data[!incl_prob_missing, ]

    inclusion_probabilities_name <- deparse(substitute(inclusion_pr_variable_name))

    ## Get block variable
    if (!is.null(substitute(block_variable_name))) {
      blocks <- eval(substitute(block_variable_name), data)
      if (is.factor(blocks)) {
        blocks <- droplevels(blocks)
      }

      blocks_missing <- find_warn_missing(blocks, "blocking")

      blocks <- blocks[!blocks_missing]
      data <- data[!blocks_missing, ]

    } else {
      blocks <- NULL
    }

    ## Get weights variable
    if (!is.null(substitute(weights))) {
      weights <- eval(substitute(weights), data)

      weights_missing <- find_warn_missing(weights, "weights")

      weights <- weights[!weights_missing]
      data <- data[!weights_missing, ]

    }

    ## Get cluster variable
    if (!is.null(substitute(cluster_variable_name))) {
      cluster <- droplevels(as.factor(eval(substitute(cluster_variable_name), data)))

      cluster_missing <- find_warn_missing(cluster, "cluster")

      cluster <- cluster[!cluster_missing]
      data <- data[!cluster_missing, ]

      # check inclusion ps constant within cluster

      if(any(!tapply(inclusion_probabilities, cluster, function(x) all(x == x[1])))) {
        stop("Inclusion probabilities must be constant within cluster.")
      }

      cluster_variable_name <- deparse(substitute(cluster_variable_name))

    } else {
      cluster <- NULL
    }

    if (is.null(blocks)){

      return_frame <- horvitz_thompson_internal(
        formula,
        inclusion_probabilities = inclusion_probabilities,
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
          inclusion_probabilities_name = inclusion_probabilities_name,
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
           inclusion_probabilities = NULL,
           inclusion_probabilities_name = NULL,
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

    if (!is.null(inclusion_probabilities_name)) {
      inclusion_probabilities <- data[, inclusion_probabilities_name]
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

    ps2 <- inclusion_probabilities[t == condition2]
    ps1 <- 1 - inclusion_probabilities[t == condition1]

    Y2 <- Y[t == condition2] / ps2
    Y1 <- Y[t == condition1] / ps1

    if (is.null(weights)) {

      diff <- (sum(Y2) - sum(Y1)) / N

      if (is.null(cluster)) {

        if(constant_effects) {
          se <-
            sqrt(
              (var_ht_total_no_cov(c(Y1 + diff/(1-ps1), Y2), c(1 - ps1, ps2)) +
                 var_ht_total_no_cov(c(Y1, Y2 - diff/(1-ps2)), c(ps1, 1 - ps2)))
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

        if(constant_effects)
          stop("constant effects not currently supported")
        k <- length(unique(cluster))

        se <-
          sqrt(
            ht_var_total_clusters(Y2, ps2, cluster[t == condition2]) / (length(Y2)^2) +
              ht_var_total_clusters(Y1, ps1, cluster[t == condition1]) / (length(Y1)^2)
          )

      }


    } else {
      stop("Other weights not supported for now.")
    }

    return_frame <- data.frame(
      est = diff,
      se = se,
      N = N
    )

    return(return_frame)

  }
