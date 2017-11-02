#' Horvitz-Thompson estimator of treatment effects
#'
#' @param formula An object of class "formula", such as Y ~ Z
#' @param condition_pr_variable_name A bare (unquoted) name of the variable with the condition 2 (treatment) probabilities.
#' @param block_variable_name An optional bare (unquoted) name of the block variable. Use for blocked designs only.
#' @param cluster_variable_name An optional bare (unquoted) name of the variable that corresponds to the clusters in the data; used for cluster randomized designs. For blocked designs, clusters must be within blocks.
#' @param condition_pr_matrix A 2n * 2n matrix of marginal and joint probabilities of all units in condition1 and condition2
#' @param declaration An object of class "ra_declaration", from the randomizr package that is an alternative way of specifying the design. Cannot be used along with any of \code{condition_pr_variable_name}, \code{block_variable_name}, \code{cluster_variable_name}, or \code{condition_pr_matrix}.
#' @param data A data.frame.
#' @param subset An optional bare (unquoted) expression specifying a subset of observations to be used.
#' @param alpha The significance level, 0.05 by default.
#' @param condition1 names of the conditions to be compared. Effects are estimated with condition1 as control and condition2 as treatment. If unspecified, condition1 is the "first" condition and condition2 is the "second" according to r defaults.
#' @param condition2 names of the conditions to be compared. Effects are estimated with condition1 as control and condition2 as treatment. If unspecified, condition1 is the "first" condition and condition2 is the "second" according to r defaults.
#' @param se_type can be one of \code{c("youngs", "constant")} and correspond's to estimating the standard errors using Young's inequality (default, conservative), or the constant effects assumption.
#'
#' @details This function implements the Horvitz-Thompson estimator for treatment effects.
#'
#' @export
#'
horvitz_thompson <-
  function(formula,
           condition_pr_variable_name = NULL,
           block_variable_name = NULL,
           cluster_variable_name = NULL,
           condition_pr_matrix = NULL,
           declaration = NULL,
           se_type = c('youngs', 'constant'),
           data,
           subset = NULL,
           alpha = .05,
           condition1 = NULL,
           condition2 = NULL,
           estimator = 'ht') {

    #-----
    # Check arguments
    #-----
    if (length(all.vars(formula[[3]])) > 1) {
      stop(
        "The formula should only include one variable on the right-hand side: the treatment variable."
      )
    }

    se_type <- match.arg(se_type)

    #-----
    # Parse arguments, clean data
    #-----
    ## User can either use declaration or the arguments, not both!
    if (!is.null(declaration) &
        (!is.null(cluster_variable_name) |
         !is.null(condition_pr_variable_name) |
         !is.null(block_variable_name) |
         !is.null(condition_pr_matrix))) {
      stop("Cannot use declaration with any of cluster_variable_name, condition_pr_variable_name, block_variable_name, condition_pr_matrix.")
    }

    ## Clean data
    where <- parent.frame()
    model_data <- eval(substitute(
      clean_model_data(
        formula = formula,
        data = data,
        subset = subset,
        cluster = cluster_variable_name,
        condition_pr = condition_pr_variable_name,
        block = block_variable_name,
        where = where
      )
    ))

    ## condition_pr_matrix, if supplied, must be same length
    if (!is.null(condition_pr_matrix) && (2*length(model_data$outcome) != nrow(condition_pr_matrix))) {
      stop(sprintf("After cleaning the data, it has %d rows while condition_pr_matrix has %d. condition_pr_matrix should have twice the rows.", length(model_data$outcome), nrow(condition_pr_matrix)))
    }

    data <- data.frame(y = model_data$outcome,
                       t = model_data$design_matrix[, ncol(model_data$design_matrix)])

    if (!is.null(declaration)) {

      ## TODO deal with declarations and missingness
      ## TODO warn about colliding arguments
      declare_out <- read_declaration(declaration,
                                      'horvitz_thompson')

      ## Check that returned model data is same length, error if different
      if (length(declare_out$condition_probabilities) != length(data$y)) {
        stop(sprintf("After cleaning the data, it has %d rows while the declaration has %d. The declaration should be th same length as the cleaned data.", length(data$y), length(declare_out$condition_probabilities)))
      }
      data$clusters <- declare_out$clusters
      data$blocks <- declare_out$blocks
      data$condition_probabilities <- declare_out$condition_probabilities

      condition_pr_matrix <- declare_out$condition_pr_matrix
      rm(declare_out)

      # if (!is.null(blocks)) {
      #   blocks_missing <- find_warn_missing(blocks, "blocking")
      # }
      #
      # if (!is.null(clusters)) {
      #   clusters_missing <- find_warn_missing(clusters, "cluster")
      # }
      #
      # if (!is.null(clusters) & is.null(blocks)) {
      #   any_missing <- any_missing | clusters_missing
      #   clusters <- clusters[!any_missing]
      # } else if (is.null(clusters) & !is.null(blocks)) {
      #   any_missing <- any_missing | blocks_missing
      #   blocks <- blocks[!any_missing]
      # } else if (!is.null(clusters) & !is.null(clusters)) {
      #   any_missing <- any_missing | blocks_missing | clusters_missing
      #   blocks <- blocks[!any_missing]
      #   clusters <- clusters[!any_missing]
      # }
      #
      # data <- data[!any_missing, ]

    } else {
      data$clusters <- model_data$cluster
      data$condition_probabilities <- model_data$condition_pr
      data$blocks <- model_data$block

      rm(model_data)

      if (is.null(data$condition_probabilities)) {
        data$condition_probabilities <- diag(condition_pr_matrix)[(length(data$y)+1):(2*length(data$y))]
      }
    }

    if (!is.null(data$clusters)) {
      # check condition ps constant within cluster
      if(any(!tapply(data$condition_probabilities, data$clusters, function(x) all(x == x[1])))) {
        stop("condition probabilities must be constant within cluster.")
      }
    }

    #-----
    # Estimation
    #-----
    if (is.null(data$blocks)){
      return_list <-
        horvitz_thompson_internal(
          condition_pr_matrix = condition_pr_matrix,
          condition1 = condition1,
          condition2 = condition2,
          data = data,
          estimator = estimator,
          se_type = se_type,
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

      if (!is.null(data$blocks)) {

        ## Check that clusters nest within blocks
        if (!all(tapply(data$blocks, data$blocks, function(x)
          all(x == x[1])))) {
          stop("All units within a cluster must be in the same block.")
        }

        ## get number of clusters per block
        clust_per_block <- tapply(data$blocks,
                                  data$blocks,
                                  function(x) length(unique(x)))
      } else {
        clust_per_block <- tabulate(as.factor(data$blocks))
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

      block_dfs <- split(data, data$blocks)

      block_estimates <- lapply(block_dfs, function(x) {
        horvitz_thompson_internal(
          data = x,
          condition1 = condition1,
          condition2 = condition2,
          condition_pr_matrix = condition_pr_matrix,
          se_type = se_type,
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
  function(condition_pr_matrix = NULL,
           condition1 = NULL,
           condition2 = NULL,
           data,
           pair_matched = FALSE,
           estimator = 'ht',
           se_type,
           alpha = .05) {

    if (is.factor(data$t)) {
      condition_names <- levels(data$t)
    } else{
      condition_names <- sort(unique(data$t))
    }

    if (length(condition_names) == 1) {
      stop("Must have units with both treatment conditions within each block.")
    }

    if (is.null(condition1) & is.null(condition2)) {
      condition1 <- condition_names[1]
      condition2 <- condition_names[2]
    }

    # Check that treatment status is uniform within cluster, checked here
    # so that the treatment vector t doesn't have to be built anywhere else
    if (!is.null(data$clusters)) {

      if (!all(tapply(data$t, data$clusters, function(x) all(x == x[1])))) {
        stop(
          "All units within a cluster must have the same treatment condition."
        )
      }
    }

    # print(table(condition_pr_matrix))

    N <- length(data$y)

    t2 <- which(data$t == condition2)
    t1 <- which(data$t == condition1)

    ps2 <- data$condition_probabilities[t2]
    ps1 <- 1 - data$condition_probabilities[t1]

    Y2 <- data$y[t2] / ps2
    Y1 <- data$y[t1] / ps1

    # Trying out estimator from Middleton & Aronow 2015 page 51
    if (!is.null(data$clusters) & estimator == 'ma') {

      print('ma')
      k <- length(unique(data$clusters))
      c_2 <- data$clusters[t2]
      c_1 <- data$clusters[t1]
      k_2 <- length(unique(c_2))
      k_1 <- length(unique(c_1))

      totals_2 <- tapply(data$y[t2], c_2, sum)
      totals_1 <- tapply(data$y[t1], c_1, sum)

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

      if (!is.null(data$clusters) & estimator == 'am') {

        # print('am')

        k <- length(unique(data$clusters))

        y2_totals <- tapply(data$y[t2], data$clusters[t2], sum)
        y1_totals <- tapply(data$y[t1], data$clusters[t1], sum)

        to_drop <- which(duplicated(data$clusters))
        t2 <- which(data$t[-to_drop] == condition2)
        t1 <- which(data$t[-to_drop] == condition1)

        # reorder totals because tapply will sort on cluster
        y2_totals <- y2_totals[as.character(data$clusters[-to_drop][t2])]
        y1_totals <- y1_totals[as.character(data$clusters[-to_drop][t1])]

        condition_pr_matrix <- condition_pr_matrix[-c(to_drop, N + to_drop), -c(to_drop, N + to_drop)]

        Y2 <- y2_totals / diag(condition_pr_matrix)[k + t2]# for now rescale, with joint pr need squared top alone
        Y1 <- y1_totals / diag(condition_pr_matrix)[t1]

        se <-
          sqrt(
            sum(Y2^2) +
            sum(Y1^2) +
            ht_var_partial(
              Y2,
              condition_pr_matrix[(k + t2), (k + t2)]
            ) +
            ht_var_partial(
              Y1,
              condition_pr_matrix[t1, t1]
            ) -
            2 * ht_covar_partial(Y2,
                                 Y1,
                                 condition_pr_matrix[(k + t2), t1],
                                 ps2,
                                 ps1)
          ) / N

      } else if (is.null(condition_pr_matrix)) {
        # Simple random assignment
        # joint inclusion probabilities = product of marginals
        if (se_type ==  'constant') {
          # rescaled potential outcomes
          y0 <- ifelse(data$t==condition1,
                       data$y / (1-data$condition_probabilities),
                       (data$y - diff) / (1-data$condition_probabilities))
          y1 <- ifelse(data$t==condition2,
                       data$y / data$condition_probabilities,
                       (data$y + diff) / data$condition_probabilities)

          # print(var_ht_total_no_cov(y1, data$condition_probabilities) +
          #       var_ht_total_no_cov(y0, 1 - data$condition_probabilities))
          #
          # print(-2*sum(c(data$y[t2], data$y[t1] + diff) * c(data$y[t2] - diff, data$y[t1])))

          se <-
            sqrt(
              (
                var_ht_total_no_cov(y1, data$condition_probabilities) +
                var_ht_total_no_cov(y0, 1 - data$condition_probabilities) +
                # TODO why is it +2 instead of - (looking at old samii/aronow)
                2 * sum(c(data$y[t2], data$y[t1] + diff) * c(data$y[t2] - diff, data$y[t1]))

              )
            ) / N
        } else {
          se <-
            sqrt(
              sum(Y2^2) + sum(Y1^2)
            ) / N
        }
      } else {
        # Complete random assignment
        if (se_type ==  'constant') {
          # rescaled potential outcomes
          y0 <- ifelse(data$t==condition1,
                       data$y / (1-data$condition_probabilities),
                       (data$y - diff) / (1-data$condition_probabilities))
          y1 <- ifelse(data$t==condition2,
                       data$y / data$condition_probabilities,
                       (data$y + diff) / data$condition_probabilities)
          se <-
            sqrt(
              ht_var_total(
                c(-y0, y1),
                condition_pr_matrix
              )
            ) / N
        } else {
          print('full youngs')
          # Young's inequality
          se <-
            sqrt(
              sum(Y2^2) +
              sum(Y1^2) +
              ht_var_partial(
                Y2,
                condition_pr_matrix[(N + t2), (N + t2)]
              ) +
              ht_var_partial(
                Y1,
                condition_pr_matrix[t1, t1]
              ) -
              2 * ht_covar_partial(Y2,
                                   Y1,
                                   condition_pr_matrix[(N + t2), t1],
                                   ps2,
                                   ps1)
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
