#' Built-in Estimators: Difference-in-means
#'
#' @param formula An object of class "formula", such as Y ~ Z
#' @param block_variable_name An optional bare (unquote) name of the block variable. Use for blocked designs only.
#' @param cluster_variable_name An optional bare (unquoted) name of the variable that corresponds to the clusters in the data; used for cluster randomized designs. For blocked designs, clusters must be within blocks.
#' @param data A data.frame.
#' @param weights An optional vector of weights (not yet implemented).
#' @param subset An optional bare (unquoted) expression specifying a subset of observations to be used.
#' @param alpha The significance level, 0.05 by default.
#' @param condition1 names of the conditions to be compared. Effects are estimated with condition1 as control and condition2 as treatment. If unspecified, condition1 is the "first" condition and condition2 is the "second" according to r defaults.
#' @param condition2 names of the conditions to be compared. Effects are estimated with condition1 as control and condition2 as treatment. If unspecified, condition1 is the "first" condition and condition2 is the "second" according to r defaults.
#'
#'
#' @details This function implements difference-in-means estimation, with and without blocking. Standard errors are estimated as the square root of the sum of the within-group variances, divided by their respective sample sizes (Equation 3.6 in Gerber and Green 2012). If blocked, the difference in means estimate is taken in each block, then averaged together according to block size.
#'
#' @export
#'
#' @examples
#'
#'  df <- data.frame(Y = rnorm(100),
#'                   Z = sample(1:3, 100, replace = TRUE),
#'                   block = sample(c("A", "B", "C"), 100, replace = TRUE))
#'
#'  difference_in_means(Y ~ Z, data = df)
#'  difference_in_means(Y ~ Z, condition1 = 3, condition2 = 2, data = df)
#'
#'  difference_in_means(Y ~ Z, block_variable_name = block, data = df)
#'  difference_in_means(Y ~ Z, block_variable_name = block, condition1 = 3, condition2 = 2, data = df)
#'
difference_in_means <-
  function(formula,
           block_variable_name = NULL,
           cluster_variable_name = NULL,
           condition1 = NULL,
           condition2 = NULL,
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

    ## Get block variable
    if (!is.null(substitute(block_variable_name))) {
      blocks <- eval(substitute(block_variable_name), data)
      if (is.factor(blocks)) {
        blocks <- droplevels(blocks)
      }
    } else {
      blocks <- NULL
    }

    ## Get weights variable
    if (!is.null(substitute(weights))) {
      weights <- eval(substitute(weights), data)
    }

    ## Get cluster variable
    if (!is.null(substitute(cluster_variable_name))) {
      cluster <- eval(substitute(cluster_variable_name), data)
      if (is.factor(cluster)) {
        cluster <- droplevels(cluster)
      }

      cluster_variable_name <- deparse(substitute(cluster_variable_name))
    } else {
      cluster <- NULL
    }

    if (is.null(blocks)){

      return_df <- difference_in_means_internal(
        formula,
        condition1 = condition1,
        condition2 = condition2,
        data = data,
        weights = weights,
        cluster = cluster,
        alpha = alpha
      )

      ## todo: add inflation from GG fn 20 ch 3
      return_df$df <- with(return_df, N - 2)
      return_df$p <- with(return_df, 2 * pt(abs(est / se), df = df, lower.tail = FALSE))
      return_df$ci_lower <- with(return_df, est - qt(1 - alpha / 2, df = df) * se)
      return_df$ci_upper <- with(return_df, est + qt(1 - alpha / 2, df = df) * se)

      return_df <- return_df[,c("est", "se", "p", "ci_lower", "ci_upper", "df")]

      return(return_df)

    } else {

      pair_matched = FALSE

      if (!is.null(cluster)) {

        # Check that clusters nest within blocks
        if (!all(tapply(blocks, cluster, function(x)
          all(x == x[1])))) {
          stop("All units within a cluster must be in the same block.")
        }

        # get number of clusters per block
        clust_per_block <- tapply(cluster, blocks, function(x) length(unique(x)))
      } else {
        clust_per_block <- tabulate(as.factor(blocks))
      }

      if (any(clust_per_block == 1)) {
        stop(
          "Some blocks have only one unit or cluster. Blocks must have multiple units or clusters."
        )
      } else if (all(clust_per_block == 2)) {
        pair_matched = TRUE
      } else if (any(clust_per_block == 2) & any(clust_per_block > 2)) {
        stop(
          "Some blocks have two units or clusters, while others have more units or clusters. Design must either be paired or all blocks must be of size 3 or greater."
        )
      }

      block_dfs <- split(data, blocks)

      block_estimates <- lapply(block_dfs, function(x) {
        difference_in_means_internal(
          formula,
          data = x,
          condition1 = condition1,
          condition2 = condition2,
          cluster_variable_name = cluster_variable_name,
          pair_matched = pair_matched,
          weights = weights,
          alpha = alpha
        )
      })

      block_estimates <- do.call(rbind, block_estimates)

      N_overall <- with(block_estimates, sum(N))

      # Blocked design, (Gerber and Green 2012, p. 73, eq. 3.10)
      diff <- with(block_estimates, sum(est * N/N_overall))

      if (pair_matched) {

        n_blocks <- nrow(block_estimates)

        if (is.null(cluster)) {
          # Pair matched, cluster randomized (Gerber and Green 2012, p. 77, eq. 3.16)
          se <-
            with(
              block_estimates,
              sqrt( (1 / (n_blocks * (n_blocks - 1))) * sum((est - diff)^2) )
            )
        } else {
          # Pair matched, unit randomized (Imai, King, Nall 2009, p. 36, eq. 6)
          se <-
            with(
              block_estimates,
              sqrt(
                ( n_blocks / ((n_blocks - 1) * N_overall^2) ) *
                  sum( (N * est - (N_overall * diff)/n_blocks)^2 )
              )
            )
        }

      } else {
        # Block randomized (Gerber and Green 2012, p. 74, footnote 17)
        se <- with(block_estimates, sqrt(sum(se^2 * (N/N_overall)^2)))
      }

      ## we don't know if this is correct!
      df <- N_overall - 2
      p <- 2 * pt(abs(diff / se), df = df, lower.tail = FALSE)
      ci_lower <- diff - qt(1 - alpha / 2, df = df) * se
      ci_upper <- diff + qt(1 - alpha / 2, df = df) * se

      return_df <- data.frame(
        est = diff,
        se = se,
        p = p,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        df = df,
        stringsAsFactors = FALSE
      )

      return(return_df)
    }
  }



weighted_var_internal <- function(w, x, xWbar){
  wbar <- mean(w)
  n <- length(w)
  return(n / ((n - 1) * sum(w) ^ 2) * (sum((w * x - wbar * xWbar) ^ 2) -
                                         2 * xWbar * sum((w - wbar) * (w * x - wbar * xWbar)) + xWbar ^ 2 * sum((w -
                                                                                                                   wbar) ^ 2)))
}



difference_in_means_internal <-
  function(formula,
           condition1 = NULL,
           condition2 = NULL,
           data,
           weights = NULL,
           cluster = NULL,
           cluster_variable_name = NULL,
           pair_matched = FALSE,
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
      cluster <- data[, cluster_variable_name]
      if (is.factor(cluster)) {
        cluster <- droplevels(cluster)
      }
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

    Y2 <- Y[t == condition2]
    Y1 <- Y[t == condition1]

    if (is.null(weights)) {

      diff <- mean(Y2) - mean(Y1)

      if (pair_matched) {
        # Pair matched designs
        se <- NA
      } else if (is.null(cluster)) {
        # Non-pair matched designs, unit level randomization
        se <- sqrt(var(Y2) / length(Y2) + var(Y1) / length(Y1))
      } else {
        # Non-pair matched designs, cluster randomization
        # (Gerber and Green 2012, p. 83, eq. 3.23)
        k <- length(unique(cluster))

        se <- sqrt(
          (var(tapply(Y2, cluster[t == condition2], mean)) * N) /
            (k * length(Y2)) +
          (var(tapply(Y1, cluster[t == condition1], mean)) * N) /
            (k * length(Y1))
        )
      }


    } else {

      # TODO: weights and clusters
      # TODO: weights and matched pair

      w2 <- weights[t == condition2]
      w1 <- weights[t == condition1]

      mean2 <- weighted.mean(Y2, w2)
      mean1 <- weighted.mean(Y1, w1)
      diff <-  mean2 - mean1

      se <- sqrt(weighted_var_internal(w2, Y2, mean2) + weighted_var_internal(w1, Y1, mean1))
    }

    return_df <- data.frame(
      est = diff,
      se = se,
      N = N
    )

    return(return_df)

  }
