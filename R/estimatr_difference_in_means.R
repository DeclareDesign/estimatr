#' Built-in Estimators: Difference-in-means
#'
#' @param formula An object of class "formula", such as Y ~ Z
#' @param blocks An optional bare (unquoted) name of the block variable. Use for blocked designs only.
#' @param clusters An optional bare (unquoted) name of the variable that corresponds to the clusters in the data; used for cluster randomized designs. For blocked designs, clusters must be within blocks.
#' @param data A data.frame.
#' @param weights An optional bare (unquoted) name of the weights variable.
#' @param subset An optional bare (unquoted) expression specifying a subset of observations to be used.
#' #' @param ci A boolean for whether to compute and return pvalues and confidence intervals, TRUE by default.
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
#'  library(fabricatr)
#'  library(randomizr)
#'  # Get appropriate standard errors for simple designs
#'  dat <- fabricate(
#'    N = 100,
#'    Y = rnorm(100),
#'    Z_simp = simple_ra(N, prob = 0.4),
#'  )
#'
#'  table(dat$Z_simp)
#'  difference_in_means(Y ~ Z_simp, data = dat)
#'
#'  # Accurates estimates and standard errors for clustered designs
#'  dat$clust <- sample(20, size = nrow(dat), replace = TRUE)
#'  dat$Z_clust <- cluster_ra(dat$clust, prob = 0.6)
#'
#'  table(dat$Z_clust, dat$clust)
#'  difference_in_means(Y ~ Z_clust, clusters = clust, data = dat)
#'
#'  # Accurate estimates and standard errors for blocked designs
#'  dat$block <- rep(1:10, each = 10)
#'  dat$Z_block <- block_ra(dat$block, prob = 0.5)
#'
#'  table(dat$Z_block, dat$block)
#'  difference_in_means(Y ~ Z_block, blocks = block, data = dat)
#'
#'  # Matched-pair estimates and standard errors are also accurate
#'  # Specified same as blocked design, function learns that
#'  # it is matched pair from size of blocks!
#'  dat$pairs <- rep(1:50, each = 2)
#'  dat$Z_pairs <- block_ra(dat$pairs, prob = 0.5)
#'
#'  table(dat$pairs, dat$Z_pairs)
#'  difference_in_means(Y ~ Z_pairs, blocks = pairs, data = dat)
#'
#'  # Also works with multi-valued treatments if users specify
#'  # comparison of interest
#'  dat$Z_multi <- simple_ra(
#'    nrow(dat),
#'    condition_names = c("Treatment 2", "Treatment 1", "Control"),
#'    prob_each = c(0.4, 0.4, 0.2)
#'  )
#'
#'  # Only need to specify which condition is treated "condition2" and
#'  # which is control "condition1"
#'  difference_in_means(
#'    Y ~ Z_multi,
#'    condition1 = "Treatment 2",
#'    condition2 = "Control",
#'    data = dat
#'  )
#'  difference_in_means(
#'    Y ~ Z_multi,
#'    condition1 = "Treatment 1",
#'    condition2 = "Control",
#'    data = dat
#'  )
#'
difference_in_means <-
  function(formula,
           blocks,
           clusters,
           condition1 = NULL,
           condition2 = NULL,
           data,
           weights,
           subset,
           ci = TRUE,
           alpha = .05) {

    if (length(all.vars(formula[[3]])) > 1) {
      stop(
        "'formula' must have only one variable on the right-hand side: the ",
        "treatment variable."
      )
    }

    where <- parent.frame()
    model_data <- eval(substitute(
      clean_model_data(
        formula = formula,
        data = data,
        subset = subset,
        cluster = clusters,
        block = blocks,
        weights = weights,
        where = where
      )
    ))

    data <- data.frame(y = model_data$outcome,
                       t = model_data$original_treatment)
    data$cluster <- model_data$cluster
    data$weights <- model_data$weights
    data$block <- model_data$block

    rm(model_data)

    # parse condition names
    if (is.null(condition1) || is.null(condition2)) {
      condition_names <- parse_conditions(
        treatment = data$t,
        condition1 = condition1,
        condition2 = condition2,
        estimator = "difference_in_means"
      )
      condition2 <- condition_names[[2]]
      condition1 <- condition_names[[1]]
    }

    if (is.null(data$block)){

      return_frame <- difference_in_means_internal(
        condition1 = condition1,
        condition2 = condition2,
        data = data,
        alpha = alpha
      )

      # For clustered cases Gerber & Green suggest footnote 20 on ch 3
      # Instead we use CR2 in lm_robust
      # For un-clustered cases we just do the following
      if (is.na(return_frame$df)) {
        return_frame$df <- with(return_frame,
                                N - 2)
      }

    } else {

      pair_matched <- FALSE

      # When learning whether it is matched pairs, should only use relevant conditions
      data <- subset.data.frame(data, t %in% c(condition1, condition2))

      clust_per_block <- check_clusters_blocks(data)

      # Check if design is pair matched
      if (any(clust_per_block == 1)) {
        stop("All blocks must have multiple units (or clusters)")
      } else if (all(clust_per_block == 2)) {
        pair_matched <- TRUE
      } else if (any(clust_per_block == 2) & any(clust_per_block > 2)) {
        stop(
          "Blocks must either all have two units (or clusters) or all have ",
          "more than two units. You cannot mix blocks of size two with ",
          "blocks of a larger size."
        )
      }

      block_dfs <- split(data, data$block)

      block_estimates <- lapply(block_dfs, function(x) {
        difference_in_means_internal(
          data = x,
          condition1 = condition1,
          condition2 = condition2,
          pair_matched = pair_matched,
          alpha = alpha
        )
      })

      block_estimates <- do.call(rbind, block_estimates)

      N_overall <- with(block_estimates, sum(N))

      # Blocked design, (Gerber Green 2012, p73, eq3.10)
      diff <- with(block_estimates, sum(est * N/N_overall))

      df <- NA
      n_blocks <- nrow(block_estimates)

      if (pair_matched) {

        if (is.null(data$cluster)) {
          # Pair matched, unit randomized (Gerber Green 2012, p77, eq3.16)
          se <-
            with(
              block_estimates,
              sqrt( (1 / (n_blocks * (n_blocks - 1))) * sum((est - diff)^2) )
            )

        } else {
          # Pair matched, cluster randomized (Imai, King, Nall 2009, p36, eq6)
          se <-
            with(
              block_estimates,
              sqrt(
                ( n_blocks / ((n_blocks - 1) * N_overall^2) ) *
                  sum( (N * est - (N_overall * diff)/n_blocks)^2 )
              )
            )
        }

        # For pair matched, cluster randomized Imai et al. 2009 recommend (p. 37)
        df <- n_blocks - 1

      } else {
        # Block randomized (Gerber and Green 2012, p. 74, footnote 17)
        se <- with(block_estimates, sqrt(sum(se^2 * (N/N_overall)^2)))


        ## we don't know if this is correct!
        ## matches lm_lin, two estimates per block
        if (is.null(data$cluster)) {
          df <- N_overall - 2 * n_blocks
        } else {
          # Also matches lm_lin for even sized clusters, should be conservative
          df <- sum(clust_per_block) - 2 * n_blocks
        }

      }

      return_frame <- data.frame(
        est = diff,
        se = se,
        df = df,
        N = N_overall
      )

    }

    return_list <- add_cis_pvals(return_frame, alpha, ci)

    #print(c("Pair Matched? ", pair_matched))

    return_list <- dim_like_return(return_list,
                                   alpha = alpha,
                                   formula = formula,
                                   conditions = list(condition1, condition2))

    attr(return_list, "class") <- "difference_in_means"

    return(return_list)

  }



weighted_var_internal <- function(w, x, xWbar){
  wbar <- mean(w)
  n <- length(w)
  return(
    n / ((n - 1) * sum(w) ^ 2) *
      (
        sum((w * x - wbar * xWbar) ^ 2)
        - 2 * xWbar * sum((w - wbar) * (w * x - wbar * xWbar))
        + xWbar ^ 2 * sum((w - wbar) ^ 2)
      )
  )
}



difference_in_means_internal <-
  function(condition1 = NULL,
           condition2 = NULL,
           data,
           pair_matched = FALSE,
           alpha = .05) {

    # Check that treatment status is uniform within cluster, checked here
    # so that the treatment vector t doesn't have to be built anywhere else
    if (!is.null(data$cluster)) {

      if (is.factor(data$cluster)) {
        data$cluster <- droplevels(data$cluster)
      }

      if (!all(tapply(data$t, data$cluster, function(x) all(x == x[1])))) {
        stop(
          "All units within a cluster must have the same treatment condition."
        )
      }
    }

    Y2 <- data$y[data$t == condition2]
    Y1 <- data$y[data$t == condition1]

    N2 <- length(Y2)
    N1 <- length(Y1)
    N <- N2 + N1

    if ((N1 == 0) || (N2 == 0)) {
      stop("Must have units with both treatment conditions within each block.")
    }

    ## Check to make sure multiple in each group if pair matched is false
    if (!pair_matched & (N2 == 1 | N1 == 1)) {
      stop(
        "Not a pair matched design and one treatment condition only has one ",
        "value, making standard errors impossible to calculate."
      )
    }

    df <- NA

    if (!is.null(data$cluster) && !pair_matched) {

      # For now, all clustered cases go to lm_robust
      # CR2 nests Gerber and Green 2012, p. 83, eq. 3.23 when clusters are
      # equal sizes (we think) and is more appropriate when clusters are different sizes

      X <- cbind(1, t = as.numeric(data$t == condition2))

      # print("Using lm_robust")
      # TODO currently lm_robust_fit does too much, need to refactor it
      # if it will be used here in the long run
      cr2_out <- lm_robust_fit(
        y = data$y,
        X = cbind(1, t = as.numeric(data$t == condition2)),
        cluster = data$cluster,
        se_type = "CR2",
        weights = data$weights,
        ci = TRUE,
        coefficient_name = "t",
        try_cholesky = TRUE,
        alpha = alpha,
        return_vcov = FALSE
      )

      diff <- cr2_out$est[2]
      se <- cr2_out$se[2]
      df <- cr2_out$df[2]

    } else {

      if (is.null(data$weights)) {

        diff <- mean(Y2) - mean(Y1)

        if (pair_matched) {
          # Pair matched designs
          se <- NA
        } else {
          # Non-pair matched designs, unit level randomization
          var_Y2 <- var(Y2)
          var_Y1 <- var(Y1)

          se <- sqrt(var_Y2 / N2 + var_Y1 / N1)

          df <- se^4 /
            (
              (var_Y2 / N2)^2 / (N2 - 1) +
                (var_Y1 / N1)^2 / (N1 - 1)
            )

        }

      } else {

        # TODO: weights and matched pair
        if (pair_matched) {
          stop("Weights not supported with matched-pairs.")
        }

        w2 <- data$weights[data$t == condition2]
        w1 <- data$weights[data$t == condition1]

        mean2 <- weighted.mean(Y2, w2)
        mean1 <- weighted.mean(Y1, w1)
        diff <-  mean2 - mean1

        var2 <- weighted_var_internal(w2, Y2, mean2)
        var1 <- weighted_var_internal(w1, Y1, mean1)

        se <- sqrt(var2 + var1)

        # todo: check welch approximation with weights
        df <- se^4 /
          (
            (var2^2 / (N2-1)) +
              (var1^2 / (N1-1))
          )
      }

    }

    return_frame <-
      data.frame(
        est = diff,
        se = se,
        N = N,
        df = df
      )

    return(return_frame)

  }
