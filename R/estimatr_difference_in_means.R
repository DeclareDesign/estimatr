#' Design-based difference-in-means estimator
#'
#' @description Difference-in-means estimators that selects the appropriate
#' point estimate, standard errors, and degrees of freedom for a variety of
#' designs: unit randomized, cluster randomized, block randomized,
#' block-cluster randomized, matched-pairs, and matched-pair cluster
#' randomized designs
#'
#' @param formula an object of class formula, as in \code{\link{lm}}, such as
#' \code{Y ~ Z} with only one variable on the right-hand side, the treatment.
#' @param data A \code{data.frame}.
#' @param blocks An optional bare (unquoted) name of the block variable. Use
#' for blocked designs only.
#' @param clusters An optional bare (unquoted) name of the variable that
#' corresponds to the clusters in the data; used for cluster randomized
#' designs. For blocked designs, clusters must nest within blocks.
#' @param weights the bare (unquoted) names of the weights variable in the
#' supplied data.
#' @param subset An optional bare (unquoted) expression specifying a subset of
#' observations to be used.
#' @param se_type An optional string that can be one of \code{c("default", "none")}. If "default" (the default), it will use the default standard error estimator for the design, and if "none" then standard errors will not be computed which may speed up run time if only the point estimate is required.
#' @param condition1 value in the treatment vector of the condition
#' to be the control. Effects are
#' estimated with \code{condition1} as the control and \code{condition2} as the
#' treatment. If unspecified, \code{condition1} is the "first" condition and
#' \code{condition2} is the "second" according to levels if the treatment is a
#' factor or according to a sortif it is a numeric or character variable (i.e
#' if unspecified and the treatment is 0s and 1s, \code{condition1} will by
#' default be 0 and \code{condition2} will be 1). See the examples for more.
#' @param condition2 value in the treatment vector of the condition to be the
#' treatment. See \code{condition1}.
#' @param ci logical. Whether to compute and return p-values and
#' confidence intervals, TRUE by default.
#' @param alpha The significance level, 0.05 by default.
#'
#' @details This function implements a difference-in-means estimator, with
#' support for blocked, clustered, matched-pairs, block-clustered, and
#' matched-pair clustered designs. One specifies their design by passing
#' the blocks and clusters in their data and this function chooses which
#' estimator is most appropriate.
#'
#' If you pass only \code{blocks}, if all blocks are of size two, we will
#' infer that the design is a matched-pairs design. If they are all size four
#' or larger, we will infer that it is a regular blocked design. If you pass
#' both \code{blocks} and \code{clusters}, we will similarly
#' infer whether it is a matched-pairs clustered design or a block-clustered
#' design the number of clusters per block. If the user passes only
#' \code{clusters}, we will infer that the design was cluster-randomized. If
#' the user specifies neither the \code{blocks} nor the \code{clusters},
#' a regular Welch's t-test will be performed.
#'
#' Importantly, if the user specifies weights, the estimation is handed off
#' to \code{\link{lm_robust}} with the appropriate robust standard errors
#' as weighted difference-in-means estimators are not implemented here.
#' More details of the about each of the estimators can be found in the
#' \href{https://declaredesign.org/r/estimatr/articles/mathematical-notes.html}{mathematical notes}.
#'
#' @return Returns an object of class \code{"difference_in_means"}.
#'
#' The post-estimation commands functions \code{summary} and \code{\link{tidy}}
#' return results in a \code{data.frame}. To get useful data out of the return,
#' you can use these data frames, you can use the resulting list directly, or
#' you can use the generic accessor functions \code{coef} and
#' \code{confint}.
#'
#' An object of class \code{"difference_in_means"} is a list containing at
#' least the following components:
#'   \item{coefficients}{the estimated difference in means}
#'   \item{std.error}{the estimated standard error}
#'   \item{statistic}{the t-statistic}
#'   \item{df}{the estimated degrees of freedom}
#'   \item{p.value}{the p-value from a two-sided t-test using \code{coefficients}, \code{std.error}, and \code{df}}
#'   \item{conf.low}{the lower bound of the \code{1 - alpha} percent confidence interval}
#'   \item{conf.high}{the upper bound of the \code{1 - alpha} percent confidence interval}
#'   \item{term}{a character vector of coefficient names}
#'   \item{alpha}{the significance level specified by the user}
#'   \item{N}{the number of observations used}
#'   \item{outcome}{the name of the outcome variable}
#'   \item{design}{the name of the design learned from the arguments passed}
#'
#' @seealso \code{\link{lm_lin}}
#'
#' @references
#' Gerber, Alan S, and Donald P Green. 2012. Field Experiments: Design, Analysis, and Interpretation. New York: W.W. Norton.
#'
#' Imai, Kosuke, Gary King, Clayton Nall. 2009. "The Essential Role of Pair Matching in Cluster-Randomized Experiments, with Application to the Mexican Universal Health Insurance Evaluation." Statistical Science 24 (1). Institute of Mathematical Statistics: 29-53. \url{https://doi.org/10.1214/08-STS274}.
#'
#' @examples
#'
#' library(fabricatr)
#' library(randomizr)
#' # Get appropriate standard errors for unit-randomized designs
#'
#' # ----------
#' # 1. Unit randomized
#' # ----------
#' dat <- fabricate(
#'   N = 100,
#'   Y = rnorm(100),
#'   Z_comp = complete_ra(N, prob = 0.4),
#' )
#'
#' table(dat$Z_comp)
#' difference_in_means(Y ~ Z_comp, data = dat)
#'
#' # ----------
#' # 2. Cluster randomized
#' # ----------
#' # Accurates estimates and standard errors for clustered designs
#' dat$clust <- sample(20, size = nrow(dat), replace = TRUE)
#' dat$Z_clust <- cluster_ra(dat$clust, prob = 0.6)
#'
#' table(dat$Z_clust, dat$clust)
#' summary(difference_in_means(Y ~ Z_clust, clusters = clust, data = dat))
#'
#' # ----------
#' # 3. Block randomized
#' # ----------
#' dat$block <- rep(1:10, each = 10)
#' dat$Z_block <- block_ra(dat$block, prob = 0.5)
#'
#' table(dat$Z_block, dat$block)
#' difference_in_means(Y ~ Z_block, blocks = block, data = dat)
#'
#' # ----------
#' # 4. Block cluster randomized
#' # ----------
#' # Learns this design if there are two clusters per block
#' dat$small_clust <- rep(1:50, each = 2)
#' dat$big_blocks <- rep(1:5, each = 10)
#'
#' dat$Z_blcl <- block_and_cluster_ra(
#'   blocks = dat$big_blocks,
#'   clusters = dat$small_clust
#'  )
#'
#' difference_in_means(
#'   Y ~ Z_blcl,
#'   blocks = big_blocks,
#'   clusters = small_clust,
#'   data = dat
#'  )
#'
#' # ----------
#' # 5. Matched-pairs
#' # ----------
#' # Matched-pair estimates and standard errors are also accurate
#' # Specified same as blocked design, function learns that
#' # it is matched pair from size of blocks!
#' dat$pairs <- rep(1:50, each = 2)
#' dat$Z_pairs <- block_ra(dat$pairs, prob = 0.5)
#'
#' table(dat$pairs, dat$Z_pairs)
#' difference_in_means(Y ~ Z_pairs, blocks = pairs, data = dat)
#'
#' # ----------
#' # 6. Matched-pair cluster randomized
#' # ----------
#' # Learns this design if there are two clusters per block
#' dat$small_clust <- rep(1:50, each = 2)
#' dat$cluster_pairs <- rep(1:25, each = 4)
#' table(dat$cluster_pairs, dat$small_clust)
#'
#' dat$Z_mpcl <- block_and_cluster_ra(
#'   blocks = dat$cluster_pairs,
#'   clusters = dat$small_clust
#'  )
#'
#' difference_in_means(
#'   Y ~ Z_mpcl,
#'   blocks = cluster_pairs,
#'   clusters = small_clust,
#'   data = dat
#'  )
#'
#' # ----------
#' # Other examples
#' # ----------
#'
#' # Also works with multi-valued treatments if users specify
#' # comparison of interest
#' dat$Z_multi <- simple_ra(
#'   nrow(dat),
#'   conditions = c("Treatment 2", "Treatment 1", "Control"),
#'   prob_each = c(0.4, 0.4, 0.2)
#' )
#'
#' # Only need to specify which condition is treated `condition2` and
#' # which is control `condition1`
#' difference_in_means(
#'   Y ~ Z_multi,
#'   condition1 = "Treatment 2",
#'   condition2 = "Control",
#'   data = dat
#' )
#' difference_in_means(
#'   Y ~ Z_multi,
#'   condition1 = "Treatment 1",
#'   condition2 = "Control",
#'   data = dat
#' )
#'
#' # Specifying weights will result in estimation via lm_robust()
#' dat$w <- runif(nrow(dat))
#' difference_in_means(Y ~ Z_comp, weights = w, data = dat)
#' lm_robust(Y ~ Z_comp, weights = w, data = dat)
#'
#' @export
difference_in_means <- function(formula,
                                data,
                                blocks,
                                clusters,
                                weights,
                                subset,
                                se_type = c("default", "none"),
                                condition1 = NULL,
                                condition2 = NULL,
                                ci = TRUE,
                                alpha = .05) {
  if (length(all.vars(f_rhs(eval_tidy(formula)))) > 1) {
    stop(
      "'formula' must have only one variable on the right-hand side: the ",
      "treatment variable."
    )
  }

  se_type <- match.arg(se_type)

  datargs <- enquos(
    formula = formula,
    weights = weights,
    subset = subset,
    block = blocks,
    cluster = clusters
  )
  data <- enquo(data)
  model_data <- clean_model_data(data = data, datargs, estimator = "dim")

  data <- data.frame(
    y = model_data$outcome,
    t = model_data$original_treatment,
    stringsAsFactors = FALSE
  )
  data$cluster <- model_data$cluster
  # rescale weights for convenience
  if (is.numeric(model_data$weights)) {
    data$weights <- model_data$weights / mean(model_data$weights)
  }
  data$block <- model_data$block

  if (!is.null(data$weights) && length(unique(data$weights)) == 1
  && is.null(data$cluster) && is.null(data$block)) {
    message(
      "Constant `weights` passed to `difference_in_means` will ",
      "unnecessarily trigger `lm_robust()` and the Welch-Satterthwaite ",
      "approximation will not be used for the degrees of freedom."
    )
  }

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

  # subset data
  data <- subset.data.frame(data, t %in% c(condition1, condition2))

  if (is.null(data$block)) {
    return_frame <- difference_in_means_internal(
      condition1 = condition1,
      condition2 = condition2,
      data = data,
      alpha = alpha,
      se_type = se_type
    )

    if (is.null(data$cluster)) {
      design <- "Standard"
    } else {
      design <- "Clustered"
    }
  } else {
    pair_matched <- FALSE

    # When learning whether design is matched pairs,
    # should only use rows in relevant conditions
    data <- subset.data.frame(data, t %in% c(condition1, condition2))

    clust_per_block <- check_clusters_blocks(data)

    # Check if design is pair matched
    if (any(clust_per_block == 1)) {
      stop("All `blocks` must have multiple units (or `clusters`)")
    } else if (all(clust_per_block == 2)) {
      pair_matched <- TRUE
    } else if (any(clust_per_block == 2) & any(clust_per_block > 2)) {
      pair_matched <- TRUE
      warning(
        "Some `blocks` have two units/`clusters` while other blocks ",
        "have more units/`clusters`. As standard variance estimates ",
        "cannot be computed within blocks with two units, we use the ",
        "matched pairs estimator of the variance."
      )
    }

    block_dfs <- split(data, data$block)

    block_estimates <- lapply(block_dfs, function(x) {
      difference_in_means_internal(
        data = x,
        condition1 = condition1,
        condition2 = condition2,
        pair_matched = pair_matched,
        alpha = alpha,
        se_type = se_type
      )
    })

    block_estimates <- do.call(rbind, block_estimates)

    N_overall <- with(block_estimates, sum(N))

    # Blocked design, (Gerber Green 2012, p73, eq3.10)
    diff <- with(block_estimates, sum(coefficients * N / N_overall))

    df <- NA
    std.error <- NA
    n_blocks <- nrow(block_estimates)

    if (pair_matched) {
      if (is.null(data$cluster)) {
        design <- "Matched-pair"

        # Pair matched, unit randomized (Gerber Green 2012, p77, eq3.16)
        if (se_type != "none") {
          std.error <-
            with(
              block_estimates,
              sqrt(
                (1 / (n_blocks * (n_blocks - 1))) *
                  sum((coefficients - diff)^2)
              )
            )
        }
      } else {
        design <- "Matched-pair clustered"
        # Pair matched, cluster randomized (Imai, King, Nall 2009, p36, eq6)
        if (se_type != "none") {
          std.error <-
            with(
              block_estimates,
              sqrt(
                (n_blocks / ((n_blocks - 1) * N_overall^2)) *
                  sum((N * coefficients - (N_overall * diff) / n_blocks)^2)
              )
            )
        }
      }

      # For pair matched, cluster randomized Imai et al. 2009 recommend (p. 37)
      df <- n_blocks - 1
    } else {
      # Block randomized (Gerber and Green 2012, p. 74, footnote 17)
      if (se_type != "none") {
        std.error <- with(
          block_estimates,
          sqrt(sum(std.error^2 * (N / N_overall)^2))
        )
      }


      ## we don't know if this is correct!
      ## matches lm_lin, two estimates per block
      if (is.null(data$cluster)) {
        design <- "Blocked"
        df <- nrow(data) - 2 * n_blocks
      } else {
        design <- "Block-clustered"
        # Also matches lm_lin for even sized clusters, should be conservative
        df <- sum(clust_per_block) - 2 * n_blocks
      }
    }

    return_frame <- data.frame(
      coefficients = diff,
      std.error = std.error,
      df = df,
      N = N_overall,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(data$weights)) {
    design <- paste0(design, " (weighted)")
  }

  return_list <- add_cis_pvals(return_frame, alpha, ci)

  # print(c("Pair Matched? ", pair_matched))

  return_list <- dim_like_return(
    return_list,
    alpha = alpha,
    formula = formula,
    conditions = list(condition1, condition2)
  )

  return_list[["design"]] <- design

  attr(return_list, "class") <- "difference_in_means"

  return(return_list)
}


difference_in_means_internal <- function(condition1 = NULL,
                                         condition2 = NULL,
                                         data,
                                         pair_matched = FALSE,
                                         alpha = .05,
                                         se_type = "default") {

  # Check that treatment status is uniform within cluster, checked here
  # so that the treatment vector t doesn't have to be built anywhere else
  if (!is.null(data$cluster)) {
    if (is.factor(data$cluster)) {
      data$cluster <- droplevels(data$cluster)
    }

    if (any(!tapply(data$t, data$cluster, function(x) all(x == x[1])))) {
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

  # Check to make sure multiple in each group if pair matched is false
  if (!pair_matched & (N2 == 1 | N1 == 1)) {
    stop(
      "If design is not pair-matched, every block must have at least two ",
      "treated and control units."
    )
  }

  df <- NA

  if (!is.null(data$cluster) && !pair_matched) {

    # For now, all clustered cases go to lm_robust
    # CR2 nests Gerber and Green 2012, p. 83, eq. 3.23 when clusters are
    # equal sizes (we think) and is more appropriate when clusters are of
    # different sizes

    X <- cbind(1, t = as.numeric(data$t == condition2))

    # TODO currently lm_robust_fit does too much, need to refactor it
    # if it will be used here in the long run
    cr2_out <- lm_robust_fit(
      y = data$y,
      X = cbind(1, t = as.numeric(data$t == condition2)),
      cluster = data$cluster,
      se_type = if(se_type == "none") "none" else "CR2",
      weights = data$weights,
      ci = TRUE,
      try_cholesky = TRUE,
      alpha = alpha,
      return_vcov = FALSE,
      has_int = TRUE,
      iv_stage = list(0)
    )

    diff <- coef(cr2_out)[2]
    std.error <- cr2_out$std.error[2]
    df <- cr2_out$df[2]
  } else {
    if (is.null(data$weights)) {
      diff <- mean(Y2) - mean(Y1)

      if (pair_matched || se_type == "none") {
        # Pair matched designs
        std.error <- NA
      } else {
        # Non-pair matched designs, unit level randomization
        var_Y2 <- var(Y2)
        var_Y1 <- var(Y1)

        std.error <- sqrt(var_Y2 / N2 + var_Y1 / N1)

        df <- std.error^4 /
          (
            (var_Y2 / N2)^2 / (N2 - 1) +
              (var_Y1 / N1)^2 / (N1 - 1)
          )
      }
    } else {
      if (pair_matched) {
        stop(
          "Cannot use `weights` with matched pairs design at the moment"
        )
      }

      X <- cbind(1, t = as.numeric(data$t == condition2))

      # print("Using lm_robust")
      # TODO currently lm_robust_fit does too much, need to refactor it
      # if it will be used here in the long run
      w_hc2_out <- lm_robust_fit(
        y = data$y,
        X = cbind(1, t = as.numeric(data$t == condition2)),
        se_type = ifelse(se_type == "none", "none", "HC2"),
        weights = data$weights,
        cluster = NULL,
        ci = TRUE,
        try_cholesky = TRUE,
        alpha = alpha,
        return_vcov = FALSE,
        has_int = TRUE,
        iv_stage = list(0)
      )

      diff <- coef(w_hc2_out)[2]
      std.error <- w_hc2_out$std.error[2]
      df <- w_hc2_out$df[2]
    }
  }

  return_frame <-
    data.frame(
      coefficients = diff,
      std.error = std.error,
      df = df,
      stringsAsFactors = FALSE
    )

  if (is.numeric(data$weights)) {
    return_frame$N <- sum(data$weights)
  } else {
    return_frame$N <- N
  }

  return(return_frame)
}
