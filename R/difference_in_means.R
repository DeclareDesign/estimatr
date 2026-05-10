#' Design-Based Difference-in-Means Estimator
#'
#' @param formula an object of class formula with one variable on the RHS
#' @param data A `data.frame`
#' @param blocks An optional bare (unquoted) name of the block variable
#' @param clusters An optional bare (unquoted) name of the cluster variable
#' @param weights the bare (unquoted) name of the weights variable
#' @param subset An optional bare (unquoted) expression specifying a subset
#' @param se_type `"default"` or `"none"`
#' @param condition1 value in treatment for the control condition
#' @param condition2 value in treatment for the treatment condition
#' @param ci logical. Whether to compute p-values and confidence intervals.
#' @param alpha The significance level, 0.05 by default.
#'
#' @details Selects the appropriate point estimate, standard errors, and degrees
#'   of freedom for unit randomized, cluster randomized, block randomized,
#'   block-cluster randomized, matched-pairs, and matched-pair cluster randomized
#'   designs.
#'
#'   If weights are specified, estimation is handed to [lm_robust()] with HC2
#'   standard errors.
#'
#' @return An object of class `"difference_in_means"`.
#'
#' @references Gerber, Alan S. and Donald P. Green. 2012. \emph{Field Experiments:
#'   Design, Analysis, and Interpretation}. New York: W.W. Norton.
#'
#'   Imai, Kosuke, Gary King, and Clayton Nall. 2009. "The Essential Role of
#'   Pair Matching in Cluster-Randomized Experiments." \emph{Statistical Science}
#'   24(1): 29-53. \doi{10.1214/08-STS274}.
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
  if (length(all.vars(rlang::f_rhs(rlang::eval_tidy(formula)))) > 1) {
    stop(
      "'formula' must have only one variable on the right-hand side: the ",
      "treatment variable."
    )
  }

  se_type <- match.arg(se_type)

  datargs <- rlang::enquos(
    formula = formula,
    weights = weights,
    subset = subset,
    block = blocks,
    cluster = clusters
  )
  data <- rlang::enquo(data)
  model_data <- clean_model_data(data = data, datargs, estimator = "dim")

  data <- data.frame(
    y = model_data$outcome,
    t = model_data$original_treatment,
    stringsAsFactors = FALSE
  )
  data$cluster <- model_data$cluster
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

  data <- subset.data.frame(data, t %in% c(condition1, condition2))

  nblocks <- NULL
  nclusters <- NULL

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
      nclusters <- return_frame[["nclusters"]]
      design <- "Clustered"
    }

  } else {
    pair_matched <- FALSE

    data <- subset.data.frame(data, t %in% c(condition1, condition2))

    clust_per_block <- check_clusters_blocks(data)

    if (any(clust_per_block == 1)) {
      stop("All `blocks` must have multiple units (or `clusters`)")
    } else if (all(clust_per_block == 2)) {
      pair_matched <- TRUE
    } else if (any(clust_per_block == 2) & any(clust_per_block > 2)) {
      pair_matched <- TRUE
      warning(
        "Some `blocks` have two units/`clusters` while other blocks ",
        "have more units/`clusters`. Using matched pairs variance estimator."
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

    N_overall <- with(block_estimates, sum(nobs))
    nclusters <- with(block_estimates, sum(nclusters))

    # Blocked design weighted combination (Gerber & Green 2012, p. 73, eq. 3.10)
    diff <- with(block_estimates, sum(coefficients * nobs / N_overall))

    df <- NA
    std.error <- NA
    nblocks <- nrow(block_estimates)

    if (pair_matched) {
      if (is.null(data$cluster)) {
        design <- "Matched-pair"

        # Pair matched, unit randomized (Gerber & Green 2012, p. 77, eq. 3.16)
        if (se_type != "none") {
          std.error <-
            with(
              block_estimates,
              sqrt(
                (1 / (nblocks * (nblocks - 1))) *
                  sum((coefficients - diff)^2)
              )
            )
        }
      } else {
        design <- "Matched-pair clustered"
        # Pair matched, cluster randomized (Imai, King & Nall 2009, p. 36, eq. 6)
        if (se_type != "none") {
          std.error <-
            with(
              block_estimates,
              sqrt(
                (nblocks / ((nblocks - 1) * N_overall^2)) *
                  sum((nobs * coefficients - (N_overall * diff) / nblocks)^2)
              )
            )
        }
      }

      # Imai et al. (2009, p. 37)
      df <- nblocks - 1
    } else {
      # Block randomized (Gerber & Green 2012, p. 74, fn. 17; Pashley & Miratrix 2021)
      if (se_type != "none") {
        std.error <- with(
          block_estimates,
          sqrt(sum(std.error^2 * (nobs / N_overall)^2))
        )
      }

      if (is.null(data$cluster)) {
        design <- "Blocked"
        # n - 2*nblocks: one DF used per block for control mean, one for treatment mean
        df <- nrow(data) - 2 * nblocks
      } else {
        design <- "Block-clustered"
        # Same logic applied to cluster-collapsed data
        df <- nclusters - 2 * nblocks
      }
    }

    return_frame <- data.frame(
      coefficients = diff,
      std.error = std.error,
      df = df,
      nobs = N_overall,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(data$weights)) {
    design <- paste0(design, " (weighted)")
  }

  return_list <- add_cis_pvals(return_frame, alpha, ci)

  return_list <- dim_like_return(
    return_list,
    alpha = alpha,
    formula = formula,
    conditions = list(condition1, condition2)
  )

  return_list[["design"]] <- design
  if (is.numeric(nblocks)) {
    return_list[["nblocks"]] <- nblocks
  }
  if (is.numeric(nclusters)) {
    return_list[["nclusters"]] <- nclusters
  }

  attr(return_list, "class") <- "difference_in_means"

  return(return_list)
}


difference_in_means_internal <- function(condition1 = NULL,
                                         condition2 = NULL,
                                         data,
                                         pair_matched = FALSE,
                                         alpha = .05,
                                         se_type = "default") {

  clustered <- !is.null(data$cluster)
  if (clustered) {
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

  if (!pair_matched & (N2 == 1 | N1 == 1)) {
    stop(
      "If design is not pair-matched, every block must have at least two ",
      "treated and control units."
    )
  }

  df <- NA
  nclusters <- NA

  if (clustered && !pair_matched) {

    cr2_out <- lm_robust_fit(
      y = data$y,
      X = cbind(1, t = as.numeric(data$t == condition2)),
      cluster = data$cluster,
      se_type = if (se_type == "none") "none" else "CR2",
      weights = data$weights,
      ci = TRUE,
      try_cholesky = TRUE,
      alpha = alpha,
      return_vcov = FALSE,
      has_int = TRUE,
      iv_stage = list(0)
    )

    diff <- coef(cr2_out)[2]
    std.error <- cr2_out[["std.error"]][2]
    df <- cr2_out[["df"]][2]
    nclusters <- cr2_out[["nclusters"]]
  } else {
    if (is.null(data$weights)) {
      diff <- mean(Y2) - mean(Y1)

      if (pair_matched || se_type == "none") {

        if (clustered) {
          nclusters <- 2
        }

        std.error <- NA
      } else {
        var_Y2 <- var(Y2)
        var_Y1 <- var(Y1)

        std.error <- sqrt(var_Y2 / N2 + var_Y1 / N1)

        # Welch-Satterthwaite degrees of freedom
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
    return_frame$nobs <- sum(data$weights)
  } else {
    return_frame$nobs <- N
  }

  if (is.numeric(nclusters)) {
    return_frame$nclusters <- nclusters
  }

  return(return_frame)
}
