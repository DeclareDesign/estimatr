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
#'   **Blocks of different sizes.** For unit randomized blocks, blocks are
#'   classified by how many units each arm holds rather than by how large the
#'   block is. A block with at least two treated and two control units carries
#'   its own Neyman variance. A block with a single treated or single control
#'   unit has no estimable within-block variance, so the variation across such
#'   blocks stands in for it. A design containing both kinds combines the two
#'   parts by squared share of the sample, following Pashley and Miratrix
#'   (2021). The `design` element of the returned object reports which case
#'   applied: `"Blocked"`, `"Matched-pair"`, `"Small blocks"`, or
#'   `"Hybrid blocked"`.
#'
#'   Two designs are refused, because the variance genuinely cannot be
#'   estimated: exactly one block with a singleton arm, and a set of
#'   different-sized such blocks in which one holds half or more of their
#'   units. Both messages suggest merging blocks or using [lm_robust()] with
#'   block fixed effects.
#'
#'   If weights are specified, estimation is handed to [lm_robust()] with HC2
#'   standard errors.
#'
#'   **Blocks of clusters.** Pashley and Miratrix treat treatment assigned
#'   within blocks, not blocks of clusters, so blocked designs with `clusters`
#'   use the earlier estimators. Every block must have at least two treated and
#'   two control clusters, unless the design is matched-pair clustered, where
#'   the variance is estimated across blocks. A block with a single treated or
#'   control cluster is refused: its within-block variance is not estimable,
#'   and estimating it anyway understates the standard error by roughly the
#'   block's cluster count.
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
#'   Pashley, Nicole E. and Luke W. Miratrix. 2021. "Insights on Variance
#'   Estimation for Blocked and Matched Pairs Designs." \emph{Journal of
#'   Educational and Behavioral Statistics} 46(3): 271-296.
#'   \doi{10.3102/1076998620946272}.
#'
#' @examples
#' set.seed(30)
#' dat <- data.frame(y = rnorm(100), z = rep(0:1, 50))
#'
#' # Unblocked, unclustered: the Welch-corrected two-sample difference
#' fit <- difference_in_means(y ~ z, data = dat)
#' fit
#' fit$design
#'
#' # Blocked designs use the Neyman variance within each block
#' dat_bl <- data.frame(
#'   bl = rep(1:10, each = 10),
#'   z  = rep(rep(0:1, each = 5), times = 10)
#' )
#' dat_bl$y <- rnorm(100) + 0.3 * dat_bl$z
#' difference_in_means(y ~ z, data = dat_bl, blocks = bl)
#'
#' # Matched pairs are recognised as such
#' dat_pr <- data.frame(pr = rep(1:50, each = 2), z = rep(c(0, 1), 50))
#' dat_pr$y <- rnorm(100) + 0.3 * dat_pr$z
#' difference_in_means(y ~ z, data = dat_pr, blocks = pr)$design
#'
#' # Blocks of unequal shape, which earlier versions refused, use the
#' # Pashley and Miratrix (2021) estimators. `design` reports which case
#' # applied rather than leaving it to be inferred from the block sizes.
#' dat_hy <- rbind(dat_bl[c("bl", "z", "y")],
#'                 transform(dat_pr[c("pr", "z", "y")], bl = pr + 100)[c("bl", "z", "y")])
#' difference_in_means(y ~ z, data = dat_hy, blocks = bl)$design
#'
#' # Clustered assignment
#' dat_cl <- data.frame(cl = rep(1:20, each = 5))
#' dat_cl$z <- rep(rep(0:1, each = 5), times = 10)
#' dat_cl$y <- rnorm(100) + 0.3 * dat_cl$z
#' difference_in_means(y ~ z, data = dat_cl, clusters = cl)
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

  } else if (is.null(data$cluster) && is.null(data$weights)) {

    # Unit-randomized blocks. Blocks are classified by how many units each arm
    # has, not by how large the block is: a block with one treated unit has no
    # estimable within-block variance whatever its size. See
    # blocked_variance_pm() for the estimators.
    blocked <- blocked_estimates(data, condition1, condition2)
    nblocks <- length(blocked$tau_k)
    N_overall <- sum(blocked$n_k)
    diff <- sum(blocked$tau_k * blocked$n_k) / N_overall

    if (se_type == "none") {
      std.error <- NA_real_
      df <- NA_real_
    } else {
      vc <- blocked_variance_pm(blocked$tau_k, blocked$n_k, blocked$var_k,
                                blocked$small, diff)
      std.error <- sqrt(vc$variance)
      df <- vc$df
    }

    design <- if (!any(blocked$small)) {
      "Blocked"
    } else if (all(blocked$small)) {
      if (all(blocked$n_k == 2)) "Matched-pair" else "Small blocks"
    } else {
      "Hybrid blocked"
    }

    return_frame <- data.frame(
      coefficients = diff,
      std.error = std.error,
      df = df,
      nobs = N_overall,
      stringsAsFactors = FALSE
    )

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

    # A block with one treated or one control cluster has no estimable
    # within-block variance: the singleton arm contributes nothing to the
    # between-cluster variability, so the estimate captures only the other arm
    # and comes out roughly a factor of the block's cluster count too small.
    # Exact enumeration: 6 clusters with 1 treated returns 0.17 of the true
    # variance. The matched-pair path is exempt because it estimates the
    # variance across blocks rather than within them.
    if (!pair_matched && !is.null(data$cluster)) {
      singleton_arm <- vapply(split(data, data$block), function(x) {
        min(length(unique(x$cluster[x$t == condition2])),
            length(unique(x$cluster[x$t == condition1]))) == 1L
      }, logical(1))

      if (any(singleton_arm)) {
        stop(
          "Some `blocks` have only one treated or one control `cluster`: ",
          paste(names(singleton_arm)[singleton_arm], collapse = ", "),
          ".\nThe within-block variance cannot be estimated from a single ",
          "cluster, and doing it anyway understates the standard error ",
          "severely. Merge those blocks with others so each has at least two ",
          "clusters per arm, or use `lm_robust()` with block fixed effects ",
          "and `clusters`."
        )
      }
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


# Per-block difference in means for a unit-randomized blocked design, with the
# within-block Neyman variance where both arms have at least two units. A block
# with a singleton arm is flagged `small`: its variance is not estimable from
# within the block, and blocked_variance_pm() estimates that part across blocks.
blocked_estimates <- function(data, condition1, condition2) {
  blocks <- split(data, data$block)

  out <- lapply(blocks, function(x) {
    Y2 <- x$y[x$t == condition2]
    Y1 <- x$y[x$t == condition1]
    n2 <- length(Y2)
    n1 <- length(Y1)

    if (n2 + n1 == 1L) {
      stop("All `blocks` must have multiple units (or `clusters`)")
    }
    if (n2 == 0L || n1 == 0L) {
      stop("Must have units with both treatment conditions within each block.")
    }

    small <- n2 == 1L || n1 == 1L
    c(
      tau = mean(Y2) - mean(Y1),
      n_k = n2 + n1,
      small = small,
      var_k = if (small) NA_real_ else var(Y2) / n2 + var(Y1) / n1
    )
  })

  out <- do.call(rbind, out)
  list(
    tau_k = out[, "tau"],
    n_k = out[, "n_k"],
    var_k = out[, "var_k"],
    small = as.logical(out[, "small"])
  )
}

# Variance of the blocked difference-in-means estimator, following Pashley and
# Miratrix (2021). Blocks where both arms have at least two units contribute
# their own Neyman variance (their equation 4). Blocks with a singleton arm
# have no within-block variance estimate, so the variation across those blocks
# stands in for it (their equations 5 and 8). A design with both kinds is the
# "hybrid" of their section 3.3, and the two parts combine by squared share of
# the sample.
#
# `tau_bk` is the overall estimate, the sample-weighted mean of the per-block
# estimates.
blocked_variance_pm <- function(tau_k, n_k, var_k, small, tau_bk) {
  n <- sum(n_k)

  # ---- big blocks: equation 4, weighted to their own share of the sample ----
  n_big <- sum(n_k[!small])
  k_big <- sum(!small)
  if (k_big > 0L) {
    var_big <- sum(n_k[!small]^2 * var_k[!small]) / n_big^2
    df_big <- n_big - 2 * k_big
  } else {
    var_big <- 0
    df_big <- 0
  }

  # ---- small blocks: equations 5 and 8 ----
  n_sb <- sum(n_k[small])
  k_small <- sum(small)
  if (k_small > 0L) {
    small_names <- names(tau_k)[small]
    if (k_small == 1L) {
      stop(
        "Only one block has a single treated or control unit (block ",
        small_names, "). The variance contributed by such blocks is ",
        "estimated from the variation across them, so at least two are ",
        "needed. Merge that block with another, drop it, or use ",
        "`lm_robust()` with block fixed effects."
      )
    }
    n_small <- n_k[small]
    tau_small <- sum(tau_k[small] * n_small) / n_sb
    dev2 <- (tau_k[small] - tau_small)^2

    if (length(unique(n_small)) == 1L) {
      # Equal sizes: equation 5, the usual matched-pairs estimator. Kept
      # separate because equation 8 is undefined at two equal-sized blocks.
      var_small <- sum(dev2) / (k_small * (k_small - 1))
    } else {
      # Equation 8, the unified estimator: no two blocks need share a size.
      # It requires every block to be under half the small-block sample, which
      # is what keeps the weights positive and the estimator conservative.
      if (any(n_small >= n_sb / 2)) {
        stop(
          "Blocks with a singleton treated or control unit cannot be handled ",
          "here: block ",
          paste(small_names[n_small >= n_sb / 2], collapse = ", "),
          " holds half or more of the units in such blocks, which makes the ",
          "variance estimator undefined. Merge blocks so that none dominates, ",
          "or use `lm_robust()` with block fixed effects."
        )
      }
      w <- n_small^2 / (n_sb - 2 * n_small)
      var_small <- sum(w * dev2) / (n_sb + sum(w))
    }
    df_small <- k_small - 1
  } else {
    var_small <- 0
    df_small <- 0
  }

  # ---- combine (section 3.3) ----
  v_big <- (n_big / n)^2 * var_big
  v_small <- (n_sb / n)^2 * var_small
  variance <- v_big + v_small

  # The paper stops at the variance. Welch-Satterthwaite across the two
  # components gives a df that reduces to n - 2K for an all-big design and to
  # K - 1 for an all-small one, matching what each literature uses on its own.
  denom <- (if (df_big > 0) v_big^2 / df_big else 0) +
    (if (df_small > 0) v_small^2 / df_small else 0)
  df <- if (denom > 0) variance^2 / denom else NA_real_

  list(variance = variance, df = df)
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
