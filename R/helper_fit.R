#' Internal method that creates linear fits
#'
#' @param y numeric outcome vector or matrix
#' @param X numeric design matrix
#' @param weights numeric weights vector
#' @param cluster numeric cluster vector
#' @param ci boolean, whether to return confidence intervals and p-values
#' @param se_type character denoting which kind of SEs to return
#' @param has_int logical, whether the model has an intercept
#' @param alpha numeric, test size for confidence intervals
#' @param return_vcov logical, whether to return the vcov matrix
#' @param return_fit logical, whether to return fitted values
#' @param try_cholesky logical. Solve by Cholesky decomposition of `X'X`
#'   rather than by the pivoted QR, falling back to the QR where the design is
#'   rank deficient. See [lm_robust()] for when the fast path is safe.
#' @param iv_stage list of length one or two for 2SLS stages
#' @param fe_rank integer, degrees of freedom absorbed by fixed effects
#' @param femat optional numeric matrix of fixed-effect dummies for the
#'   estimation sample. Only `CR2` needs it: its adjustment is built from
#'   cluster-level blocks of the hat matrix, which do not decompose the way the
#'   diagonal does. HC2 and HC3 take `fe_leverage` instead. `NULL` unless the
#'   requested `se_type` requires it.
#' @param fe_leverage numeric vector of per-observation leverage contributed by
#'   the absorbed fixed effects, or `NULL`. `h_ii` of the full design splits
#'   exactly into the demeaned-X leverage plus this term, for any number of FE
#'   factors, which is what makes HC2 and HC3 available under `fixed_effects`
#'   without building the dummy matrix.
#' @param linear_hypothesis optional hypotheses, in the form [lh_robust()]
#'   takes, whose `CR2` Satterthwaite degrees of freedom are returned as
#'   `hypothesis_df`. Ignored for every other `se_type`, where a combination of
#'   coefficients has the same degrees of freedom as each of them.
#' @param fe_rank_k integer, the absorbed rank counted in K by the clustered
#'   `se_type = "stata_nested"` factor: `fe_rank` less the fixed effects nested
#'   in the cluster variable. Read under that `se_type` alone; `fe_rank` still
#'   sets the residual degrees of freedom.
#' @param df_type character or `NULL`, the degrees of freedom; see
#'   [lm_robust()]. `"satterthwaite"` with `se_type = "HC2"` needs `femat`
#'   under absorbed fixed effects.
#'
#' @examples
#' # The fitter behind lm_robust(), exported for packages that have already
#' # built their own design matrix. Most users want lm_robust().
#' set.seed(45)
#' X <- cbind(`(Intercept)` = 1, x = rnorm(50))
#' y <- X[, "x"] + rnorm(50)
#'
#' lm_robust_fit(
#'   y = y, X = X,
#'   weights = NULL, cluster = NULL,
#'   se_type = "HC2", has_int = TRUE
#' )
#'
#' @export
lm_robust_fit <- function(y,
                          X,
                          weights,
                          cluster,
                          ci = TRUE,
                          se_type,
                          has_int,
                          alpha = 0.05,
                          return_vcov = TRUE,
                          return_fit = TRUE,
                          try_cholesky = FALSE,
                          iv_stage = list(0),
                          fe_rank = 0L,
                          fe_leverage = NULL,
                          femat = NULL,
                          linear_hypothesis = NULL,
                          fe_rank_k = fe_rank,
                          df_type = NULL) {

  # ----------
  # Check se type
  # ----------

  clustered <- !is.null(cluster)
  weighted <- !is.null(weights)
  se_type <- check_se_type(se_type, clustered, has_fe = (fe_rank > 0L),
                           has_fe_leverage = !is.null(fe_leverage),
                           fe_dummies = !is.null(femat), weighted = weighted)
  df_type <- check_df_type(df_type, se_type, has_fe = (fe_rank > 0L),
                           weighted = weighted)

  # HC2 with Bell-McCaffrey degrees of freedom: the dummies are for the
  # degrees-of-freedom pass alone, and the HC2 fit itself runs exactly as it
  # does by default, on the leverage decomposition.
  bm_df <- ci && identical(df_type, "satterthwaite") && se_type == "HC2"
  femat_df <- femat
  if (bm_df) femat <- NULL
  # Welch's degrees of freedom are read off the data as passed: prep_data()
  # rescales X and y by the weights.
  welch_df <- ci && identical(df_type, "welch")
  if (welch_df) {
    if (iv_stage[[1]] != 0) {
      stop("`df_type = \"welch\"` is not available for `iv_robust()`.")
    }
    welch_data <- list(X = X, y = as.matrix(y), weights = weights)
  }

  # -----------
  # Prep data for fitting
  # -----------

  data <- list(
    y = as.matrix(y),
    X = X
  )
  if (!is.null(femat)) {
    data[["femat"]] <- femat
    if (!is.double(data[["femat"]])) storage.mode(data[["femat"]]) <- "double"
  }

  # lm_solver maps these straight into Eigen, which requires doubles, and
  # tidy() needs a name for every outcome column. Callers passing matrices
  # directly (this function is exported) hit both (estimatr #269).
  if (!is.double(data[["y"]])) storage.mode(data[["y"]]) <- "double"
  if (!is.double(data[["X"]])) storage.mode(data[["X"]]) <- "double"

  # The model frame's row names are a character vector as long as the data, and
  # every object derived from the design matrix inherits a copy: the fitted
  # values from `X %*% beta`, the residuals formed from those, and the row
  # subset the cluster sort takes. Only `fitted.values` keeps them in the
  # returned object. Carrying them here and attaching them once at the end is
  # the same output for a third of the call at n = 100,000. `prep_data()`
  # reorders rows but never drops any, and the unsort below restores the
  # original order, so names captured now still line up at the end.
  ny <- ncol(data[["y"]])
  ynames <- colnames(data[["y"]])
  if (is.null(ynames)) {
    ynames <- if (ny == 1L) "y" else paste0("y", seq_len(ny))
  }
  multivariate <- ny > 1
  if (weighted) {
    data[["weights"]] <- weights
  }
  if (iv_stage[[1]] == 2) {
    data[["X_first_stage"]] <- iv_stage[[2]]
  }
  if (clustered) {
    # The C++ takes integer codes. A character cluster used to reach it and
    # come back as an Rcpp type error, or abort the process outright under a
    # debug build. Coercing here also makes `nclusters` right for a character
    # vector, which `length(unique())` would have got by luck rather than by
    # construction.
    if (!is.integer(cluster)) {
      cluster <- as.integer(as.factor(cluster))
    }
    data[["cluster"]] <- cluster
  }

  k <- ncol(data[["X"]])

  if (is.null(colnames(data[["X"]]))) {
    colnames(data[["X"]]) <- paste0("X", 1:k)
  }
  variable_names <- colnames(data[["X"]])

  data <- prep_data(
    data = data,
    se_type = se_type,
    clustered = clustered,
    weighted = weighted,
    iv_stage = iv_stage
  )

  # -----------
  # Estimate coefficients
  # -----------

  fit <-
    lm_solver(
      X = data[["X"]],
      y = data[["y"]],
      try_cholesky = try_cholesky
    )

  fit$beta_hat <- as.matrix(fit$beta_hat)
  dimnames(fit$beta_hat) <- list(variable_names, ynames)

  est_exists <- !is.na(fit$beta_hat)
  covs_used <- which(est_exists[, 1])
  # A zero-weight row contributes nothing to the fit, so it is not an
  # observation for the residual degrees of freedom or the HC1/stata scale
  # factors. `lm()` counts it the same way: ten zero weights among 100 rows
  # gave df.residual 98 here against `lm()`'s 88, which moved every classical
  # and HC1 standard error and every p-value. The row stays in `residuals` and
  # `fitted.values`, as `lm()` keeps it.
  N <- if (weighted) sum(weights > 0) else nrow(data[["X"]])

  x_rank <- length(covs_used)
  tot_rank <- x_rank + fe_rank

  if (multivariate) {
    return_list <- list(
      coefficients = fit$beta_hat,
      std.error = matrix(NA, k, ny),
      df = matrix(NA, k, ny)
    )
  } else {
    return_list <- list(
      coefficients = setNames(as.vector(fit$beta_hat), variable_names),
      std.error = NA,
      df = NA
    )
  }

  # ----------
  # Estimate variance
  # ----------

  # Set by the degenerate-variance rule below and read again by the F
  # statistic, which is computed after this block closes. Length 0 wherever
  # the rule does not run, which is what both readers test.
  no_variation <- logical(0)

  if (se_type != "none" || return_fit) {

    if (x_rank < ncol(data[["X"]])) {
      data <- drop_collinear(data, covs_used, weighted, iv_stage)
      fit$beta_hat <- fit$beta_hat[covs_used, , drop = FALSE]
    }

    fit_vals <- list()
    if (iv_stage[[1]] == 2) {
      X_name <- "X_first_stage"
      X_name_unweighted <- "X_first_stage_unweighted"
    } else {
      X_name <- "X"
      X_name_unweighted <- "Xunweighted"
    }

    # Subset only when a column was actually dropped for collinearity. The
    # subset is a full copy of the design matrix, and in the ordinary
    # full-rank case it copies it to itself.
    X_fit <- data[[X_name]]
    if (x_rank < ncol(X_fit)) {
      X_fit <- X_fit[, seq_len(x_rank), drop = FALSE]
    }
    fit_vals[["fitted.values"]] <- as.matrix(X_fit %*% fit$beta_hat)

    # `a - b` takes the first non-NULL dimnames of its operands. `data[["y"]]`
    # used to supply them, through the row names it carried, which left the
    # residuals' colnames those of the outcome: NULL for one outcome, the
    # outcome names for several. With the row names gone the fitted values win
    # that race instead and every fit would pick up `ynames`, which reaches
    # `res_var`, `r.squared` and `adj.r.squared` through `colSums()`. Set them
    # from the outcome, as before, rather than inherit them.
    ei <- data[["y"]] - fit_vals[["fitted.values"]]
    colnames(ei) <- colnames(data[["y"]])
    fit_vals[["ei"]] <- ei

    if (weighted) {
      fit_vals[["fitted.values.unweighted"]] <- as.matrix(
        data[[X_name_unweighted]] %*% fit$beta_hat
      )

      ei_unw <- data[["yunweighted"]] - fit_vals[["fitted.values.unweighted"]]
      colnames(ei_unw) <- colnames(data[["yunweighted"]])
      fit_vals[["ei.unweighted"]] <- ei_unw

      if (se_type %in% c("CR2", "CR3")) {
        data[["X"]] <- data[["weights"]] * data[["X"]]
        if (!is.null(data[["femat"]])) {
          data[["femat"]] <- data[["weights"]] * data[["femat"]]
        }
      }
    }

    if (iv_stage[[1]] == 2) {
      fit_vals[["fitted.values.iv"]] <- as.matrix(data[["X"]] %*% fit$beta_hat)
      ei_iv <- data[["y"]] - fit_vals[["fitted.values.iv"]]
      colnames(ei_iv) <- colnames(data[["y"]])
      fit_vals[["ei.iv"]] <- ei_iv
      if (weighted) {
        fit_vals[["ei.iv"]] <- data[["weights"]] * fit_vals[["ei.iv"]]
      }
      return_list[["ei.iv"]] <- fit_vals[["ei.iv"]]
    }

    if (se_type != "none") {

      # One row per hypothesis over the kept coefficients. car appends the
      # right-hand side as a last column, which has no place in the degrees of
      # freedom, and returns a single hypothesis as a vector, which rbind()
      # makes a one-row matrix.
      hypotheses <- NULL
      if (!is.null(linear_hypothesis) && se_type %in% c("CR2", "CR3") && ci && !multivariate) {
        hypotheses <- if (is.character(linear_hypothesis)) {
          rbind(car::makeHypothesis(variable_names, linear_hypothesis))[, seq_len(k), drop = FALSE]
        } else {
          rbind(linear_hypothesis)[, seq_len(k), drop = FALSE]
        }
        hypotheses <- hypotheses[, covs_used, drop = FALSE]
      }

      # A zero-weight row is not in the fit and contributes nothing to any
      # sandwich: its scaled X row is 0 and its leverage is 0. It is dropped
      # here rather than left for the C++ to ignore because CR2's per-cluster
      # adjustment reads every row of a cluster's block, and a zero-weight row
      # left in the block moved the weighted CR2 standard errors by 1.4e-3
      # relative, where clubSandwich, which drops such rows, did not move at
      # all. Dropping them also makes `J` the positive-weight cluster count
      # that the cluster degrees of freedom and the stata factor need.
      rows <- if (weighted) which(data[["weights"]] > 0) else seq_len(nrow(data[["X"]]))
      all_rows <- length(rows) == nrow(data[["X"]])
      take <- function(m) if (all_rows || is.null(m)) m else m[rows, , drop = FALSE]

      vcov_fit <- lm_variance(
        # Widening the design here is what makes HC2, HC3 and CR2 available
        # under `fixed_effects`. The C++ detects X.cols() > ncol(XtX_inv) and
        # rebuilds the meat from the full design, dropping the rank-deficient
        # dummy columns itself.
        X = take(if (!is.null(data[["femat"]]))
          cbind(data[["X"]], data[["femat"]])
          else data[["X"]]),
        Xunweighted = take(if (!is.null(data[["femat"]]) && weighted)
          cbind(data[["Xunweighted"]], data[["fematunweighted"]])
          else data[["Xunweighted"]]),
        XtX_inv = fit$XtX_inv,
        ei = take(if (se_type %in% c("CR2", "CR3") && weighted)
          fit_vals[["ei.unweighted"]]
          else fit_vals[["ei"]]),
        weight_mean = data[["weight_mean"]],
        cluster = if (all_rows) data[["cluster"]] else data[["cluster"]][rows],
        J = data[["J_eff"]],
        ci = ci,
        # The C++ knows one Stata factor; the nested convention differs from
        # it only in K, which is the rank passed below.
        se_type = if (se_type == "stata_nested") "stata" else se_type,
        fe_rank = if (se_type == "stata_nested") fe_rank_k else fe_rank,
        # Only HC2/HC3 consume this; it is NULL for every other se_type and for
        # multi-way FE, so the C++ falls back to the plain hat value.
        fe_leverage = if (se_type %in% c("HC2", "HC3") && !is.null(fe_leverage))
          fe_leverage[rows]
          else NULL,
        n_eff = N,
        hypotheses = hypotheses
      )

      # One rule for what HC2 and HC3 report at a full-leverage observation.
      # lm_variance() discards an observation from the variance at
      # h > 1 - sqrt(eps), returns the size of that set as
      # `n_leverage_near_one`, and returns for each coefficient the share of
      # its classical variance those rows carry (`var_not_estimable`). A
      # coefficient whose share exceeds sqrt(eps) is identified by the
      # discarded rows alone, so its standard error is NA; by
      # Frisch-Waugh-Lovell every other coefficient keeps the standard error of
      # the design without those rows, to machine precision. The fit warns if
      # and only if some standard error is NA for that reason. Where none is,
      # the count is a message rather than a warning: a singleton fixed-effect
      # group is the ordinary case, its demeaned row is zero, and every
      # reported standard error is the one the design without it gives.
      #
      # Separately, a variance diagonal can come back negative on a design
      # close enough to singular for the sandwich, a difference of large and
      # nearly equal quantities, to lose the difference to rounding. `sqrt()`
      # would then emit R's bare "NaNs produced", which names neither the fit
      # nor the reason. The NaN is kept and the reason is attached.
      var_hat <- diag(vcov_fit$Vcov_hat)
      neg_var <- !is.na(var_hat) & var_hat < 0
      var_hat[neg_var] <- NaN
      return_list$std.error[est_exists] <- sqrt(var_hat)

      not_estimable <- logical(length(var_hat))
      frac_dropped <- vcov_fit[["var_not_estimable"]]
      if (length(frac_dropped) == length(var_hat)) {
        not_estimable <- frac_dropped > sqrt(.Machine$double.eps)
        if (any(not_estimable)) {
          se_vec <- return_list$std.error[est_exists]
          se_vec[not_estimable] <- NA_real_
          return_list$std.error[est_exists] <- se_vec
        }
      }
      # A variance of numerically zero is rounding error, not precision, and
      # nothing above catches it. The leverage rule reaches an observation the
      # fit reproduces exactly; this reaches a COEFFICIENT whose variance is
      # assembled only from such observations, which needs no observation at
      # leverage 1 and so is invisible to it. Both shapes occur: an outcome the
      # regressors fit exactly returns standard errors of 6e-17, a t of 1e16
      # and p = 0 for everything, and a sparse factor level confined to one arm
      # returns 8e-18 for that arm's contrast alone, with three stars on an
      # estimate of 2e-17, while the rest of the fit is ordinary.
      #
      # The test is one ratio, against the classical variance the coefficient
      # would carry if the regressors explained nothing: tss/df times
      # XtX_inv(j,j), which is dimensionless and so comparable across columns
      # of unlike scale, as `var_not_estimable` above is. The ratio is about
      # 1 - R^2 for a classical fit, so a fit has to reproduce its outcome to
      # the last bit to be caught; the tolerance is `.Machine$double.eps`
      # rather than its square root for that reason. It applies to every
      # se_type, since a classical fit returns the same 6e-17 here.
      dfres <- N - tot_rank
      if (dfres > 0 && length(var_hat) == x_rank * ny) {
        y_fit <- data[["y"]]
        rss <- colSums(fit_vals[["ei"]]^2 * data[["weight_mean"]])
        # `ei` and `y` are both on the weighted scale, so the ratio they form
        # is scale free; the centring uses the mean of the weighted outcome
        # rather than the weighted mean, which moves a tolerance reference and
        # nothing that is reported.
        tss_scale <- if (has_int) {
          ym <- y_fit - rep(colMeans(y_fit), each = nrow(y_fit))
          colSums(ym * ym * data[["weight_mean"]])
        } else {
          colSums(y_fit * y_fit * data[["weight_mean"]])
        }
        reference <- rep(tss_scale / dfres, each = x_rank) *
          rep(diag(fit$XtX_inv), times = ny)
        no_variation <- !is.na(var_hat) &
          var_hat <= .Machine$double.eps * reference
        # An outcome that is constant has no variation to estimate anything
        # from and a reference of 0, which no ratio can flag. Its variances
        # come back as exactly 0 rather than as rounding error.
        flat <- !(tss_scale > 0)
        if (any(flat)) no_variation[rep(flat, each = x_rank)] <- TRUE

        if (any(no_variation)) {
          se_vec <- return_list$std.error[est_exists]
          se_vec[no_variation] <- NA_real_
          return_list$std.error[est_exists] <- se_vec

          # Two readings, and the difference is what the user should do next.
          # An outcome reproduced exactly is a statement about the fit; a
          # handful of coefficients is a statement about those columns.
          exact <- flat | rss <= .Machine$double.eps * tss_scale
          if (any(exact)) {
            message(
              "The model reproduces ",
              if (ny > 1) paste0(paste(ynames[exact], collapse = ", "), " ") else "the outcome ",
              "exactly, so there is no residual variation left to estimate a ",
              "variance from and every standard error for ",
              if (ny > 1 && sum(exact) == 1) "it" else if (ny > 1) "them" else "it",
              " is NA. What the variance calculation returns here is zero or ",
              "rounding error, not precision."
            )
          }
          named <- no_variation & !rep(exact, each = x_rank)
          if (any(named)) {
            # Two mechanisms reach a coefficient whose variance is numerically
            # zero while the fit around it is ordinary, and a sentence that
            # named only the first was false on a paired conjoint. Without
            # clusters, the rows that identify the coefficient are each fitted
            # exactly. With clusters, the variance is built from within-cluster
            # sums of scores, and those can cancel to zero while every residual
            # in the cluster is far from it. The clustered wording also covers
            # a clustered fit whose identifying rows are exact, since residuals
            # that are all zero sum to zero.
            one <- sum(named) == 1
            message(
              "The standard error is NA for ",
              paste(unique(rep(variable_names[covs_used], ny)[named]),
                    collapse = ", "),
              ": ",
              if (clustered) {
                paste0(
                  "the cluster-level scores for ",
                  if (one) "that coefficient sum " else "those coefficients sum ",
                  "to numerically zero, so the cluster-robust variance is ",
                  "assembled from nothing, although the residuals themselves ",
                  "are not zero. That happens when a regressor is constant ",
                  "within cluster and the residuals cancel over the cluster."
                )
              } else {
                paste0(
                  "every observation that identifies ",
                  if (one) "that coefficient is " else "those coefficients is ",
                  "fitted exactly, so the variance is built from residuals ",
                  "that are all numerically zero."
                )
              },
              " Reported as a number it would be zero or rounding error, and ",
              "its t statistic would put p at 0."
            )
          }
        }
      }

      n_lev <- vcov_fit[["n_leverage_near_one"]]
      if (!isTRUE(n_lev > 0)) n_lev <- 0L
      return_list[["n_leverage_near_one"]] <- as.integer(n_lev)

      if (any(neg_var)) {
        warning(
          sum(neg_var), " of ", length(var_hat), " variance estimates came out ",
          "negative, so those standard errors are NaN. The design is close to ",
          "singular, and the sandwich estimator loses the difference between ",
          "two nearly equal quantities to rounding. Drop covariates, or use ",
          "`se_type = \"classical\"`."
        )
      }

      if (n_lev > 0) {
        lev_count <- paste0(
          n_lev, if (n_lev == 1) " observation has " else " observations have ",
          "a computed leverage at or near 1, which happens when the design is ",
          "close to saturated and the observation is fitted exactly or nearly ",
          "so. `se_type = \"", se_type, "\"` divides by (1 - leverage), so ",
          if (n_lev == 1) "it is " else "they are ",
          "dropped from the variance."
        )
        if (any(not_estimable)) {
          warning(
            lev_count, " The standard error is NA for ",
            # The design is shared across outcomes, so a multivariate fit
            # carries the same non-estimable coefficient once per outcome and
            # would otherwise name it that many times. One name per
            # coefficient is what the collinear-drop message already gives.
            paste(unique(rep(variable_names, ncol(est_exists))[est_exists][not_estimable]),
                  collapse = ", "),
            ", whose coefficient ",
            if (n_lev == 1) "that observation alone identifies" else "those observations alone identify",
            ": dropping ", if (n_lev == 1) "it " else "them ",
            "leaves nothing to estimate that variance from. Every other ",
            "standard error in this fit is the one the design without ",
            if (n_lev == 1) "it " else "them ",
            "gives. Use `se_type = \"HC1\"` or `\"classical\"`, or drop ",
            "covariates, to use every observation."
          )
        } else {
          message(
            lev_count, " No standard error is affected: each is the one the ",
            "design without ", if (n_lev == 1) "it " else "them ", "gives."
          )
        }
      }

      if (ci) {
        return_list$df[est_exists] <-
          ifelse(vcov_fit$dof == -99, NA, vcov_fit$dof)
      }
      if (!is.null(hypotheses)) {
        return_list[["hypothesis_df"]] <- vcov_fit$hypothesis_dof
      }

      # Bell and McCaffrey's degrees of freedom for HC2, from the design alone.
      # Only the degrees of freedom are replaced, so every coefficient and
      # variance above is the default HC2 fit's.
      if (bm_df || welch_df) {
        bm_hypotheses <- NULL
        if (!is.null(linear_hypothesis) && !multivariate) {
          bm_hypotheses <- if (is.character(linear_hypothesis)) {
            rbind(car::makeHypothesis(variable_names, linear_hypothesis))[, seq_len(k), drop = FALSE]
          } else {
            rbind(linear_hypothesis)[, seq_len(k), drop = FALSE]
          }
          bm_hypotheses <- bm_hypotheses[, covs_used, drop = FALSE]
        }
      }
      if (welch_df) {
        welch <- hc2_welch_df(
          X = welch_data[["X"]],
          y = welch_data[["y"]],
          weights = welch_data[["weights"]],
          hypotheses = bm_hypotheses
        )
        return_list$df[est_exists] <- welch[seq_len(k), ]
        if (!is.null(bm_hypotheses)) {
          return_list[["hypothesis_df"]] <- welch[-seq_len(k), 1]
        }
      }
      if (bm_df) {
        bm <- hc2_bm_df(
          X = take(data[["X"]]),
          femat = femat_df,
          hypotheses = bm_hypotheses
        )
        return_list$df[est_exists] <- rep(bm$df, ny)
        if (!is.null(bm_hypotheses)) {
          return_list[["hypothesis_df"]] <- bm$hypothesis_df
        }
      }
    }
  }

  # ----------
  # Augment return object
  # ----------

  return_list <- add_cis_pvals(return_list, alpha, ci && se_type != "none")

  if (return_fit) {

    # Weighted fits carry residuals on both scales; the unweighted pair is the
    # one on the scale of the data, so it is what gets returned.
    if (weighted) {
      fitted.vals_name <- "fitted.values.unweighted"
      resid_name <- "ei.unweighted"
    } else {
      fitted.vals_name <- "fitted.values"
      resid_name <- "ei"
    }

    return_list[["fitted.values"]] <- as.matrix(fit_vals[[fitted.vals_name]])
    # The residuals are returned unnamed, as 1.0.6 returned them, and until now
    # lm_return() dropped the names at the very end. They are the model frame's
    # row names, inherited from the design matrix through `X %*% beta`, and
    # letting them travel this far means they are copied through the cluster
    # unsort and again by `drop()` before being thrown away. Cutting them here
    # is the same output for 6 ms less at n = 100,000, a third of a plain fit.
    # Guarded: the assignment duplicates the matrix, and an absorbed
    # fixed-effects fit runs on a demeaned design that never had row names.
    resid_mat <- as.matrix(fit_vals[[resid_name]])
    if (!is.null(rownames(resid_mat))) rownames(resid_mat) <- NULL
    return_list[["residuals"]] <- resid_mat

    if (clustered && se_type != "none") {
      # prep_data sorted the rows by cluster; put them back
      unsort <- order(data[["cl_ord"]])
      return_list[["fitted.values"]] <- return_list[["fitted.values"]][unsort, , drop = FALSE]
      return_list[["residuals"]] <- return_list[["residuals"]][unsort, , drop = FALSE]
    }

    colnames(return_list[["fitted.values"]]) <- ynames
    colnames(return_list[["residuals"]]) <- ynames
  }

  return_list[["term"]] <- variable_names
  return_list[["outcome"]] <- ynames
  return_list[["alpha"]] <- alpha
  return_list[["se_type"]] <- se_type
  return_list[["weighted"]] <- weighted
  return_list[["fes"]] <- (fe_rank > 0L)
  return_list[["clustered"]] <- clustered
  return_list[["df.residual"]] <- N - tot_rank
  return_list[["nobs"]] <- N
  if (clustered) {
    return_list[["nclusters"]] <- data[["J_eff"]]
  }
  return_list[["k"]] <- k
  return_list[["rank"]] <- x_rank

  if (se_type != "none") {

    # `ei` is the weighted residual sqrt(w / mean(w)) * e, so this is the
    # weighted residual variance; `weight_mean` is 1 for an unweighted fit.
    # Unnamed, as it always was for one outcome: a multivariate `ei` carries
    # the outcome names, and `r.squared` is built from this.
    return_list[["res_var"]] <-
      unname(colSums(fit_vals[["ei"]]^2 * data[["weight_mean"]])) / return_list[["df.residual"]]

    tss_r2s <- get_r2s(
      y = data[["y"]],
      return_list = return_list,
      has_int = has_int,
      yunweighted = data[["yunweighted"]],
      weights = data[["weights"]],
      weight_mean = data[["weight_mean"]]
    )

    # Absorbing fixed effects removes the intercept column from the design, so
    # for the F statistic there is no intercept to skip: every coefficient is a
    # regressor under test. Treating one as an intercept anyway tested a
    # hypothesis one restriction short of the intended one, and for a
    # single-regressor FE model it took nomdf to zero and dropped the statistic
    # altogether. Only the F statistic is affected; get_r2s() still needs the
    # formula's own intercept flag.
    fstat_int <- if (fe_rank > 0L) 0L else has_int
    nomdf <- x_rank - fstat_int
    if (clustered) {
      dendf <- data[["J_eff"]] - 1
    } else {
      dendf <- return_list[["df.residual"]]
    }

    if (nomdf > 0) {
      f <- get_fstat(
        tss_r2s = tss_r2s,
        return_list = return_list,
        iv_ei = fit_vals[["ei.iv"]],
        nomdf = nomdf,
        dendf = dendf,
        vcov_fit = vcov_fit,
        has_int = fstat_int,
        iv_stage = iv_stage,
        no_variation = no_variation
      )
    } else {
      f <- NULL
    }

    return_list <- c(return_list, tss_r2s)
    return_list[["fstatistic"]] <- f

    if (return_vcov) {
      return_list[["vcov"]] <- vcov_fit$Vcov_hat
      if (multivariate) {
        coef_names <- lapply(
          seq_len(ncol(est_exists)),
          function(j) return_list$term[est_exists[, j]]
        )

        outcome_coef_names <- paste0(
          rep(paste0(return_list[["outcome"]], ":"),
              times = vapply(coef_names, length, integer(1))),
          unlist(coef_names, FALSE, FALSE)
        )

        dimnames(return_list[["vcov"]]) <- list(
          outcome_coef_names,
          outcome_coef_names
        )
      } else {
        dimnames(return_list[["vcov"]]) <- list(
          return_list$term[est_exists],
          return_list$term[est_exists]
        )
      }
    }
  }

  attr(return_list, "class") <- "lm_robust"

  return(return_list)
}

# The default warning below is `rlang::warn(.frequency = "once")`: absorbed
# fixed effects are usually fitted in a loop or a simulation, and a per-call
# warning would be noise rather than a notice. Once per session, keyed by
# `.frequency_id`. rlang always signals these under testthat, so the tests
# still see them on every call.
#' @importFrom rlang warn
#' @noRd
check_se_type <- function(se_type, clustered, has_fe = FALSE,
                          has_fe_leverage = FALSE,
                          fe_dummies = FALSE, weighted = FALSE) {

  stata_se_types <- c("stata", "stata_nested")
  cl_se_types  <- c("CR0", "CR2", "CR3", stata_se_types)
  rob_se_types <- c("HC0", "HC1", "HC2", "HC3", "classical", stata_se_types)

  if (clustered) {
    if (is.null(se_type)) {
      # CR2 under `fixed_effects` is computable, but only by expanding the
      # dummies, which is O(g^3) in the number of levels and would silently
      # undo the reason for absorbing them. The default stays on the fast
      # estimator and says so, rather than quietly returning a different
      # number from 1.0.6 (which defaulted to CR2 here).
      if (has_fe) {
        warn(
          paste0(
            "With `fixed_effects` and `clusters`, `se_type` defaults to ",
            "\"CR0\"; estimatr 1.0.6 defaulted to \"CR2\" here.\n",
            "CR2 is available by asking for it, `se_type = \"CR2\"`, at the ",
            "cost of expanding the fixed effects into dummies.\n",
            "Write `se_type = \"CR0\"` to accept this default and remove this ",
            "warning."
          ),
          .frequency = "once",
          .frequency_id = "estimatr_fe_clustered_default"
        )
      }
      se_type <- if (has_fe) "CR0" else "CR2"
    } else if (!(se_type %in% c(cl_se_types, "none"))) {
      stop(
        "`se_type` must be either 'CR0', 'stata', 'stata_nested', 'CR2', ",
        "'CR3', or 'none' when ",
        "`clusters` are specified.\nYou passed: ", se_type
      )
    }
  } else {
    if (is.null(se_type)) {
      # HC2 for every unclustered fit, with or without fixed effects and
      # whatever their number: fe_leverage() makes it exact and cheap, so
      # there is nothing to trade away and nothing to warn about. This is also
      # what 1.0.6 returned.
      se_type <- "HC2"
    } else if (se_type %in% setdiff(cl_se_types, stata_se_types)) {
      stop(
        "`se_type` must be either 'HC0', 'HC1', 'stata', 'HC2', 'HC3', ",
        "'classical' or 'none' with no `clusters`.\nYou passed: ", se_type,
        " which is reserved for a case with clusters."
      )
    } else if (!(se_type %in% c(rob_se_types, "none"))) {
      stop(
        "`se_type` must be either 'HC0', 'HC1', 'stata', 'HC2', 'HC3', ",
        "'classical' or 'none' with no `clusters`.\nYou passed: ", se_type
      )
    }
    # Without clusters nothing is nested in anything, and both Stata
    # conventions are HC1.
    if (se_type %in% stata_se_types) se_type <- "HC1"
  }

  # Only now the fixed-effects branches. They used to run first and return
  # before the menu above, so `fixed_effects` with `clusters` and
  # `se_type = "HC2"` was accepted and silently ignored the clusters, and
  # `fixed_effects` with `se_type = "CR2"` and no clusters reached the C++ with
  # a NULL cluster vector. 1.0.6 errored on both. Whether an se_type is on the
  # menu does not depend on the fixed effects, so the menu is checked first.

  # HC2 and HC3 need the hat values of the full [dummies | X] design. Those
  # split exactly into the demeaned-X leverage plus fe_leverage(), for ANY
  # number of fixed-effect factors, so they are exact and cheap here and carry
  # no restriction. The identity is stated and verified in fe_leverage().
  if (has_fe && has_fe_leverage && se_type %in% c("HC2", "HC3")) {
    return(se_type)
  }

  # CR2 is the one case left that has to materialise the dummies: it needs the
  # whole per-cluster block of the hat matrix, which does not decompose.
  if (has_fe && se_type %in% c("HC2", "HC3", "CR2", "CR3")) {
    # 1.0.6 refused this combination too: the CR2 adjustment needs the design
    # weighted by w rather than by sqrt(w), and the absorbed dummies are not
    # carried in a form that supports it.
    if (se_type %in% c("CR2", "CR3") && weighted) {
      stop(
        "`se_type = \"", se_type, "\"` cannot be combined with `weights` and ",
        "`fixed_effects`.\nUse `se_type = \"stata\"`, or replace ",
        "`fixed_effects` with explicit dummies:\n",
        "  lm_robust(y ~ x + factor(fe_var), weights = w, clusters = cl, ",
        "se_type = \"", se_type, "\")"
      )
    }
    if (!fe_dummies) {
      stop(
        "`se_type = \"", se_type, "\"` requires hat values from the full ",
        "[X | FE dummies] design matrix, which were not supplied.\n",
        "Call `lm_robust()` or `iv_robust()`, which build them, or pass ",
        "`femat` to `lm_robust_fit()`."
      )
    }
  }

  return(se_type)
}

# The Bell and McCaffrey (2002) degrees of freedom for HC2. A linear
# combination c'b = sum_i a_i y_i of the coefficients has HC2 variance
# V = sum_i a_i^2 e_i^2 / (1 - h_i), and its degrees of freedom are
# 2 E[V]^2 / Var[V] with the errors normal under the working covariance the
# weights imply, W^(-1), the identity when unweighted. In the sqrt(w)-scaled
# design, with Z an orthonormal basis of it, H = ZZ' and s_i = a_i^2 / (1 - h_i),
# V is a quadratic form in B = S^(1/2) (I - H) S^(1/2), so
#     df = tr(B)^2 / ||B||_F^2,
#     tr(B) = sum s_i (1 - h_i),
#     ||B||_F^2 = sum s_i^2 (1 - 2 h_i) + ||Z'SZ||_F^2,
# which is O(N p^2) per combination rather than the N by N matrix B. These are
# clubSandwich's CR2 Satterthwaite degrees of freedom with every observation
# its own cluster and `inverse_var = TRUE`, the convention HC2's variance
# follows under weights. An observation HC2 discards for leverage at or near 1
# carries s_i = 0, as it carries nothing in the variance.
#
# Under absorbed fixed effects X is demeaned and so orthogonal to the dummies
# D, and H = H_X + H_D with H_D = D P D', P the pseudo-inverse of D'D. An
# orthonormal basis of [X, D] is N by (number of levels) dense, which at 500 levels and 1e5 rows
# took 70 seconds. Every H_D term reduces instead to sums over the nonzeros of
# D, which number N times the factor count: h_i gains d_i' P d_i, and
#     ||Z'SZ||_F^2 = ||Z_X'SZ_X||_F^2 + 2 tr((D'SZ_X)' P (D'SZ_X)) + tr(KPKP),
# with K = D'SD. P and K are diagonal under one factor.
#' @noRd
hc2_bm_df <- function(X, femat = NULL, hypotheses = NULL) {
  r <- ncol(X)
  qx <- qr(X)
  stopifnot(qx$rank == r)
  Z <- qr.Q(qx)
  # Column j of Q R^(-T) is the a for coefficient j.
  A <- Z %*% t(backsolve(qr.R(qx), diag(r)))
  h <- rowSums(Z^2)

  if (!is.null(femat)) {
    G <- ncol(femat)
    nz <- which(femat != 0, arr.ind = TRUE)
    nz_row <- nz[, 1]
    nz_col <- nz[, 2]
    nz_val <- femat[nz]
    # Every pair of nonzeros that share a row, for D' diag(v) D.
    pairs <- merge(
      data.frame(row = nz_row, c1 = nz_col, v1 = nz_val),
      data.frame(row = nz_row, c2 = nz_col, v2 = nz_val),
      by = "row"
    )
    pair_ix <- pairs$c1 + (pairs$c2 - 1L) * G
    pair_val <- pairs$v1 * pairs$v2
    dvd <- function(v) {
      out <- matrix(0, G, G)
      agg <- rowsum(v[pairs$row] * pair_val, pair_ix)
      out[as.integer(rownames(agg))] <- agg
      out
    }

    DtD <- dvd(rep(1, nrow(X)))
    one_way <- all(DtD[row(DtD) != col(DtD)] == 0)
    if (one_way) {
      d <- diag(DtD)
      P <- diag(ifelse(d > 0, 1 / d, 0), G)
    } else {
      e <- eigen(DtD, symmetric = TRUE)
      pos <- e$values > max(e$values) * 1e-10
      P <- e$vectors[, pos, drop = FALSE] %*%
        (t(e$vectors[, pos, drop = FALSE]) / e$values[pos])
    }
    h <- h + as.vector(rowsum(P[pair_ix] * pair_val, pairs$row,
                              reorder = TRUE)[as.character(seq_len(nrow(X))), 1])
  }

  keep <- h <= 1 - sqrt(.Machine$double.eps)

  bm <- function(a) {
    s <- ifelse(keep, a^2 / (1 - h), 0)
    tr <- sum(s * (1 - h))
    fro <- sum(s^2 * (1 - 2 * h)) + sum(crossprod(Z, Z * s)^2)
    if (!is.null(femat)) {
      DSZ <- matrix(0, G, r)
      agg <- rowsum(s[nz_row] * nz_val * Z[nz_row, , drop = FALSE], nz_col)
      DSZ[as.integer(rownames(agg)), ] <- agg
      K <- dvd(s)
      fro <- fro + 2 * sum(DSZ * (P %*% DSZ)) +
        if (one_way) sum((diag(K) * diag(P))^2) else {
          KP <- K %*% P
          sum(KP * t(KP))
        }
    }
    if (fro > 0) tr^2 / fro else NA_real_
  }

  list(
    df = apply(A, 2, bm),
    hypothesis_df = if (!is.null(hypotheses))
      apply(A %*% t(hypotheses), 2, bm)
  )
}

# Welch-Satterthwaite degrees of freedom for HC2 on a two-group design: an
# intercept and one regressor taking two values on the positive-weight rows.
# Every coefficient, and every hypothesis, is then a combination of the two
# groups' weighted means, whose HC2 variances are
#     v_g = sum(w^2 e^2 / (1 - w / W_g)) / W_g^2,
# so its df is Welch's on its two terms. With constant weights v_g is
# var(y_g) / n_g, and the df of the slope is the textbook Welch df. Returns one
# row per coefficient, then one per hypothesis, and one column per outcome.
#' @noRd
hc2_welch_df <- function(X, y, weights = NULL, hypotheses = NULL) {
  if (is.null(weights)) weights <- rep(1, nrow(X))
  keep <- weights > 0
  X <- X[keep, , drop = FALSE]
  int_col <- apply(X, 2, function(x) all(x == 1))
  vals <- if (ncol(X) == 2 && sum(int_col) == 1) sort(unique(X[, !int_col])) else NULL
  if (length(vals) != 2) {
    stop(
      "`df_type = \"welch\"` needs a two-group design: an intercept and one ",
      "regressor taking two values, as in `y ~ z` with a binary `z`, and no ",
      "`fixed_effects`. Use `df_type = \"satterthwaite\"` for other designs."
    )
  }
  a <- vals[1]
  b <- vals[2]
  # Each row is the weight a coefficient puts on the two group means.
  A <- matrix(0, 2, 2)
  A[int_col, ] <- c(b, -a) / (b - a)
  A[!int_col, ] <- c(-1, 1) / (b - a)
  if (!is.null(hypotheses)) A <- rbind(A, hypotheses %*% A)

  w <- weights[keep]
  in_b <- X[, !int_col] == b
  n_g <- c(sum(!in_b), sum(in_b))
  if (any(n_g < 2)) return(matrix(NA_real_, nrow(A), ncol(y)))
  apply(y[keep, , drop = FALSE], 2, function(y_k) {
    v <- vapply(list(!in_b, in_b), function(in_g) {
      w_g <- w[in_g]
      w_sum <- sum(w_g)
      e <- y_k[in_g] - sum(w_g * y_k[in_g]) / w_sum
      sum(w_g^2 * e^2 / (1 - w_g / w_sum)) / w_sum^2
    }, numeric(1))
    terms <- sweep(A^2, 2, v, "*")
    df <- rowSums(terms)^2 / rowSums(sweep(terms^2, 2, n_g - 1, "/"))
    df[!is.finite(df)] <- NA_real_
    df
  })
}

# `df_type = NULL` keeps each se_type's own degrees of freedom: the residual
# N - K for the HC family and classical, Satterthwaite for CR2 and CR3, and
# J - 1 for CR0 and the Stata factors. Naming that default is accepted. The
# alternatives are for HC2: Satterthwaite, where Bell and McCaffrey derived
# it, and Welch, on a two-group design.
#' @noRd
check_df_type <- function(df_type, se_type, has_fe = FALSE, weighted = FALSE) {
  if (is.null(df_type) || se_type == "none") return(NULL)
  if (!(identical(df_type, "residual") || identical(df_type, "satterthwaite") ||
        identical(df_type, "welch"))) {
    stop("`df_type` must be NULL, 'residual', 'satterthwaite', or 'welch'.\nYou passed: ",
         paste(df_type, collapse = ", "))
  }
  if (df_type == "welch" && se_type == "HC2") {
    if (has_fe) {
      stop(
        "`df_type = \"welch\"` is not available with `fixed_effects`: it ",
        "needs a two-group design, an intercept and one binary regressor."
      )
    }
    return(df_type)
  }
  own <- if (se_type %in% c("CR2", "CR3")) {
    "satterthwaite"
  } else if (se_type %in% c("CR0", "stata", "stata_nested")) {
    "clusters"
  } else {
    "residual"
  }
  if (df_type == own) return(df_type)
  if (df_type == "satterthwaite" && se_type == "HC2") {
    if (has_fe && weighted) {
      stop(
        "`df_type = \"satterthwaite\"` with `se_type = \"HC2\"` cannot be ",
        "combined with `weights` and `fixed_effects`, for the reason `\"CR2\"` ",
        "cannot: it expands the fixed effects into dummies.\nReplace ",
        "`fixed_effects` with explicit dummies in the formula."
      )
    }
    return(df_type)
  }
  stop(
    "`df_type = \"", df_type, "\"` is not available with `se_type = \"",
    se_type, "\"`.\n`\"satterthwaite\"` is available with \"HC2\", \"CR2\", ",
    "and \"CR3\"; `\"welch\"` with \"HC2\"; `\"residual\"` with \"HC0\", ",
    "\"HC1\", \"HC2\", \"HC3\", and \"classical\"."
  )
}


get_r2s <- function(y, return_list, has_int, yunweighted, weights, weight_mean) {

  N <- nrow(y)
  if (return_list[["weighted"]]) {
    if (has_int) {
      # weights = sqrt(w/μ_w) from prep_data, so weights^2 = w/μ_w.
      # Standard weighted mean uses w (= weights^2 × μ_w), not sqrt(w).
      wmean <- drop(crossprod(weights^2, yunweighted) / sum(weights^2))
      # Grouped exactly as `sweep(...)^2 * weights^2` was, so the sum comes
      # out bit-identical; `sweep()`'s array-plus-aperm pair is what goes.
      ymw <- yunweighted - rep(wmean, each = N)
      tss <- colSums((ymw * ymw) * (weights * weights)) * weight_mean
    } else {
      tss <- colSums(y^2 * weight_mean)
    }
  } else {
    if (has_int) {
      # Centre and square directly. `sweep()` builds the subtrahend with
      # `array()` and then `aperm()`s it, so it materialises two extra n-by-k
      # copies to do what recycling does in one.
      ym <- y - rep(colMeans(y), each = N)
      tss <- colSums(ym * ym)
    } else {
      tss <- colSums(y^2)
    }
  }

  tss <- as.vector(tss)

  r.squared <-
    1 - (
      return_list[["df.residual"]] * return_list[["res_var"]] /
        tss
    )

  # `N` above is the row count, which the centring needs; the adjustment
  # counts observations, and a zero-weight row is not one.
  adj.r.squared <-
    1 - (
      (1 - r.squared) *
        ((return_list[["nobs"]] - has_int) / return_list[["df.residual"]])
    )

  return(list(
    tss = tss,
    r.squared = r.squared,
    adj.r.squared = adj.r.squared
  ))
}

get_fstat <- function(tss_r2s,
                      return_list,
                      iv_ei,
                      nomdf,
                      dendf,
                      vcov_fit,
                      has_int,
                      iv_stage,
                      no_variation = logical(0)) {

  coefs <- as.matrix(return_list$coefficients)

  if (length(return_list[["outcome"]]) > 1) {
    fstat_names <- paste0(return_list[["outcome"]], ":value")
  } else {
    fstat_names <- "value"
  }

  if (iv_stage[[1]] != 2 && return_list[["se_type"]] == "classical") {
    fstat <- tss_r2s$r.squared * return_list[["df.residual"]] /
        ((1 - tss_r2s$r.squared) * (nomdf))
  } else if (return_list[["se_type"]] == "classical" &&
             iv_stage[[1]] == 2 &&
             !return_list[["weighted"]]) {
    ivrss <- colSums(iv_ei^2)
    fstat <- ((tss_r2s$tss - ivrss) / nomdf) / return_list[["res_var"]]
  } else {
    # `coefs` is the full-length vector, carrying an NA wherever a column was
    # dropped, while `indices` and the variance matrix both count positions in
    # the kept set. Where the dropped column is the last one the two agree, and
    # where it is not they do not: `y ~ x1 + mid + en` with `mid` collinear
    # read the NA and every robust se_type reported the F statistic as NA,
    # where lm() and this function's classical branch both return it. Dropping
    # the NA rows first puts `coefs` on the same footing as the variance.
    coefs <- coefs[!is.na(coefs[, 1]), , drop = FALSE]
    indices <-
      seq.int(has_int + 1, return_list[["rank"]], by = 1)

    fstat <- compute_fstat(
      coef_matrix = coefs,
      coef_indices = indices,
      vcov_fit = vcov_fit$Vcov_hat,
      rank = return_list[["rank"]],
      nomdf = nomdf
    )

  }

  # A coefficient whose variance is numerically zero makes the joint test over
  # any set containing it as meaningless as its own standard error, and the
  # classical branch reaches the same place through an R-squared of exactly 1.
  # Reporting NA for one and a number for the other is the inconsistency this
  # rule exists to remove.
  rank <- return_list[["rank"]]
  if (length(no_variation) == rank * length(fstat)) {
    tested <- seq.int(has_int + 1, rank, by = 1)
    for (i in seq_along(fstat)) {
      if (any(no_variation[tested + (i - 1) * rank])) fstat[i] <- NA_real_
    }
  }

  f <- c(
    setNames(fstat, fstat_names),
    numdf = nomdf,
    dendf = dendf
  )

  return(f)
}

compute_fstat <- function(coef_matrix, coef_indices, vcov_fit, rank, nomdf) {

  fstat <- numeric(ncol(coef_matrix))

  for (i in seq_along(fstat)) {
    vcov_indices <- coef_indices + (i - 1) * rank
    fstat[i] <- tryCatch(
      {
        V <- vcov_fit[vcov_indices, vcov_indices, drop = FALSE]
        # `chol()` errors only on a matrix it can see is not positive
        # definite, and a variance block at rcond 1e-32 is not one: the
        # factorization succeeds and `chol2inv()` returns an inverse built
        # from cancellation, which has given an F of 389 on a covariate
        # balance test where the contrast has no variance at all. A matrix
        # whose reciprocal condition number is below the double precision
        # floor has no numerical inverse, so there is no statistic to report.
        if (rcond(V) < .Machine$double.eps) {
          NA_real_
        } else {
          crossprod(
            coef_matrix[coef_indices, i],
            chol2inv(chol(V)) %*% coef_matrix[coef_indices, i]
          ) / nomdf
        }
      },
      error = function(e) {
        NA_real_
      }
    )
  }

  fstat
}

prep_data <- function(data,
                      se_type,
                      clustered,
                      weighted,
                      iv_stage) {

  if (clustered && se_type != "none") {
    data[["cl_ord"]] <- order(data[["cluster"]])
    data[["cluster"]] <- data[["cluster"]][data[["cl_ord"]]]
    data[["y"]] <- data[["y"]][data[["cl_ord"]], , drop = FALSE]
    data[["X"]] <- data[["X"]][data[["cl_ord"]], , drop = FALSE]
    if (!is.null(data[["femat"]])) {
      data[["femat"]] <- data[["femat"]][data[["cl_ord"]], , drop = FALSE]
    }

    if (weighted) {
      data[["weights"]] <- data[["weights"]][data[["cl_ord"]]]
    }
    if (iv_stage[[1]] == 2) {
      data[["X_first_stage"]] <-
        data[["X_first_stage"]][data[["cl_ord"]], , drop = FALSE]
    }

    data[["J"]] <- length(unique(data[["cluster"]]))
    # A cluster whose every row has weight 0 is not a cluster: it contributes
    # nothing to the fit or the meat, so it must not count toward the cluster
    # degrees of freedom, the stata factor, or `nclusters`, any more than a
    # zero-weight row counts toward `nobs`. `J` stays the number of distinct
    # codes because the C++ sizes its per-cluster arrays by it.
    data[["J_eff"]] <- if (weighted) {
      length(unique(data[["cluster"]][data[["weights"]] > 0]))
    } else {
      data[["J"]]
    }

    # Every cluster-robust estimator here divides by J - 1 somewhere, except
    # CR2, whose degrees of freedom are Satterthwaite and so never hit the
    # J - 1 = 0 guard. With one cluster CR2 returned standard errors of order
    # 1e-17 and said nothing, which reads as a precisely estimated zero.
    if (data[["J_eff"]] < 2L) {
      stop(
        "`clusters` has only one level, so there is no between-cluster ",
        "variation for a cluster-robust variance to estimate.\nDrop ",
        "`clusters`, or use `se_type = \"classical\"` if the clustering is ",
        "not the level you meant to inference over."
      )
    }
  } else {
    data[["J"]] <- 1
    data[["J_eff"]] <- 1
  }

  if (weighted) {
    data[["Xunweighted"]] <- data[["X"]]
    data[["yunweighted"]] <- data[["y"]]
    data[["weight_mean"]] <- mean(data[["weights"]])
    data[["weights"]] <- sqrt(data[["weights"]] / data[["weight_mean"]])
    data[["X"]] <- data[["weights"]] * data[["X"]]
    data[["y"]] <- data[["weights"]] * data[["y"]]
    if (!is.null(data[["femat"]])) {
      data[["fematunweighted"]] <- data[["femat"]]
      data[["femat"]] <- data[["weights"]] * data[["femat"]]
    }
    if (iv_stage[[1]] == 2) {
      data[["X_first_stage_unweighted"]] <- data[["X_first_stage"]]
      data[["X_first_stage"]] <- data[["weights"]] * data[["X_first_stage"]]
    }

  } else {
    data[["weight_mean"]] <- 1
  }

  return(data)
}

drop_collinear <- function(data, covs_used, weighted, iv_stage) {
  data[["X"]] <- data[["X"]][, covs_used, drop = FALSE]
  if (weighted) {
    data[["Xunweighted"]] <- data[["Xunweighted"]][, covs_used, drop = FALSE]
  }
  if (iv_stage[[1]] == 2) {
    data[["X_first_stage"]] <- data[["X_first_stage"]][, covs_used, drop = FALSE]
    if (weighted) {
      data[["X_first_stage_unweighted"]] <-
        data[["X_first_stage_unweighted"]][, covs_used, drop = FALSE]
    }
  }

  return(data)
}
