#' Internal method that creates linear fits (no fixed effects)
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
#' @param try_cholesky logical, whether to try Cholesky decomposition
#' @param iv_stage list of length one or two for 2SLS stages
#' @param fe_rank integer, degrees of freedom absorbed by fixed effects
#' @param femat optional numeric matrix of fixed-effect dummies for the
#'   estimation sample, supplying the columns HC2, HC3 and CR2 need from the
#'   full design. `NULL` unless the requested `se_type` requires it.
#' @param fe_leverage numeric vector of per-observation leverage contributed by
#'   the absorbed fixed effects, or `NULL`. `h_ii` of the full design splits
#'   exactly into the demeaned-X leverage plus this term, for any number of FE
#'   factors, which is what makes HC2 and HC3 available under `fixed_effects`
#'   without building the dummy matrix.
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
                          femat = NULL) {

  # ----------
  # Check se type
  # ----------

  clustered <- !is.null(cluster)
  weighted <- !is.null(weights)
  se_type <- check_se_type(se_type, clustered, has_fe = (fe_rank > 0L),
                           has_fe_leverage = !is.null(fe_leverage),
                           fe_dummies = !is.null(femat), weighted = weighted)

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
    data[["cluster"]] <- cluster
  }

  k <- ncol(data[["X"]])

  if (is.null(colnames(data[["X"]]))) {
    colnames(data[["X"]]) <- paste0("X", 1:k)
  }
  variable_names <- colnames(data[["X"]])

  which_covs <- setNames(rep(TRUE, k), variable_names)

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
  N <- nrow(data[["X"]])

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

      if (se_type == "CR2") {
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

      vcov_fit <- lm_variance(
        # Widening the design here is what makes HC2, HC3 and CR2 available
        # under `fixed_effects`. The C++ detects X.cols() > ncol(XtX_inv) and
        # rebuilds the meat from the full design, dropping the rank-deficient
        # dummy columns itself.
        X = if (!is.null(data[["femat"]]))
          cbind(data[["X"]], data[["femat"]])
          else data[["X"]],
        Xunweighted = if (!is.null(data[["femat"]]) && weighted)
          cbind(data[["Xunweighted"]], data[["fematunweighted"]])
          else data[["Xunweighted"]],
        XtX_inv = fit$XtX_inv,
        ei = if (se_type == "CR2" && weighted)
          fit_vals[["ei.unweighted"]]
          else fit_vals[["ei"]],
        weight_mean = data[["weight_mean"]],
        cluster = data[["cluster"]],
        J = data[["J"]],
        ci = ci,
        se_type = se_type,
        which_covs = which_covs[covs_used],
        fe_rank = fe_rank,
        # Only HC2/HC3 consume this; it is NULL for every other se_type and for
        # multi-way FE, so the C++ falls back to the plain hat value.
        fe_leverage = if (se_type %in% c("HC2", "HC3")) fe_leverage else NULL
      )

      return_list$std.error[est_exists] <- sqrt(diag(vcov_fit$Vcov_hat))

      # HC2 and HC3 divide by (1 - h_ii). A near-saturated design produces
      # observations that are fitted exactly, whose computed hat value can land
      # marginally above 1; lm_variance() drops those from the meat rather than
      # divide by a negative number, and this is where that gets said. It is
      # worth saying: those observations contribute nothing, so the standard
      # error is built from fewer rows than the fit used (estimatr #395).
      n_lev <- vcov_fit[["n_leverage_above_one"]]
      if (isTRUE(n_lev > 0)) {
        warning(
          n_lev, if (n_lev == 1) " observation has " else " observations have ",
          "a computed leverage above 1, which happens when the design is close ",
          "to saturated and the observation is fitted exactly. `se_type = \"",
          se_type, "\"` divides by (1 - leverage), so those observations are ",
          "dropped from the variance rather than divided by a negative number. ",
          "Use `se_type = \"HC1\"` or `\"classical\"`, or drop covariates, to ",
          "use every observation."
        )
      } else if (any(is.nan(return_list$std.error)) &&
                 se_type %in% c("HC2", "HC3", "CR2")) {
        warning(
          "Some standard errors are NaN. `se_type = \"", se_type, "\"` divides ",
          "by the observation's leverage, which is at or near 1 for some ",
          "observations here, as happens when the design is close to ",
          "saturated. Use `se_type = \"HC1\"` or `\"classical\"`, or drop ",
          "covariates, to get finite standard errors."
        )
      }

      if (ci) {
        return_list$df[est_exists] <-
          ifelse(vcov_fit$dof == -99, NA, vcov_fit$dof)
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
    return_list[["nclusters"]] <- data[["J"]]
  }
  return_list[["k"]] <- k
  return_list[["rank"]] <- x_rank

  if (se_type != "none") {

    return_list[["res_var"]] <- get_resvar(
      data = data,
      ei = fit_vals[["ei"]],
      df.residual = return_list[["df.residual"]],
      vcov_fit = vcov_fit,
      weighted = weighted
    )

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
      dendf <- data[["J"]] - 1
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
        iv_stage = iv_stage
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

  cl_se_types  <- c("CR0", "CR2", "stata")
  rob_se_types <- c("HC0", "HC1", "HC2", "HC3", "classical", "stata")

  # HC2 and HC3 need the hat values of the full [dummies | X] design. Those
  # split exactly into the demeaned-X leverage plus fe_leverage(), for ANY
  # number of fixed-effect factors, so they are exact and cheap here and carry
  # no restriction. The identity is stated and verified in fe_leverage().
  if (has_fe && has_fe_leverage && !is.null(se_type) &&
      se_type %in% c("HC2", "HC3")) {
    return(se_type)
  }

  # CR2 is the one case left that has to materialise the dummies: it needs the
  # whole per-cluster block of the hat matrix, which does not decompose.
  if (has_fe && !is.null(se_type) && se_type %in% c("HC2", "HC3", "CR2")) {
    # 1.0.6 refused this combination too: the CR2 adjustment needs the design
    # weighted by w rather than by sqrt(w), and the absorbed dummies are not
    # carried in a form that supports it.
    if (se_type == "CR2" && weighted) {
      stop(
        "`se_type = \"CR2\"` cannot be combined with `weights` and ",
        "`fixed_effects`.\nUse `se_type = \"stata\"`, or replace ",
        "`fixed_effects` with explicit dummies:\n",
        "  lm_robust(y ~ x + factor(fe_var), weights = w, clusters = cl, ",
        "se_type = \"CR2\")"
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
    return(se_type)
  }

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
        "`se_type` must be either 'CR0', 'stata', 'CR2', or 'none' when ",
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
    } else if (se_type %in% setdiff(cl_se_types, "stata")) {
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
    if (se_type == "stata") se_type <- "HC1"
  }

  return(se_type)
}

get_resvar <- function(data, ei, df.residual, vcov_fit, weighted) {
  res_var <-
    if (weighted)
      colSums(ei^2 * data[["weight_mean"]]) / df.residual
    else
      as.vector(ifelse(vcov_fit[["res_var"]] < 0, NA, vcov_fit[["res_var"]]))
  return(res_var)
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

  adj.r.squared <-
    1 - (
      (1 - r.squared) *
        ((N - has_int) / return_list[["df.residual"]])
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
                      iv_stage) {

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
        crossprod(
          coef_matrix[coef_indices, i],
          chol2inv(chol(vcov_fit[vcov_indices, vcov_indices])) %*%
            coef_matrix[coef_indices, i]
        ) / nomdf
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
  } else {
    data[["J"]] <- 1
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
