#' Internal method that creates linear fits
#'
#' @param y numeric outcome vector
#' @param X numeric design matrix
#' @param weights numeric weights vector
#' @param cluster numeric cluster vector
#' @param ci boolean that when T returns confidence intervals and p-values
#' @param se_type character denoting which kind of SEs to return
#' @param has_int logical, whether the model has an intercept, used for \eqn{R^2}
#' @param alpha numeric denoting the test size for confidence intervals
#' @param return_vcov logical, whether to return the vcov matrix for later usage
#' @param return_fit logical, whether to return fitted values
#' @param return_unweighted_fit logical, whether to return unweighted fitted values in place of weighted fitted values if regression is weighted
#' @param try_cholesky logical, whether to try using a cholesky decomposition to solve LS instead of a QR decomposition
#' @param X_first_stage numeric matrix of the first stage design matrix, only use for second stage of 2SLS IV regression, otherwise leave as \code{NULL}
#'
#' @export
#'
lm_robust_fit <- function(y,
                          X,
                          weights,
                          cluster,
                          ci,
                          se_type,
                          has_int, # TODO get this out of here
                          alpha = 0.05,
                          return_vcov = TRUE,
                          return_fit = FALSE,
                          return_unweighted_fit = FALSE,
                          try_cholesky = FALSE,
                          X_first_stage = NULL) {
  y <- as.matrix(y)
  ny <- ncol(y)
  multivariate <- ny > 1
  weighted <- !is.null(weights)
  iv_second_stage <- !is.null(X_first_stage)
  clustered <- !is.null(cluster)

  # ----------
  # Check se type
  # ----------

  # Allowable se_types with clustering
  cl_se_types <- c("CR0", "CR2", "stata")
  rob_se_types <- c("HC0", "HC1", "HC2", "HC3", "classical", "stata")

  # Parse cluster variable
  if (clustered) {

    # set/check se_type
    if (is.null(se_type)) {
      se_type <- "CR2"
    } else if (!(se_type %in% c(cl_se_types, "none"))) {
      stop(
        "`se_type` must be either 'CR0', 'stata', 'CR2', or 'none' when ",
        "`clusters` are specified.\nYou passed: ", se_type
      )
    }
  } else {

    # set/check se_type
    if (is.null(se_type)) {
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
    } else if (se_type == "stata") {
      se_type <- "HC1" # In IV this is true with small option
    }
  }

  # -----------
  # Prep data for fitting
  # -----------

  k <- ncol(X)

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", 1:k)
  }
  variable_names <- colnames(X)

  # Legacy, in case we want to only get some covs in the future
  which_covs <- rep(TRUE, ncol(X))

  # Reorder if there are clusters and you need the SE or to return the fit
  if (clustered && se_type != "none") {
    cl_ord <- order(cluster)
    y <- y[cl_ord, , drop = FALSE]
    X <- X[cl_ord, , drop = FALSE]
    cluster <- cluster[cl_ord]
    J <- length(unique(cluster))
    if (weighted) {
      weights <- weights[cl_ord]
    }
    if (iv_second_stage) {
      X_first_stage <- X_first_stage[cl_ord, , drop = FALSE]
    }
  } else {
    J <- 1
  }

  # Weight if there are weights
  if (weighted) {
    Xunweighted <- X
    yunweighted <- y
    weight_mean <- mean(weights)
    weights <- sqrt(weights / weight_mean)
    X <- weights * X
    y <- weights * y
    if (iv_second_stage) {
      X_first_stage_unweighted <- X_first_stage
      X_first_stage <- weights * X_first_stage
    }
  } else {
    weight_mean <- 1
    Xunweighted <- NULL
  }

  # -----------
  # Estimate coefficients
  # -----------

  fit <-
    lm_solver(
      X = X,
      y = y,
      try_cholesky = try_cholesky
    )

  fit$beta_hat <- as.matrix(fit$beta_hat)
  dimnames(fit$beta_hat) <- list(variable_names, colnames(y))

  # Use first model to get linear dependencies
  est_exists <- !is.na(fit$beta_hat)
  covs_used <- est_exists[, 1]
  N <- nrow(X)
  rank <- sum(covs_used)

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
    if (rank < ncol(X)) {
      X <- X[, covs_used, drop = FALSE]
      if (weighted) {
        Xunweighted <- Xunweighted[, covs_used, drop = FALSE]
      }
      if (iv_second_stage) {
        X_first_stage <- X_first_stage[, covs_used, drop = FALSE]
        if (weighted) {
          X_first_stage_unweighted <- X_first_stage_unweighted[, covs_used, drop = FALSE]
        }
      }

      fit$beta_hat <- fit$beta_hat[covs_used, ]
    }

    # compute fitted.values and residuals
    # need unweighted for CR2, as well as X weighted by weights again
    # so that instead of having X * sqrt(W) we have X * W
    if (se_type == "CR2" && weighted) {
      if (iv_second_stage) {
        fitted.values <- X_first_stage_unweighted %*% fit$beta_hat
      } else {
        fitted.values <- Xunweighted %*% fit$beta_hat
      }
      ei <- as.matrix(yunweighted - fitted.values)
      X <- weights * X
    } else {
      if (iv_second_stage) {
        fitted.values <- X_first_stage %*% fit$beta_hat
      } else {
        fitted.values <- X %*% fit$beta_hat
      }
      ei <- as.matrix(y - fitted.values)
    }

    if (se_type != "none") {
      if (se_type == "CR2") {
        vcov_fit <- lm_variance_cr2(
          X = X,
          Xunweighted = Xunweighted,
          XtX_inv = fit$XtX_inv,
          ei = ei,
          weight_mean = weight_mean,
          clusters = cluster,
          J = J,
          ci = ci,
          which_covs = which_covs[covs_used]
        )
        vcov_fit[["res_var"]] <-
          colSums((y - X %*% fit$beta_hat)^2) /
            (N - rank)
      } else {
        vcov_fit <- lm_variance(
          X = X,
          XtX_inv = fit$XtX_inv,
          ei = ei,
          cluster = cluster,
          J = J,
          ci = ci,
          type = se_type,
          which_covs = which_covs[covs_used]
        )
      }
      # print(est_exists)
      # print(vcov_fit)
      return_list$std.error[est_exists] <- sqrt(diag(vcov_fit$Vcov_hat))

      if (ci) {
        # If any not computed in variance fn, replace with NA
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
    if ((se_type == "CR2" && weighted) || iv_second_stage) {
      # Have to get weighted fits as original fits were unweighted for
      # variance estimation or used wrong regressors in IV
      return_list[["fitted.values"]] <- as.matrix(X %*% fit$beta_hat)
    } else {
      return_list[["fitted.values"]] <- as.matrix(fitted.values)
    }

    if (weighted && return_unweighted_fit) {
      return_list[["fitted.values"]] <- as.matrix(Xunweighted %*% fit$beta_hat)
    }

    # If we reordered to get SEs earlier, have to fix order
    if (clustered && se_type != "none") {
      return_list[["fitted.values"]] <- return_list[["fitted.values"]][order(cl_ord), ]
    }
  }

  return_list[["term"]] <- variable_names
  return_list[["outcome"]] <- colnames(y)
  return_list[["alpha"]] <- alpha
  return_list[["se_type"]] <- se_type
  return_list[["weighted"]] <- weighted

  # return_list[["fitted.values"]] <- fit$fit
  # return_list[["residuals"]] <- fit$residuals

  return_list[["df.residual"]] <- N - rank

  # return_list[["XtX_inv"]] <- fit$XtX_inv
  return_list[["N"]] <- N
  return_list[["k"]] <- k
  return_list[["rank"]] <- rank

  if (se_type != "none") {
    if (weighted) {
      if (has_int) {
        return_list[["tot_var"]] <-
          colSums(
            weights^2 *
              (yunweighted - weighted.mean(yunweighted, weights^2))^2
          ) * weight_mean
      } else {
        return_list[["tot_var"]] <- colSums(y^2 * weight_mean)
      }
      return_list[["res_var"]] <- diag(as.matrix(colSums(ei^2 * weight_mean) / (N - rank)))
    } else {
      if (has_int) {
        return_list[["tot_var"]] <- .rowSums(apply(y, 1, `-`, colMeans(y))^2, ny, N)
      } else {
        return_list[["tot_var"]] <- colSums(y^2)
      }
      return_list[["res_var"]] <- diag(as.matrix(ifelse(vcov_fit$res_var < 0, NA, vcov_fit$res_var)))
    }
    return_list[["tot_var"]] <- as.vector(return_list[["tot_var"]])

    return_list[["r.squared"]] <-
      1 - (
        return_list[["df.residual"]] * return_list[["res_var"]] /
          return_list[["tot_var"]]
      )

    return_list[["adj.r.squared"]] <-
      1 - (
        (1 - return_list[["r.squared"]]) *
          ((N - has_int) / return_list[["df.residual"]])
      )

    if (ny > 1) {
      fstat_names <- paste0(colnames(y), ":value")
    } else {
      fstat_names <- "value"
    }
    if (iv_second_stage && se_type != "none") {
      indices <- seq.int(has_int + 1, rank, by = 1)
      fstat <- setNames(
        crossprod(
          fit$beta_hat[indices],
          solve(vcov_fit$Vcov_hat[indices, indices], fit$beta_hat[indices])
        ) / (rank - has_int),
        fstat_names
      )
    } else {
      fstat <- setNames(
        return_list[["r.squared"]] * return_list[["df.residual"]] /
          ((1 - return_list[["r.squared"]]) * (rank - has_int)),
        fstat_names
      )
    }
    return_list[["fstatistic"]] <- c(
      fstat,
      numdf = rank - has_int,
      dendf = return_list[["df.residual"]]
    )

    if (return_vcov) {
      # return_list$residuals <- fit$residuals
      return_list[["vcov"]] <- vcov_fit$Vcov_hat
      if (multivariate) {
        coef_names <- paste0(
          rep(paste0(return_list[["outcome"]], ":"), each = rank),
          rep(return_list$term, times = ny)
        )
        # print(return_list[["vcov"]])
        # print(coef_names)
        dimnames(return_list[["vcov"]]) <- list(
          coef_names,
          coef_names
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
