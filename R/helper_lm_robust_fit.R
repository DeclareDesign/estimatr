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
                          yoriginal = NULL,
                          Xoriginal = NULL,
                          weights,
                          cluster,
                          fixed_effects = NULL,
                          ci = TRUE,
                          se_type,
                          has_int, # TODO get this out of here
                          alpha = 0.05,
                          return_vcov = TRUE,
                          return_fit = TRUE,
                          return_unweighted_fit = TRUE,
                          try_cholesky = FALSE,
                          X_first_stage = NULL) {

  y <- as.matrix(y)
  ny <- ncol(y)
  fes <- !is.null(fixed_effects)
  if (fes) {
    if (is.numeric(yoriginal)) {
      yoriginal <- as.matrix(yoriginal)
    }
    if (is.numeric(Xoriginal)) {
      Xoriginal <- as.matrix(Xoriginal)
    }
    fe_rank <- attr(fixed_effects, "fe_rank")
  } else {
    fe_rank <- 0
  }
  multivariate <- ny > 1
  weighted <- !is.null(weights)
  iv_second_stage <- !is.null(X_first_stage)
  clustered <- !is.null(cluster)

  # ----------
  # Check se type
  # ----------

  se_type <- check_se_type(se_type, clustered, iv_second_stage)

  if (weighted && se_type == "CR2" && fes) {
    stop(
      "Cannot use `fixed_effects` with weighted CR2 estimation at the moment. ",
      "Try setting `se_type` = \"stata\""
    )
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

  ## TODO build list with all data and have weight and order helper functions
  # Reorder if there are clusters and you need the SE or to return the fit
  if (clustered && se_type != "none") {
    cl_ord <- order(cluster)
    cluster <- cluster[cl_ord]
    y <- y[cl_ord, , drop = FALSE]
    X <- X[cl_ord, , drop = FALSE]

    if (fes) {
      fixed_effects <- fixed_effects[cl_ord, , drop = FALSE]
      yoriginal <- yoriginal[cl_ord, , drop = FALSE]
      Xoriginal <- Xoriginal[cl_ord, , drop = FALSE]
    }
    if (weighted) {
      weights <- weights[cl_ord]
    }
    if (iv_second_stage) {
      X_first_stage <- X_first_stage[cl_ord, , drop = FALSE]
    }

    J <- length(unique(cluster))
  } else {
    J <- 1
  }

  if (fes) {
    femat <- model.matrix(~ 0 + ., data = as.data.frame(fixed_effects))
  }

  # Weight if there are weights
  if (weighted) {
    Xunweighted <- X
    yunweighted <- y
    weight_mean <- mean(weights)
    weights <- sqrt(weights / weight_mean)
    X <- weights * X
    y <- weights * y
    if (fes) {
      if (is.numeric(yoriginal)) {
        yoriginalunweighted <- yoriginal
        yoriginal <- weights * yoriginal
      }
      fematunweighted <- femat
      femat <- weights * femat
    }
    if (iv_second_stage) {
      X_first_stage_unweighted <- X_first_stage
      X_first_stage <- weights * X_first_stage
    }
  } else {
    weight_mean <- 1
    yunweighted <- NULL
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
  covs_used <- which(est_exists[, 1])
  N <- nrow(X)

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

    # Drop NA columns from data and from beta_hat
    if (x_rank < ncol(X)) {
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
      if (is.numeric(Xoriginal)) {
        Xoriginal <- Xoriginal[, covs_used, drop = FALSE]
        if (weighted) {
          Xoriginal <- Xoriginal[, covs_used, drop = FALSE]
        }
      }

      fit$beta_hat <- fit$beta_hat[covs_used, , drop = FALSE]
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
      ei_weight <- as.matrix(y - X[, 1:x_rank, drop = FALSE] %*% fit$beta_hat)
      X <- weights * X
      if (fes) {
        femat <- weights * femat
      }
    } else {
      if (iv_second_stage) {
        fitted.values <- X_first_stage %*% fit$beta_hat
      } else {
        fitted.values <- X %*% fit$beta_hat
      }
      ei <- as.matrix(y - fitted.values)
    }

    if (iv_second_stage) {
      iv_fits <- X %*% fit$beta_hat
      if (weighted) {
        iv_ei <- weights * as.matrix(y - iv_fits)
      } else {
        iv_ei <- as.matrix(y - iv_fits)
      }
    } else {
      iv_fits <- NULL
      iv_ei <- NULL
    }

    if (se_type != "none") {

      vcov_fit <- lm_variance(
        X = if (se_type %in% c("HC2", "HC3", "CR2") && fes) cbind(X, femat) else X,
        Xunweighted = if (se_type %in% c("HC2", "HC3", "CR2") && fes && weighted) cbind(Xunweighted, fematunweighted) else Xunweighted,
        XtX_inv = fit$XtX_inv,
        ei = ei,
        weight_mean = weight_mean,
        cluster = cluster,
        J = J,
        ci = ci,
        se_type = se_type,
        which_covs = which_covs[covs_used],
        fe_rank = fe_rank
      )

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
    if (iv_second_stage && !fes) {
      if (weighted) {
        return_list[["fitted.values"]] <- as.matrix(X_first_stage_unweighted %*% fit$beta_hat)
      } else {
        return_list[["fitted.values"]] <- as.matrix(fitted.values)
      }
    } else if ((se_type == "CR2" && weighted)) {
      # Have to get weighted fits as original fits were unweighted for
      # variance estimation or used wrong regressors in IV
      return_list[["fitted.values"]] <- as.matrix(X %*% fit$beta_hat)
      if (fes) {
        return_list[["fitted.values"]] <- as.matrix(yoriginal - (y - return_list[["fitted.values"]]))
      }
    } else if (weighted && return_unweighted_fit) {
      if (fes && is.numeric(yoriginal)) {
        return_list[["fitted.values"]] <- as.matrix(yoriginalunweighted - ei / weights)
      } else {
        return_list[["fitted.values"]] <- as.matrix(Xunweighted %*% fit$beta_hat)
      }
    } else {
      if (fes && is.numeric(yoriginal)) {
        return_list[["fitted.values"]] <- as.matrix(yoriginal - ei)
      } else {
        return_list[["fitted.values"]] <- as.matrix(fitted.values)
      }
    }

    if (fes && (ncol(fixed_effects) == 1)) {
      return_list[["fixed_effects"]] <- setNames(
        tapply(
          return_list[["fitted.values"]] -
            Xoriginal[, variable_names, drop = FALSE] %*% fit$beta_hat,
          fixed_effects,
          `[`,
          1
        ),
        colnames(femat)
      )
    }

    # If we reordered to get SEs earlier, have to fix order
    if (clustered && se_type != "none") {
      return_list[["fitted.values"]] <- return_list[["fitted.values"]][order(cl_ord), , drop = FALSE]
    }

    colnames(return_list[["fitted.values"]]) <- colnames(y)
  }

  return_list[["term"]] <- variable_names
  return_list[["outcome"]] <- colnames(y)
  return_list[["alpha"]] <- alpha
  return_list[["se_type"]] <- se_type
  return_list[["weighted"]] <- weighted
  return_list[["fes"]] <- fes
  return_list[["clustered"]] <- clustered
  return_list[["df.residual"]] <- N - tot_rank
  return_list[["N"]] <- N
  return_list[["k"]] <- k
  return_list[["rank"]] <- x_rank

  if (se_type != "none") {

    if (se_type == "CR2" && weighted) {
      ei <- ei_weight
    }

    if (weighted) {
      return_list[["res_var"]] <-
        colSums(ei^2 * weight_mean) / (return_list[["df.residual"]])
    } else {
      return_list[["res_var"]] <-
        diag(as.matrix(ifelse(vcov_fit[["res_var"]] < 0, NA, vcov_fit[["res_var"]])))
    }

    tss_r2s <- get_r2s(
      y = y,
      return_list = return_list,
      has_int = has_int,
      yunweighted = yunweighted,
      weights = weights,
      weight_mean = weight_mean
    )

    nomdf <- x_rank - as.numeric(!fes) * has_int
    if (clustered) {
      dendf <- J - 1
    } else {
      dendf <- return_list[["df.residual"]]
    }

    if (nomdf > 0) {
      f <- get_fstat(
        tss_r2s = tss_r2s,
        return_list = return_list,
        iv_ei = iv_ei,
        nomdf = nomdf,
        dendf = dendf,
        fit = fit,
        vcov_fit = vcov_fit,
        has_int = has_int,
        iv_second_stage = iv_second_stage
      )
    } else {
      f <- NULL
    }

    if (!fes) {
      return_list <- c(return_list, tss_r2s)
      return_list[["fstatistic"]] <- f
    } else {
      return_list <- c(return_list, setNames(tss_r2s, paste0("proj_", names(tss_r2s))))
      return_list[["proj_fstatistic"]] <- f

      tss_r2s <- get_r2s(
        y = yoriginal,
        return_list = return_list,
        has_int = has_int,
        yunweighted = yoriginalunweighted,
        weights = weights,
        weight_mean = weight_mean
      )

      # nomdf <- tot_rank - has_int
      # f <- get_fstat(
      #   tss_r2s = tss_r2s,
      #   return_list = return_list,
      #   iv_ei = iv_ei,
      #   nomdf = nomdf,
      #   dendf = dendf,
      #   fit = fit,
      #   vcov_fit = vcov_fit,
      #   has_int = has_int,
      #   iv_second_stage = iv_second_stage
      # )

      return_list <- c(return_list, tss_r2s)
      # return_list[["fstatistic"]] <- f

    }

    if (return_vcov) {
      # return_list$residuals <- fit$residuals
      return_list[["vcov"]] <- vcov_fit$Vcov_hat
      if (multivariate) {
        coef_names <- lapply(seq_len(ncol(est_exists)), function(j) return_list$term[est_exists[, j]])

        outcome_coef_names <- paste0(
          rep(paste0(return_list[["outcome"]], ":"), times = sapply(coef_names, length)),
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

check_se_type <- function(se_type, clustered, iv) {

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

  return(se_type)
}

get_r2s <- function(y, return_list, has_int, yunweighted, weights, weight_mean, iv_second_stage) {

  N <- nrow(y)
  if (return_list[["weighted"]]) {
    if (has_int) {
      tss <-
        colSums(apply(yunweighted, 2, function(x) {
          weights^2 * (x - weighted.mean(x, weights^2))^2
        })) * weight_mean
    } else {
      tss <- colSums(y^2 * weight_mean)
    }
  } else {
    if (has_int) {
      tss <- .rowSums(apply(y, 1, `-`, colMeans(y))^2, ncol(y), N)
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

get_fstat <- function(tss_r2s, return_list, iv_ei, nomdf, dendf, fit, vcov_fit, has_int, iv_second_stage) {

  if (length(return_list[["outcome"]]) > 1) {
    fstat_names <- paste0(return_list[["outcome"]], ":value")
  } else {
    fstat_names <- "value"
  }

  if (!iv_second_stage && return_list[["se_type"]] == "classical") {
    fstat <- tss_r2s$r.squared * return_list[["df.residual"]] /
        ((1 - tss_r2s$r.squared) * (nomdf))
  } else if (return_list[["se_type"]] == "classical" && iv_second_stage && !return_list[["weighted"]]) {
    ivrss <- colSums(iv_ei^2)
    fstat <- ((tss_r2s$tss - ivrss) / nomdf) / return_list[["res_var"]]
  } else {
    indices <- seq.int(has_int + (!return_list[["fes"]]), return_list[["rank"]], by = 1)
    fstat <-
      tryCatch({
        sapply(seq_len(ncol(fit$beta_hat)),
               function(x) {
                 vcov_indices <- indices + (x - 1) * return_list[["rank"]]
                 crossprod(fit$beta_hat[indices, x],
                           chol2inv(chol(vcov_fit$Vcov_hat[vcov_indices, vcov_indices])) %*%
                             fit$beta_hat[indices, x]) / nomdf
               })
      }, error = function(e) {
        # warning("Unable to compute f-statistic because variance-covariance matrix cannot be inverted")
        rep(NA, length(fstat_names))
      })
  }

  f <- c(
    setNames(fstat, fstat_names),
    numdf = nomdf,
    dendf = dendf
  )

  return(f)
}
