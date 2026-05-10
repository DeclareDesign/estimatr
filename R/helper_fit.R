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
                          fe_rank = 0L) {

  # ----------
  # Check se type
  # ----------

  clustered <- !is.null(cluster)
  weighted <- !is.null(weights)
  se_type <- check_se_type(se_type, clustered, has_fe = (fe_rank > 0L))

  # -----------
  # Prep data for fitting
  # -----------

  data <- list(
    y = as.matrix(y),
    X = X
  )

  ny <- ncol(data[["y"]])
  ynames <- colnames(data[["y"]])
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

    fit_vals[["fitted.values"]] <- as.matrix(
      data[[X_name]][, seq_len(x_rank), drop = FALSE] %*% fit$beta_hat
    )

    fit_vals[["ei"]] <- as.matrix(data[["y"]] - fit_vals[["fitted.values"]])

    if (weighted) {
      fit_vals[["fitted.values.unweighted"]] <- as.matrix(
        data[[X_name_unweighted]] %*% fit$beta_hat
      )

      fit_vals[["ei.unweighted"]] <- as.matrix(
        data[["yunweighted"]] - fit_vals[["fitted.values.unweighted"]]
      )

      if (se_type == "CR2") {
        data[["X"]] <- data[["weights"]] * data[["X"]]
      }
    }

    if (iv_stage[[1]] == 2) {
      fit_vals[["fitted.values.iv"]] <- as.matrix(data[["X"]] %*% fit$beta_hat)
      fit_vals[["ei.iv"]] <- as.matrix(data[["y"]] - fit_vals[["fitted.values.iv"]])
      if (weighted) {
        fit_vals[["ei.iv"]] <- data[["weights"]] * fit_vals[["ei.iv"]]
      }
      return_list[["ei.iv"]] <- fit_vals[["ei.iv"]]
    }

    if (se_type != "none") {

      vcov_fit <- lm_variance(
        X = data[["X"]],
        Xunweighted = data[["Xunweighted"]],
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
        fe_rank = fe_rank
      )

      return_list$std.error[est_exists] <- sqrt(diag(vcov_fit$Vcov_hat))

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

    fitted.vals_name <- if (weighted)
      "fitted.values.unweighted"
      else "fitted.values"

    return_list[["fitted.values"]] <- as.matrix(fit_vals[[fitted.vals_name]])

    if (clustered && se_type != "none") {
      return_list[["fitted.values"]] <- return_list[["fitted.values"]][order(data[["cl_ord"]]), , drop = FALSE]
    }

    colnames(return_list[["fitted.values"]]) <- ynames
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

    nomdf <- x_rank - has_int
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
        has_int = has_int,
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

check_se_type <- function(se_type, clustered, has_fe = FALSE) {

  cl_se_types  <- c("CR0", "CR2", "stata")
  rob_se_types <- c("HC0", "HC1", "HC2", "HC3", "classical", "stata")

  if (clustered) {

    if (is.null(se_type)) {
      # With FE, CR2 requires the full augmented hat matrix; default to stata instead
      se_type <- if (has_fe) "stata" else "CR2"
    } else if (!(se_type %in% c(cl_se_types, "none"))) {
      stop(
        "`se_type` must be either 'CR0', 'stata', 'CR2', or 'none' when ",
        "`clusters` are specified.\nYou passed: ", se_type
      )
    }
  } else {

    if (is.null(se_type)) {
      # With FE, HC2/HC3 require the full augmented hat matrix; default to HC1 instead
      se_type <- if (has_fe) "HC1" else "HC2"
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
      se_type <- "HC1"
    }
  }

  # FE constraints: HC2/HC3/CR2 require the full augmented hat matrix over all
  # FE dummies, which estimatrZero does not form (feols-style demeaning only).
  # Warn and fall back to the FE-compatible equivalent.
  if (has_fe) {
    if (!clustered && se_type %in% c("HC2", "HC3")) {
      warning(
        "'", se_type, "' is not supported with `fixed_effects` in estimatrZero ",
        "(it requires forming the full hat matrix over FE dummies, which is ",
        "avoided for speed). Using 'HC1'. For exact HC2 with fixed effects, ",
        "see fixest::feols()."
      )
      se_type <- "HC1"
    }
    if (clustered && se_type == "CR2") {
      warning(
        "'CR2' is not supported with `fixed_effects` in estimatrZero. ",
        "Using 'stata' cluster-robust standard errors."
      )
      se_type <- "stata"
    }
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
      tss <- colSums(sweep(yunweighted, 2, wmean)^2 * weights^2) * weight_mean
    } else {
      tss <- colSums(y^2 * weight_mean)
    }
  } else {
    if (has_int) {
      # sweep replaces apply(y, 1, '-', colMeans(y)): ~80x faster for large n
      tss <- colSums(sweep(y, 2, colMeans(y))^2)
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
