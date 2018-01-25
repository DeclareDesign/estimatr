#' Internal method that creates linear fits
#'
#' @param y numeric outcome vector
#' @param X numeric design matrix
#' @param weights numeric weights vector
#' @param cluster numeric cluster vector
#' @param ci boolean that when T returns confidence intervals and p-values
#' @param se_type character denoting which kind of SEs to return
#' @param alpha numeric denoting the test size for confidence intervals
#' @param return_vcov logical, whether to return the vcov matrix for later usage
#' @param try_cholesky logical, whether to try using a cholesky decomposition to solve LS instead of a QR decomposition
#' @param has_int logical, whether the model has an intercept, used for \eqn{R^2}
#'
#' @export
#'
lm_robust_fit <- function(y,
                          X,
                          weights,
                          cluster,
                          ci,
                          se_type,
                          alpha,
                          return_vcov,
                          try_cholesky,
                          has_int) {

  ## allowable se_types with clustering
  cl_se_types <- c("CR0", "CR2", "stata")
  rob_se_types <- c("HC0", "HC1", "HC2", "HC3", "classical", "stata")

  ## Parse cluster variable
  if (!is.null(cluster)) {

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
      se_type <- "HC1"
    }
  }

  k <- ncol(X)

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", 1:k)
  }
  variable_names <- colnames(X)

  # Legacy, in case we want to only get some covs in the future
  which_covs <- rep(TRUE, ncol(X))

  if (!is.null(cluster)) {
    cl_ord <- order(cluster)
    y <- y[cl_ord]
    X <- X[cl_ord, , drop = FALSE]
    cluster <- cluster[cl_ord]
    J <- length(unique(cluster))
    if (!is.null(weights)) {
      weights <- weights[cl_ord]
    }
  } else {
    J <- 1
  }

  if (!is.null(weights)) {
    Xunweighted <- X
    weight_mean <- mean(weights)
    weights <- sqrt(weights / weight_mean)
    X <- weights * X
    y <- weights * y
  } else {
    weight_mean <- 1
    Xunweighted <- NULL
  }


  fit <-
    lm_solver(
      Xfull = X,
      y = y,
      Xunweighted = Xunweighted,
      weight = weights,
      weight_mean = weight_mean,
      cluster = cluster,
      J = J,
      ci = ci,
      type = se_type,
      which_covs = which_covs,
      try_cholesky = try_cholesky
    )


  return_frame <- data.frame(
    coefficients = setNames(as.vector(fit$beta_hat), variable_names),
    se = NA,
    df = NA,
    stringsAsFactors = FALSE
  )

  est_exists <- !is.na(return_frame$coefficients)

  N <- nrow(X)
  rank <- sum(est_exists)

  if (se_type != "none") {
    return_frame$se[est_exists] <- sqrt(diag(fit$Vcov_hat))

    if (ci) {
      if (se_type %in% cl_se_types) {

        # Replace -99 with NA, easy way to flag that we didn't compute
        # the DoF because the user didn't ask for it
        return_frame$df[est_exists] <-
          ifelse(fit$dof == -99,
            NA,
            fit$dof
          )
      } else {

        # TODO explicitly pass rank from RRQR/cholesky
        return_frame$df[est_exists] <- N - rank
      }
    }
  }

  # ----------
  # Build return object
  # ----------

  return_list <- add_cis_pvals(return_frame, alpha, ci && se_type != "none")

  return_list[["coefficient_name"]] <- variable_names
  return_list[["outcome"]] <- NA_character_
  return_list[["alpha"]] <- alpha
  return_list[["se_type"]] <- se_type

  return_list[["weighted"]] <- !is.null(weights)
  if (return_list[["weighted"]]) {
    return_list[["tot_var"]] <- ifelse(
      has_int,
      # everything correct except for this
      sum(weights ^ 2 * (y / weights - weighted.mean(y / weights, weights ^ 2)) ^ 2) * weight_mean,
      sum(y ^ 2 * weight_mean)
    )
    return_list[["res_var"]] <- sum(fit$residuals ^ 2 * weight_mean) / (N - rank)
  } else {
    return_list[["tot_var"]] <- ifelse(has_int, sum((y - mean(y)) ^ 2), sum(y ^ 2))
    return_list[["res_var"]] <- ifelse(fit$res_var < 0, NA, fit$res_var)
  }

  return_list[["df.residual"]] <- N - rank
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

  return_list[["fstatistic"]] <- c(
    value = (return_list[["r.squared"]] * return_list[["df.residual"]])
    / ((1 - return_list[["r.squared"]]) * (rank - has_int)),
    numdf = rank - has_int,
    dendf = return_list[["df.residual"]]
  )

  # return_list[["XtX_inv"]] <- fit$XtX_inv
  return_list[["N"]] <- N
  return_list[["k"]] <- k
  return_list[["rank"]] <- rank

  if (return_vcov && se_type != "none") {
    # return_list$residuals <- fit$residuals
    return_list[["vcov"]] <- fit$Vcov_hat
    dimnames(return_list[["vcov"]]) <- list(
      return_list$coefficient_name[est_exists],
      return_list$coefficient_name[est_exists]
    )
  }

  attr(return_list, "class") <- "lm_robust"

  return(return_list)
}
