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
                          has_int, # TODO get this out of here
                          alpha = 0.05,
                          return_vcov = FALSE,
                          return_fit = TRUE,
                          try_cholesky = FALSE,
                          X_first_stage = NULL) {


  weighted <- !is.null(weights)
  iv <- !is.null(X_first_stage)
  clustered <- !is.null(cluster)

  # ----------
  # Check se type
  # ----------

  # TODO what is implemented for IV?

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
    y <- as.matrix(y)[cl_ord, , drop = FALSE]
    X <- X[cl_ord, , drop = FALSE]
    cluster <- cluster[cl_ord]
    J <- length(unique(cluster))
    if (weighted) {
      weights <- weights[cl_ord]
    }
    if (iv) {
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
    if (iv) {
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

  return_frame <- data.frame(
    coefficients = setNames(as.vector(fit$beta_hat), variable_names),
    se = NA,
    df = NA,
    stringsAsFactors = FALSE
  )

  est_exists <- !is.na(return_frame$coefficients)
  N <- nrow(X)
  rank <- sum(est_exists)

  # ----------
  # Estimate variance
  # ----------

  if (se_type != "none" || return_fit) {

    if (rank < ncol(X)) {
      X <- X[, est_exists, drop = FALSE]
      if (weighted){
        Xunweighted <- Xunweighted[, est_exists, drop = FALSE]
      }
      if (iv) {
        X_first_stage <- X_first_stage[, est_exists, drop = FALSE]
        if (weighted) {
          X_first_stage_unweighted <- X_first_stage_unweighted[, est_exists, drop = FALSE]
        }
      }
      fit$beta_hat <- fit$beta_hat[est_exists]
    }

    # compute fitted.values and residuals
    # need unweighted for CR2, as well as X weighted by weights again
    # so that instead of having X * sqrt(W) we have X * W
    if (se_type == "CR2" && weighted) {
      if (iv) {
        fitted.values <- X_first_stage_unweighted %*% fit$beta_hat
      } else {
        fitted.values <- Xunweighted %*% fit$beta_hat
      }
      ei <- yunweighted - fitted.values
      X <- weights * X
    } else {
      if (iv) {
        fitted.values <- X_first_stage %*% fit$beta_hat
      } else {
        fitted.values <- X %*% fit$beta_hat
      }
      ei <- y - fitted.values
    }

    # TODO deal with multiple outcomes
    if (se_type == "CR2") {
      vcov_fit <- lm_variance_cr2(
        X = X,
        Xunweighted = Xunweighted,
        XtX_inv = fit$XtX_inv,
        beta_hat = fit$beta_hat,
        ei = ei,
        weight_mean = weight_mean,
        clusters = cluster,
        J = J,
        ci = ci,
        which_covs = which_covs[est_exists]
      )

      vcov_fit[["res_var"]] <-
        sum((y - X %*% fit$beta_hat)^2) /
        (N - rank)

    } else if (se_type != "none") {
      vcov_fit <- lm_variance(
        X = X,
        XtX_inv = fit$XtX_inv,
        beta_hat = fit$beta_hat,
        ei = ei,
        cluster = cluster,
        J = J,
        ci = ci,
        type = se_type,
        which_covs = which_covs[est_exists]
      )
    }

    if (se_type != "none") {
      return_frame$se[est_exists] <- sqrt(diag(vcov_fit$Vcov_hat))
    }

    if (ci && se_type != "none") {
      if (se_type %in% cl_se_types) {

        # Replace -99 with NA, easy way to flag that we didn't compute
        # the DoF because the user didn't ask for it
        return_frame$df[est_exists] <-
          ifelse(vcov_fit$dof == -99,
                 NA,
                 vcov_fit$dof
          )
      } else {

        # TODO explicitly pass rank from RRQR/cholesky
        return_frame$df[est_exists] <- N - rank
      }
    }
  }

  # ----------
  # Augment return object
  # ----------

  return_list <- add_cis_pvals(return_frame, alpha, ci && se_type != "none")

  if (return_fit) {
    if ((se_type == "CR2" && weighted) || iv) {
      # Have to get weighted fits as original fits were unweighted for
      # variance estimation or used wrong matrix for iv
      return_list[["fitted.values"]] <- y - X %*% fit$beta_hat
    } else {
      return_list[["fitted.values"]] <- fitted.values
    }

    # If we reordered to get SEs earlier, have to fix order
    if (clustered && se_type != "none") {
      return_list[["fitted.values"]] <- return_list[["fitted.values"]][order(cl_ord)]
    }
  }

  return_list[["coefficient_name"]] <- variable_names
  return_list[["outcome"]] <- NA_character_
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
      return_list[["tot_var"]] <- ifelse(
        has_int,
        sum(weights ^ 2 * (yunweighted - weighted.mean(yunweighted, weights ^ 2)) ^ 2) * weight_mean,
        sum(y ^ 2 * weight_mean)
      )
      return_list[["res_var"]] <- sum(ei ^ 2 * weight_mean) / (N - rank)
    } else {
      return_list[["tot_var"]] <- ifelse(has_int, sum((y - mean(y)) ^ 2), sum(y ^ 2))
      return_list[["res_var"]] <- ifelse(vcov_fit$res_var < 0, NA, vcov_fit$res_var)
    }

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

    if (return_vcov) {
      # return_list$residuals <- fit$residuals
      return_list[["vcov"]] <- vcov_fit$Vcov_hat
      dimnames(return_list[["vcov"]]) <- list(
        return_list$coefficient_name[est_exists],
        return_list$coefficient_name[est_exists]
      )
    }
  }

  attr(return_list, "class") <- "lm_robust"

  return(return_list)
}
