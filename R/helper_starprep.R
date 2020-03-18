#' Build lm_robust object from lm fit
#'
#' @param model an lm model object
#' @param se_type The sort of standard error sought. If \code{clusters} is
#' not specified the options are "HC0", "HC1" (or "stata", the equivalent),
#' "HC2" (default), "HC3", or "classical". If \code{clusters} is specified the
#' options are "CR0", "CR2" (default), or "stata". Can also specify "none",
#' which may speed up estimation of the coefficients.
#' @param clusters A vector corresponding to the clusters in the data.
#' @param ci logical. Whether to compute and return p-values and confidence
#' intervals, TRUE by default.
#' @param alpha The significance level, 0.05 by default.
#'
#' @return an \code{\link{lm_robust}} object.
#'
#' @examples
#' lmo <- lm(mpg ~ hp, data = mtcars)
#'
#' # Default HC2
#' commarobust(lmo)
#'
#' commarobust(lmo, se_type = "HC3")
#'
#' commarobust(lmo, se_type = "stata", clusters = mtcars$carb)
#'
#' @export
commarobust <- function(model,
                        se_type = NULL,
                        clusters = NULL,
                        ci = TRUE,
                        alpha = 0.05) {

  if (class(model)[1] != "lm") {
    stop("`model` must be an lm object")
  }

  coefs <- as.matrix(coef(model))
  est_exists <- !is.na(coefs)
  covs_used <- which(est_exists)
  coefs <- coefs[covs_used, , drop = FALSE]

  Qr <- qr(model)
  p1 <- seq_len(model$rank)

  XtX_inv <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  clustered <- !is.null(clusters)

  se_type <- check_se_type(se_type = se_type, clustered = clustered)

  X <- model.matrix.default(model)
  contrasts <- attr(X, "contrasts")
  N <- nrow(X)

  data <- list(
    y = as.matrix(model.response(model$model)),
    X = X,
    weights = weights(model)
  )

  weighted <- is.numeric(data[["weights"]])

  ## check clusters and add to data
  if (clustered) {
    if (is.matrix(clusters) && ncol(clusters) > 1) {
      stop("`clusters` must be a single vector or column denoting the clusters.")
    }
    if (length(clusters) != N) {
      stop("`clusters` must be the same length as the model data.")
    }
    data[["cluster"]] <- as.factor(clusters)
  }

  data <- prep_data(
    data = data,
    se_type = se_type,
    clustered = clustered,
    weighted = weighted,
    fes = FALSE,
    iv_stage = list(0)
  )

  ei <- as.matrix(resid(model))
  if (clustered) {
    ei <- ei[data[["cl_ord"]], , drop = FALSE]
  }

  if (any(!est_exists)) {
    data <- drop_collinear(data, covs_used, weighted, FALSE)
  }

  if (weighted) {
    ei <- data[["weights"]] * ei
    XtX_inv <- data[["weight_mean"]] * XtX_inv

    # Need unweighted resid and need to reweight X
    if (se_type == "CR2") {
      eiunweighted <-
        as.matrix(data[["yunweighted"]] - data[["Xunweighted"]] %*% coefs)
      data[["X"]] <- data[["weights"]] * data[["X"]]
    }
  }

  vcov_fit <- lm_variance(
    X = data[["X"]],
    Xunweighted = data[["Xunweighted"]],
    XtX_inv = XtX_inv,
    ei = if (se_type == "CR2" && weighted) eiunweighted else ei,
    weight_mean = data[["weight_mean"]],
    cluster = data[["cluster"]],
    J = data[["J"]],
    ci = ci,
    se_type = se_type,
    which_covs = rep(TRUE, model$rank),
    fe_rank = 0
  )


  ## Build return_list
  return_list <- list(
    coefficients = as.matrix(coef(model)),
    std.error = NA,
    df = NA,
    term = names(coef(model)),
    outcome = as.character(rlang::f_lhs(formula(model))),
    alpha = alpha,
    se_type = se_type,
    df.residual = df.residual(model),
    weighted = weighted,
    fes = FALSE,
    clustered = clustered,
    N = nobs(model),
    rank = model$rank,
    k = ncol(X),
    fitted.values = fitted.values(model),
    contrasts = contrasts,
    terms = model$terms,
    xlevels = model$xlevels,
    weights = weights(model)
  )
  return_list[["std.error"]][est_exists] <- sqrt(diag(vcov_fit$Vcov_hat))
  return_list[["df"]][est_exists] <- ifelse(vcov_fit$dof == -99, NA, vcov_fit$dof)
  if (clustered) {
    return_list[["N_clusters"]] <- data[["J"]]
  }

  return_list[["res_var"]] <- get_resvar(
    data = data,
    ei = ei,
    df.residual = return_list[["df.residual"]],
    vcov_fit = vcov_fit,
    weighted = weighted
  )

  return_list <- add_cis_pvals(return_list, alpha, ci && se_type != "none")

  ## Add F stat
  tss_r2s <- get_r2s(
    y = data[["y"]],
    return_list = return_list,
    yunweighted = data[["yunweighted"]],
    has_int = attr(model$terms, "intercept"),
    weights = data[["weights"]],
    weight_mean = data[["weight_mean"]]
  )

  if (clustered) {
    dendf <- data[["J"]] - 1
  } else {
    dendf <- return_list[["df.residual"]]
  }

  f <- get_fstat(
    tss_r2s = tss_r2s,
    return_list = return_list,
    iv_ei = NULL,
    nomdf = model$rank - attr(model$terms, "intercept"),
    dendf = dendf,
    vcov_fit = vcov_fit,
    has_int = attr(model$terms, "intercept"),
    iv_stage = list(0)
  )

  return_list <- c(return_list, tss_r2s)
  return_list[["fstatistic"]] <- f

  return_list[["vcov"]] <- vcov_fit$Vcov_hat
  dimnames(return_list[["vcov"]]) <- list(
    return_list$term[est_exists],
    return_list$term[est_exists]
  )

  return_list <- lm_return(return_list, model_data = NULL, formula = NULL)

  attr(return_list, "class") <- "lm_robust"

  return(return_list)
}


#' Prepare model fits for stargazer
#'
#' @param ... a list of lm_robust or lm objects
#' @param stat either "std.error" (the default), "statistic" (the t-statistic), "p.value", "ci", or "df"
#' @param se_type (optional) if any of the objects are lm objects, what standard
#' errors should be used. Must only be one type and will be used for all lm
#' objects passed to starprep. See \code{commarobust} for more.
#' @param clusters (optional) if any of the objects are lm objects, what clusters
#' should be used, if clusters should be used. Must only be one vector and will
#' be used for all lm objects passed to starprep. See \code{commarobust} for more.
#' @param alpha (optional) if any of the objects are lm objects, what significance level
#' should be used for the p-values or confidence intervals
#'
#' @details Used to help extract statistics from lists of model fits for stargazer.
#' Prefers lm_robust objects, but because \code{stargazer} does not work with \code{lm_robust}
#' objects, \code{starprep} can also take \code{lm} objects and calls \code{commarobust} to get
#' the preferred, robust statistics.
#'
#' @return a list of vectors of extracted statistics for stargazers
#'
#' @examples
#'
#' library(stargazer)
#'
#' lm1 <- lm(mpg ~ hp, data = mtcars)
#' lm2 <- lm(mpg ~ hp + wt, data = mtcars)
#'
#' # Use default "HC2" standard errors
#' stargazer(lm1, lm2,
#'           se = starprep(lm1, lm2),
#'           p = starprep(lm1, lm2, stat = "p.value"),
#'           omit.stat = "f")
#' # NB: We remove the F-stat because stargazer only can use original F-stat
#' # which uses classical SEs
#'
#' # Use default "CR2" standard errors with clusters
#' stargazer(lm1, lm2,
#'           se = starprep(lm1, lm2, clusters = mtcars$carb),
#'           p = starprep(lm1, lm2, clusters = mtcars$carb, stat = "p.value"),
#'           omit.stat = "f")
#'
#' # Can also specify significance levels and different standard errors
#' stargazer(lm1, lm2,
#'           ci.custom = starprep(lm1, lm2, se_type = "HC3", alpha = 0.1, stat = "ci"),
#'           omit.stat = "f")
#'
#' @export
starprep <- function(...,
                     stat = c("std.error", "statistic", "p.value", "ci", "df"),
                     se_type = NULL,
                     clusters = NULL,
                     alpha = 0.05) {

  if (inherits(..1, "list")) {
    if (...length() > 1) {
      stop("`...` must be one list of model fits or several comma separated model fits")
    }
    fits <- ..1
  } else {
    fits <- list(...)
  }

  is_list_of_lm <- vapply(fits, inherits, what = c("lm","lm_robust"), TRUE)

  if (any(!is_list_of_lm)) {
    stop("`...` must contain only `lm` or `lm_robust` objects.")
  }

  fitlist <- lapply(
    fits,
    function(x) {
      if (inherits(x, "lm"))
        commarobust(x, se_type = se_type, clusters = clusters, alpha = alpha)
      else
        x
    }
  )

  stat <- match.arg(stat)

  if (stat == "ci") {
    out <- lapply(fitlist, function(x) cbind(x[["conf.low"]], x[["conf.high"]]))
  } else {
    out <- lapply(fitlist, `[[`, stat)
  }

  return(out)
}
