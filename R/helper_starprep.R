#' Build lm_robust object from lm fit
#'
#' @export
commarobust <- function(model, se_type, clusters = NULL, alpha = 0.05, ci = TRUE) {

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

  se_type <- check_se_type(
    se_type = se_type,
    clustered = clustered,
    iv = FALSE
  )

  X <- model.matrix.lm(model)
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
    data[["cluster"]] <- clusters
  }

  data <- prep_data(
    data = data,
    se_type = se_type,
    clustered = clustered,
    weighted = weighted,
    fes = FALSE,
    iv_second_stage = FALSE
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

    # Need unweighted resid and need
    if (se_type == "CR2") {
      eiunweighted <- as.matrix(data[["yunweighted"]] - data[["Xunweighted"]] %*% coefs)
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
    k = model$rank,
    fitted.values = fitted.values(model),
    contrasts = contrasts,
    terms = model$terms,
    xlevels = model$xlevels,
    weights = weights(model)
  )
  return_list[["std.error"]][est_exists] <- sqrt(diag(vcov_fit$Vcov_hat))
  return_list[["df"]][est_exists] <- ifelse(vcov_fit$dof == -99, NA, vcov_fit$dof)

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
    iv_second_stage = FALSE
  )

  return_list <- c(return_list, tss_r2s)
  return_list[["fstatistic"]] <- f

  return_list[["vcov"]] <- vcov_fit$Vcov_hat
  dimnames(return_list[["vcov"]]) <- list(
    return_list$term[est_exists],
    return_list$term[est_exists]
  )

  attr(return_list, "class") <- "lm_robust"

  return(return_list)
}

#' Prepare model fits for stargazer
#'
#' @param ... a list of lm_robust objects
#' @param stat either "se" (the default) or "p"
#'
#' @return a list of vectors of extracted statistics for stargazers
#'
#' @examples
#'
#' @export
starprep <- function(..., stat = "se"){
  fitlist = list(...)
  if(stat == "se") {
    out <- lapply(fitlist, function(x) x$std.error)
  } else if (stat == "p") {
    out <- lapply(fitlist, function(x){
      ps <- x$p.value
      ps["(Intercept)"] <- 1
      ps
    })
  }
  return(out)
}
