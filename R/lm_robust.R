#' Ordinary Least Squares with Robust Standard Errors
#'
#' @param formula an object of class formula, as in [lm()]
#' @param data A `data.frame`
#' @param weights the bare (unquoted) name of the weights variable
#' @param subset An optional bare (unquoted) expression specifying a subset
#' @param clusters An optional bare (unquoted) name of the cluster variable
#' @param fixed_effects An optional right-sided formula of fixed effects to
#'   absorb, such as `~ blockID`. Uses Frisch-Waugh-Lovell demeaning (feols
#'   style). Defaults to HC1 (no clusters) or stata (with clusters). HC2, HC3,
#'   and CR2 are not available with fixed effects and will warn and fall back.
#' @param se_type The sort of standard error. Without clusters: "HC0", "HC1",
#'   "HC2" (default, or HC1 with FE), "HC3", "classical", "stata", or "none".
#'   With clusters: "CR0", "CR2" (default, or stata with FE), "stata", or "none".
#' @param ci logical. Whether to compute p-values and confidence intervals.
#' @param alpha The significance level, 0.05 by default.
#' @param return_vcov logical. Whether to return the vcov matrix.
#' @param try_cholesky logical. Whether to try Cholesky decomposition.
#'
#' @return An object of class `"lm_robust"`.
#'
#' @export
lm_robust <- function(formula,
                      data,
                      weights,
                      subset,
                      clusters,
                      fixed_effects,
                      se_type = NULL,
                      ci = TRUE,
                      alpha = .05,
                      return_vcov = TRUE,
                      try_cholesky = FALSE) {
  datargs <- rlang::enquos(
    formula = formula,
    weights = weights,
    subset = subset,
    cluster = clusters,
    fixed_effects = fixed_effects
  )
  data <- rlang::enquo(data)
  model_data <- clean_model_data(data = data, datargs)

  has_fe  <- is.character(model_data[["fixed_effects"]])
  fe_rank <- 0L
  yoriginal <- NULL

  if (has_fe) {
    yoriginal  <- as.matrix(model_data[["outcome"]])
    model_data <- demean_fes(model_data)
    # fe_rank: degrees of freedom consumed by FE (levels - 1 per variable, +1 for absorbed intercept)
    fe_rank <- sum(model_data[["fe_levels"]]) - length(model_data[["fe_levels"]]) + 1L

    if (ncol(model_data$design_matrix) == 0L) {
      n_obs <- nrow(yoriginal)
      df_r <- n_obs - fe_rank
      residuals_proj <- drop(model_data$outcome)
      rss <- sum(residuals_proj^2)
      tss_full <- sum((yoriginal - mean(yoriginal))^2)
      r2_full <- 1 - rss / tss_full
      return_list <- list(
        coefficients  = setNames(numeric(0), character(0)),
        std.error     = setNames(numeric(0), character(0)),
        statistic     = setNames(numeric(0), character(0)),
        p.value       = setNames(numeric(0), character(0)),
        conf.low      = setNames(numeric(0), character(0)),
        conf.high     = setNames(numeric(0), character(0)),
        df            = setNames(numeric(0), character(0)),
        df.residual   = df_r,
        res_var       = rss / df_r,
        vcov          = matrix(numeric(0), 0L, 0L),
        fitted.values = drop(yoriginal) - residuals_proj,
        residuals     = residuals_proj,
        weighted      = !is.null(model_data$weights),
        se_type       = "none",
        fes           = TRUE,
        nobs          = n_obs,
        rank          = 0L,
        proj_r.squared  = 0,
        r.squared       = r2_full,
        adj.r.squared   = 1 - (1 - r2_full) * (n_obs - 1L) / df_r,
        contrasts     = attr(model_data$design_matrix, "contrasts"),
        terms         = model_data$terms,
        xlevels       = model_data$xlevels,
        felevels      = model_data$fe_levels,
        weights       = model_data$weights,
        outcome       = deparse(formula[[2]], nlines = 5)
      )
      return_list[["call"]] <- match.call()
      return(return_list)
    }
  }

  return_list <-
    lm_robust_fit(
      y = model_data$outcome,
      X = model_data$design_matrix,
      weights = model_data$weights,
      cluster = model_data$cluster,
      ci = ci,
      se_type = se_type,
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      has_int = attr(model_data$terms, "intercept"),
      iv_stage = list(0),
      fe_rank = fe_rank
    )

  return_list <- lm_return(
    return_list,
    model_data = model_data,
    formula = formula
  )

  if (has_fe) {
    # Rename projected (demeaned) R2 stats
    for (nm in c("r.squared", "adj.r.squared", "tss", "fstatistic")) {
      if (!is.null(return_list[[nm]])) {
        return_list[[paste0("proj_", nm)]] <- return_list[[nm]]
        return_list[[nm]] <- NULL
      }
    }

    # Reconstruct full fitted values: by FWL, projected residuals = full residuals
    residuals_proj <- drop(model_data[["outcome"]]) - return_list[["fitted.values"]]
    return_list[["fitted.values"]] <- drop(yoriginal) - residuals_proj

    # Full model R2 using original Y
    n_obs <- nrow(yoriginal)
    y_mean <- mean(yoriginal)
    tss_full <- sum((yoriginal - y_mean)^2)
    rss_full <- sum(residuals_proj^2)
    r2_full  <- 1 - rss_full / tss_full
    return_list[["r.squared"]]     <- r2_full
    return_list[["adj.r.squared"]] <- 1 - (1 - r2_full) * (n_obs - 1L) / return_list[["df.residual"]]
  }

  return_list[["call"]] <- match.call()

  return(return_list)
}
