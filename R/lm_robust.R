#' Ordinary Least Squares with Robust Standard Errors
#'
#' @param formula an object of class formula, as in [lm()]
#' @param data A `data.frame`
#' @param weights the bare (unquoted) name of the weights variable
#' @param subset An optional bare (unquoted) expression specifying a subset
#' @param clusters An optional bare (unquoted) name of the cluster variable
#' @param fixed_effects An optional one-sided formula of fixed effects to absorb,
#'   such as `~ blockID` or `~ block + year`. Uses the Frisch-Waugh-Lovell (FWL)
#'   theorem: each variable is demeaned within FE groups before OLS is run. FWL
#'   guarantees exact coefficient and residual recovery.
#'
#'   There is no SE restriction, at any number of FE factors: `"HC2"` and
#'   `"HC3"` are exact and cost no more than the others. The projection onto
#'   the full \[dummies | X\] design splits as `P_[X|D] = P_D + P_{M_D X}`, so
#'   `h_ii` is the demeaned-X hat value plus `diag(P_D)`, and no dummy hat
#'   matrix is built. With one factor `diag(P_D)` is `w_i / sum(w in group i)`;
#'   with several it costs a factorisation the size of the design's narrowest
#'   dimension. Results are identical to writing the dummies out, and to
#'   estimatr, at a fraction of the time and memory.
#'
#'   `"CR2"` is the one exception, and needs the dummies whatever the number of
#'   factors, since its adjustment is built from cluster-level blocks of the hat
#'   matrix rather than from `h_ii`. It is available with `fixed_effects` but
#'   pays for the expansion, and is refused in combination with `weights`, as in
#'   estimatr 1.0.6.
#' @param se_type The standard error type. Defaults depend on whether clusters
#'   and/or fixed effects are present:
#'   \itemize{
#'     \item No clusters, no FE: `"HC2"` (default), `"HC0"`, `"HC1"`,
#'       `"HC3"`, `"classical"`, `"stata"`, `"none"`.
#'     \item Clusters, no FE: `"CR2"` (default), `"CR0"`, `"stata"`, `"none"`.
#'     \item No clusters, with FE (any number of factors): `"HC2"` (default),
#'       `"HC0"`, `"HC1"`, `"HC3"`, `"classical"`, `"stata"`, `"none"`. The
#'       same menu as with no FE at all.
#'     \item Clusters, with FE: `"CR0"` (default), `"CR2"`, `"stata"`,
#'       `"none"`. `"CR2"` expands the fixed effects into dummies, so it is not
#'       the default here; it is refused with `weights`.
#'   }
#'   `"stata"` is an alias for HC1 (no clusters) or CR0 (with clusters).
#' @param ci logical. Whether to compute p-values and confidence intervals.
#' @param alpha The significance level, 0.05 by default.
#' @param return_vcov logical. Whether to return the vcov matrix.
#' @param try_cholesky logical. Whether to try Cholesky decomposition.
#'
#' @return An object of class `"lm_robust"`.
#'
#' @examples
#' set.seed(15)
#' dat <- data.frame(
#'   y = rpois(40, lambda = 4),
#'   x = rnorm(40),
#'   z = rbinom(40, 1, prob = 0.4),
#'   cl = rep(1:10, each = 4),
#'   bl = rep(c("A", "B", "C", "D"), each = 10),
#'   w = runif(40)
#' )
#'
#' # HC2 is the default
#' fit <- lm_robust(y ~ x + z, data = dat)
#' fit
#' tidy(fit)
#' summary(fit)
#' confint(fit, level = 0.8)
#'
#' # Other variance estimators, including Stata's
#' lm_robust(y ~ x + z, data = dat, se_type = "classical")
#' lm_robust(y ~ x + z, data = dat, se_type = "stata")
#'
#' # Clustered inference defaults to CR2
#' lm_robust(y ~ x + z, data = dat, clusters = cl)
#' lm_robust(y ~ x + z, data = dat, clusters = cl, se_type = "stata")
#'
#' # Weights and subsets behave as they do in lm()
#' lm_robust(y ~ x + z, data = dat, weights = w, clusters = cl)
#' lm_robust(y ~ x, data = dat, subset = z == 1)
#'
#' # Fixed effects are absorbed rather than expanded into dummies. With a
#' # single factor the HC2 default is exact and costs nothing extra.
#' lm_robust(y ~ z, data = dat, fixed_effects = ~ bl)
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

  has_fe  <- !is.null(model_data[["fixed_effects"]])
  fe_rank <- 0L
  fe_lev <- NULL
  yoriginal <- NULL

  if (has_fe) {
    yoriginal  <- as.matrix(model_data[["outcome"]])
    model_data <- demean_fes(model_data)
    # fe_rank: degrees of freedom consumed by FE (levels - 1 per variable, +1 for absorbed intercept)
    fe_rank <- sum(model_data[["fe_levels"]]) - length(model_data[["fe_levels"]]) + 1L

    # The nominal count above overstates the rank whenever one factor is partly
    # spanned by the others -- a nested factor, or a disconnected design --
    # which inflates the rank correction and shrinks the residual degrees of
    # freedom. fe_leverage() returns the exact rank from the same
    # eigendecomposition that gives the leverage, so it is used for every
    # se_type; only HC2 and HC3 also need the vector.
    fe_proj <- fe_leverage(model_data[["fe_codes"]], model_data[["weights"]],
                           leverage = needs_fe_leverage(se_type))
    fe_rank <- fe_proj[["rank"]]
    fe_lev <- fe_proj[["leverage"]]

    if (ncol(model_data$design_matrix) == 0L) {
      n_obs <- nrow(yoriginal)
      df_r <- n_obs - fe_rank
      residuals_proj <- drop(model_data$outcome)
      fitted_full <- drop(yoriginal) - residuals_proj
      if (!is.null(model_data[["obs_names"]])) {
        if (is.matrix(fitted_full)) {
          rownames(fitted_full) <- model_data[["obs_names"]]
        } else {
          names(fitted_full) <- model_data[["obs_names"]]
        }
      }
      rss <- sum(diag(crossprod(residuals_proj)))
      tss_full <- sum(diag(crossprod(yoriginal - mean(yoriginal))))
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
        fitted.values = fitted_full,
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
        felevels      = model_data$fe_level_names,
        tss           = tss_full,
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
      fe_rank = fe_rank,
      fe_leverage = fe_lev,
      femat = if (has_fe && needs_fe_dummies(se_type, model_data))
        fe_dummy_matrix(model_data)
        else NULL
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
    residuals_proj <- drop(return_list[["residuals"]])
    fitted_full <- drop(yoriginal) - residuals_proj
    # yoriginal comes from the stripped model data, so the names go back on
    # here, exactly as lm_return() does it for a fit without fixed effects.
    if (!is.null(model_data[["obs_names"]])) {
      if (is.matrix(fitted_full)) {
        rownames(fitted_full) <- model_data[["obs_names"]]
      } else {
        names(fitted_full) <- model_data[["obs_names"]]
      }
    }
    return_list[["fitted.values"]] <- fitted_full

    # Absorbed group effects, so predict() can put them back.
    return_list[["fixed_effects"]] <- absorbed_group_effects(
      return_list[["fitted.values"]], return_list[["coefficients"]], model_data
    )

    # Full model R2 using original Y
    n_obs <- nrow(yoriginal)
    y_mean <- mean(yoriginal)
    # `crossprod()` reads the columns without materialising their squares, so
    # the two sums of squares cost one n-length temporary between them instead
    # of three. A multivariate outcome is n-by-k here, and both totals are
    # pooled across the k columns as they have always been, so it is the trace
    # that is wanted and not the whole k-by-k matrix.
    tss_full <- sum(diag(crossprod(yoriginal - y_mean)))
    rss_full <- sum(diag(crossprod(residuals_proj)))
    r2_full  <- 1 - rss_full / tss_full
    return_list[["r.squared"]]     <- r2_full
    return_list[["adj.r.squared"]] <- 1 - (1 - r2_full) * (n_obs - 1L) / return_list[["df.residual"]]
    return_list[["tss"]]           <- tss_full
    return_list[["felevels"]]      <- model_data[["fe_level_names"]]
  }

  return_list[["call"]] <- match.call()

  return(return_list)
}
