#' Ordinary Least Squares with Robust Standard Errors
#'
#' Fits a linear model by ordinary least squares and returns
#' heteroskedasticity-robust or cluster-robust standard errors, with the
#' small-sample corrections used in design-based work. Fixed effects can be
#' absorbed rather than expanded into dummy columns, at no cost in the
#' available standard error types.
#'
#' @param formula (required) An object of class formula, as in [lm()]
#' @param data (optional) A `data.frame`
#' @param weights (optional) The bare (unquoted) name of the weights variable
#' @param subset (optional) A bare (unquoted) expression specifying a subset
#' @param clusters (optional) A bare (unquoted) name of the cluster variable
#' @param fixed_effects (optional) A one-sided formula of fixed effects to
#'   absorb rather than expand into dummy columns, such as `~ blockID` or
#'   `~ block + year`. Each variable is demeaned within the groups before OLS
#'   is run, so by the Frisch-Waugh-Lovell theorem the coefficients and
#'   residuals are the dummy regression's exactly.
#'
#'   Absorbing costs nothing in available standard error types. `"HC2"` and
#'   `"HC3"` are exact at any number of factors, because the leverage of the
#'   full design splits into the demeaned-X leverage plus a term that is cheap
#'   to compute, so no dummy hat matrix is built.
#'
#'   `"CR2"` is the exception: its adjustment is built from cluster-level
#'   blocks of the hat matrix rather than from the diagonal, and blocks do not
#'   split that way, so it expands the dummies and pays for the expansion.
#'   That is why `fixed_effects` with `clusters` defaults to `"CR0"`. Asking
#'   for `se_type = "CR2"` still works and still gives the 1.0.6 number.
#'   Refused is the three together: `"CR2"` with both `weights` and
#'   `fixed_effects`, as in estimatr 1.0.6.
#'
#'   The projection identity, the several-factor case, the exact-rank
#'   calculation, and the weighted CR2 and HC2 conventions are derived in
#'   `vignette("mathematical-notes")`.
#' @param se_type (optional) The standard error type. Defaults depend on whether clusters
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
#'   `"stata"` means two different things. With no clusters it is exactly
#'   `"HC1"`, and the fitted object reports `se_type = "HC1"`. With clusters it
#'   is **not** an alias for `"CR0"`: it is CR0 scaled by Stata's finite-sample
#'   factor, `(J / (J - 1)) * ((N - 1) / (N - K))` on the variance, and the
#'   object reports `se_type = "stata"` to keep the distinction visible.
#' @param ci (optional) Logical. Whether to compute p-values and confidence intervals.
#' @param alpha (optional) The significance level, 0.05 by default.
#' @param return_vcov (optional) Logical. Whether to return the vcov matrix.
#' @param try_cholesky (optional) Logical. Whether to solve by Cholesky
#'   decomposition of `X'X` rather than by the default pivoted QR. `FALSE` by
#'   default, and worth turning on in most applied settings: about 1.4 times
#'   faster at n = 100,000 with two regressors, and 1.7 times faster at
#'   n = 200,000 with 60 regressors, where it is 0.15s against 0.25s. The
#'   saving is per fit, so it is worth most in a simulation that fits the same
#'   design thousands of times.
#'
#'   Rank deficiency is caught on either path. Redundant columns come back as
#'   `NA` exactly as they do from [lm()] whichever path ran, and a design that
#'   is rank deficient falls back to the QR.
#'
#'   Whether it is safe turns on one question, whether two regressors are
#'   nearly the same variable. Forming `X'X` squares the condition number, so
#'   the Cholesky path has about twice the rounding error of the QR, and
#'   only near-collinearity makes that visible. Differences of scale do not,
#'   because the columns are normalized before either decomposition, so a
#'   covariate in dollars beside one in years costs nothing. For a treatment
#'   indicator, a few covariates, block or cluster dummies, the centered
#'   interactions [lm_lin()] builds, or a factorial, the two paths agree to at
#'   least 10 significant digits, which is why [difference_in_means()] sets it
#'   to `TRUE` internally. Agreement falls to about 3 digits as the scaled
#'   condition index reaches `1e6`, and the QR fallback takes over above
#'   roughly `1e8`. Nothing interpretable lives in that range: a design at
#'   `1e6` returns a coefficient of 4.8e4 with a standard error of 4.6e4 on a
#'   regressor whose true effect is zero. To check a design directly, scale the
#'   columns first, since the unscaled condition number of a design in mixed
#'   units is large for a reason that does not affect the fit:
#'   `kappa(sweep(X, 2, sqrt(colSums(X^2)), "/"), exact = TRUE)`.
#'
#' @return An object of class `"lm_robust"`, a list holding the estimate table in `coefficients`, `std.error`, `df`, `statistic`,
#'   `p.value`, `conf.low`, `conf.high`, `term`, and `outcome`; the fit in
#'   `fitted.values`, `residuals`, `vcov`, `nobs`, `k`, `rank`, `df.residual`,
#'   and `res_var`; the summary statistics `r.squared`, `adj.r.squared`,
#'   `tss`, and `fstatistic`; and `se_type`, `weighted`, `clustered`, `fes`,
#'   `alpha`, `terms`, `xlevels`, and `call`.
#'
#'   Absorbed fits add `fixed_effects`, `felevels` (the absorbed levels of
#'   each factor), and the within-projection summaries `proj_r.squared`,
#'   `proj_adj.r.squared`, `proj_tss`, and `proj_fstatistic`.
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
    # spanned by the others (a nested factor, or a disconnected design),
    # which inflates the rank correction and shrinks the residual degrees of
    # freedom. fe_leverage() returns the exact rank from the same
    # eigendecomposition that gives the leverage, so it is used for every
    # se_type; only HC2 and HC3 also need the vector.
    fe_proj <- fe_leverage(model_data[["fe_codes"]], model_data[["weights"]],
                           leverage = needs_fe_leverage(se_type, !is.null(model_data[["cluster"]])))
    fe_rank <- fe_proj[["rank"]]
    fe_lev <- fe_proj[["leverage"]]

    if (ncol(model_data$design_matrix) == 0L) {
      n_obs <- nrow(yoriginal)
      df_r <- n_obs - fe_rank
      residuals_proj <- drop(model_data$outcome)
      fitted_full <- drop(yoriginal) - residuals_proj
      fitted_full <- attach_obs_names(fitted_full, model_data)
      ss <- fe_r2(yoriginal, residuals_proj, model_data[["weights"]])
      rss <- ss[["rss"]]
      tss_full <- ss[["tss"]]
      r2_full <- 1 - rss / tss_full
      return_list <- list(
        # `term` must be present even though it is empty: `$term` on a list
        # partially matches `terms` when it is absent, and summarize_tidy()
        # then tried to use the terms formula as row names.
        term          = character(0),
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
        k             = 0L,
        proj_r.squared  = 0,
        proj_adj.r.squared = 0,
        proj_tss        = rss,
        r.squared       = r2_full,
        adj.r.squared   = 1 - (1 - r2_full) * (n_obs - 1L) / df_r,
        alpha           = alpha,
        clustered       = !is.null(model_data[["cluster"]]),
        contrasts     = attr(model_data$design_matrix, "contrasts"),
        terms         = model_data$terms,
        xlevels       = model_data$xlevels,
        felevels      = model_data$fe_level_names,
        tss           = tss_full,
        weights       = model_data$weights,
        outcome       = deparse(formula[[2]], nlines = 5)
      )
      # Absorbed group effects, as the ordinary fixed-effects path stores them,
      # so predict() can put them back. With no regressors they are the fitted
      # values themselves.
      return_list[["fixed_effects"]] <- absorbed_group_effects(
        return_list[["fitted.values"]], return_list[["coefficients"]], model_data
      )
      return_list[["call"]] <- match.call()
      # Without this the fit came back as a bare list, so print(), tidy(),
      # summary() and every other method dispatched on nothing.
      class(return_list) <- "lm_robust"
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
      femat = if (has_fe && needs_fe_dummies(se_type))
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
    fitted_full <- attach_obs_names(fitted_full, model_data)
    return_list[["fitted.values"]] <- fitted_full

    # Absorbed group effects, so predict() can put them back.
    return_list[["fixed_effects"]] <- absorbed_group_effects(
      return_list[["fitted.values"]], return_list[["coefficients"]], model_data
    )

    # Full model R2 using original Y, weighted where the fit is and one value
    # per outcome column, as the same model with explicit dummies reports.
    n_obs <- nrow(yoriginal)
    ss <- fe_r2(yoriginal, residuals_proj, model_data[["weights"]])
    tss_full <- ss[["tss"]]
    rss_full <- ss[["rss"]]
    r2_full  <- 1 - rss_full / tss_full
    return_list[["r.squared"]]     <- r2_full
    return_list[["adj.r.squared"]] <- 1 - (1 - r2_full) * (n_obs - 1L) / return_list[["df.residual"]]
    return_list[["tss"]]           <- tss_full
    return_list[["felevels"]]      <- model_data[["fe_level_names"]]
  }

  return_list[["call"]] <- match.call()

  return(return_list)
}
