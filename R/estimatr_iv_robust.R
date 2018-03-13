#' Two-Stage Least Squares Instrumental Variables Regression
#'
#' @description This formula estimates an instrumental variables regression
#' using two-stage least squares with a variety of
#' options for robust standard errors
#'
#' @param formula an object of class formula of the regression and the instruments.
#' For example, the formula \code{y ~ x1 + x2 | z1 + z2} specifies \code{x1} and \code{x2}
#' as endogenous regressors and \code{z1} and \code{z2} as their respective instruments.
#' @param data A \code{data.frame}
#' @param weights the bare (unquoted) names of the weights variable in the
#' supplied data.
#' @param subset An optional bare (unquoted) expression specifying a subset
#' of observations to be used.
#' @param clusters An optional bare (unquoted) name of the variable that
#' corresponds to the clusters in the data.
#' @param se_type The sort of standard error sought. If `clusters` is
#' not specified the options are "HC0", "HC1" (or "stata", the equivalent),
#'  "HC2" (default), "HC3", or
#' "classical". If `clusters` is specified the options are "CR0", "CR2" (default), or "stata".
#' @param ci logical. Whether to compute and return p-values and confidence
#' intervals, TRUE by default.
#' @param alpha The significance level, 0.05 by default.
#' @param return_vcov logical. Whether to return the variance-covariance
#' matrix for later usage, TRUE by default.
#' @param try_cholesky logical. Whether to try using a Cholesky
#' decomposition to solve least squares instead of a QR decomposition,
#' FALSE by default. Using a Cholesky decomposition may result in speed gains, but should only
#' be used if users are sure their model is full-rank (i.e., there is no
#' perfect multi-collinearity)
#'
#' @details
#'
#' This function performs two-stage least squares estimation of
#'
#' @return An object of class \code{"iv_robust"}.
#'
#' For now, only use \code{tidy} on this object.
#'
#' @examples
#' library(fabricatr)
#' dat <- fabricate(
#'   N = 40,
#'   y = rpois(N, lambda = 4),
#'   Z = rbinom(N, 1, prob = 0.4),
#'   D  = Z * rbinom(N, 1, prob = 0.8),
#'   X = rnorm(N)
#' )
#'
#' # Instrument for treatment `D` with encouragement `Z`
#' tidy(iv_robust(y ~ D + X | Z + X, data = dat))
#'
#' @export
iv_robust <- function(formula,
                      data,
                      weights,
                      subset,
                      clusters,
                      se_type = NULL,
                      ci = TRUE,
                      alpha = .05,
                      return_vcov = TRUE,
                      try_cholesky = FALSE) {

  datargs <- enquos(
    formula = formula,
    weights = weights,
    subset = subset,
    cluster = clusters
  )
  data <- enquo(data)
  model_data <- clean_model_data(data = data, datargs)

  # -----------
  # First stage
  # -----------

  first_stage <-
    lm_robust_fit(
      y = model_data$design_matrix,
      X = model_data$instrument_matrix,
      weights = model_data$weights,
      cluster = model_data$cluster,
      ci = FALSE,
      se_type = "none",
      has_int = attr(model_data$terms, "intercept"),
      alpha = alpha,
      return_fit = TRUE,
      return_unweighted_fit = TRUE,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky
    )

  # ------
  # Second stage
  # ------
  colnames(first_stage$fitted.values) <- colnames(model_data$design_matrix)

  second_stage <-
    lm_robust_fit(
      y = model_data$outcome,
      X = first_stage$fitted.values,
      weights = model_data$weights,
      cluster = model_data$cluster,
      ci = ci,
      se_type = se_type,
      has_int = attr(model_data$terms, "intercept"),
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      X_first_stage = model_data$design_matrix
    )

  return_list <- lm_return(
    second_stage,
    model_data = model_data,
    formula = formula
  )

  return_list[["call"]] <- match.call()

  return_list[["terms_regressors"]] <- model_data[["terms_regressors"]]
  class(return_list) <- "iv_robust"

  return(return_list)
}
