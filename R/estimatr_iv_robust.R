#' Two-Stage Least Squares Instrumental Variables Regression
#'
#' @description This formula estimates an instrumental variables regression
#' using two-stage least squares with a variety of options for robust
#' standard errors.
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
#' This function performs two-stage least squares estimation to fit
#' instrumental variables regression. The syntax is similar to that in
#' \code{ivreg} from the \code{AER} package. Regressors and instruments
#' should be specified in a two-part formula, such as
#' \code{y ~ x1 + x2 | z1 + z2 + z3}, where \code{x1} and \code{x2} are
#' regressors and \code{z1}, \code{z2}, and \code{z3} are instruments.
#'
#' Default variance estimators
#'
#' Math notes for further information
#'
#' @return An object of class \code{"iv_robust"}.
#'
#' Methods
#' object names
#'
#' @examples
#' library(fabricatr)
#' dat <- fabricate(
#'   N = 40,
#'   Y = rpois(N, lambda = 4),
#'   Z = rbinom(N, 1, prob = 0.4),
#'   D  = Z * rbinom(N, 1, prob = 0.8),
#'   X = rnorm(N)
#' )
#'
#' # Instrument for treatment `D` with encouragement `Z`
#' tidy(iv_robust(Y ~ D + X | Z + X, data = dat))
#'
#' # Instrument with Stata's `ivregress 2sls , small rob` HC1 variance
#' tidy(iv_robust(Y ~ D | Z, data = dat, se_type = "stata"))
#'
#' # With clusters, we use CR2 errors by default
#' dat$cl <- rep(letters[1:5], length.out = nrow(dat))
#' tidy(iv_robust(Y ~ D | Z, data = dat, clusters = cl))
#'
#' # Again, easy to replicate Stata (again with `small` correction in Stata)
#' tidy(iv_robust(Y ~ D | Z, data = dat, clusters = cl, se_type = "stata"))
#'
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

  if (is.null(model_data$instrument_matrix)) {
    stop(
      "Must specify a `formula` with both regressors and instruments. For ",
      "example, `formula = y ~ x1 + x2 | x1 + z2` where x1 and x2 are the ",
      "regressors and z1 and z2 are the instruments.\n\nSee ?iv_robust."
    )
  }

  if (ncol(model_data$instrument_matrix) < ncol(model_data$design_matrix))  {
    warning("More regressors than instruments")
  }

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
