#' Two-Stage Least Squares Instrumental Variables Regression
#'
#' @description This formula estimates an instrumental variables regression
#' using two-stage least squares with a variety of options for robust
#' standard errors
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
#' @param fixed_effects An optional right-sided formula containing the fixed
#' effects that will be projected out of the data, such as \code{~ blockID}. Do not
#' pass multiple-fixed effects with intersecting groups. Speed gains are greatest for
#' variables with large numbers of groups and when using "HC1" or "stata" standard errors.
#' See 'Details'.
#' @param se_type The sort of standard error sought. If `clusters` is
#' not specified the options are "HC0", "HC1" (or "stata", the equivalent),
#'  "HC2" (default), "HC3", or
#' "classical". If `clusters` is specified the options are "CR0", "CR2" (default), or "stata". Can also specify "none", which may speed up estimation of the coefficients.
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
#' regressors and \code{z1}, \code{z2}, and \code{z3} are instruments. Unlike
#' \code{ivreg}, you must explicitly specify all exogenous regressors on
#' both sides of the bar.
#'
#' The default variance estimators are the same as in \code{\link{lm_robust}}.
#' Without clusters, we default to \code{HC2} standard errors, and with clusters
#' we default to \code{CR2} standard errors. 2SLS variance estimates are
#' computed using the same estimators as in \code{\link{lm_robust}}, however the
#' design matrix used are the second-stage regressors, which includes the estimated
#' endogenous regressors, and the residuals used are the difference
#' between the outcome and a fit produced by the second-stage coefficients and the
#' first-stage (endogenous) regressors. More notes on this can be found at
#' \href{https://declaredesign.org/R/estimatr/articles/mathematical-notes.html}{the mathematical appendix}.
#'
#' If \code{fixed_effects} are specified, both the outcome, regressors, and instruments
#' are centered using the method of alternating projections (Halperin 1962; Gaure 2013). Specifying
#' fixed effects in this way will result in large speed gains with standard error
#' estimators that do not need to invert the matrix of fixed effects. This means using
#' "classical", "HC0", "HC1", "CR0", or "stata" standard errors will be faster than other
#' standard error estimators. Be wary when specifying fixed effects that may result
#' in perfect fits for some observations or if there are intersecting groups across
#' multiple fixed effect variables (e.g. if you specify both "year" and "country" fixed effects
#' with an unbalanced panel where one year you only have data for one country).
#'
#' @return An object of class \code{"iv_robust"}.
#'
#' The post-estimation commands functions \code{summary} and \code{\link{tidy}}
#' return results in a \code{data.frame}. To get useful data out of the return,
#' you can use these data frames, you can use the resulting list directly, or
#' you can use the generic accessor functions \code{coef}, \code{vcov},
#' \code{confint}, and \code{predict}.
#'
#' An object of class \code{"iv_robust"} is a list containing at least the
#' following components:
#'   \item{coefficients}{the estimated coefficients}
#'   \item{std.error}{the estimated standard errors}
#'   \item{df}{the estimated degrees of freedom}
#'   \item{p.value}{the p-values from a two-sided t-test using \code{coefficients}, \code{std.error}, and \code{df}}
#'   \item{ci.lower}{the lower bound of the \code{1 - alpha} percent confidence interval}
#'   \item{ci.upper}{the upper bound of the \code{1 - alpha} percent confidence interval}
#'   \item{term}{a character vector of coefficient names}
#'   \item{alpha}{the significance level specified by the user}
#'   \item{se_type}{the standard error type specified by the user}
#'   \item{res_var}{the residual variance}
#'   \item{N}{the number of observations used}
#'   \item{k}{the number of columns in the design matrix (includes linearly dependent columns!)}
#'   \item{rank}{the rank of the fitted model}
#'   \item{vcov}{the fitted variance covariance matrix}
#'   \item{r.squared}{the \eqn{R^2} of the second stage regression}
#'   \item{adj.r.squared}{the \eqn{R^2} of the second stage regression, but penalized for having more parameters, \code{rank}}
#'   \item{fstatistic}{a vector with the value of the second stage F-statistic with the numerator and denominator degrees of freedom}
#'   \item{weighted}{whether or not weights were applied}
#'   \item{call}{the original function call}
#'   \item{fitted.values}{the matrix of predicted means}
#' We also return \code{terms} with the second stage terms and \code{terms_regressors} with the first stage terms, both of which used by \code{predict}. If \code{fixed_effects} are specified, then we return \code{proj_fstatistic}, \code{proj_r.squared}, and \code{proj_adj.r.squared}, which are model fit statistics that are computed on the projected model (after demeaning the fixed effects).
#'
#' @references
#'
#' Gaure, Simon. 2013. "OLS with multiple high dimensional category variables." Computational Statistics & Data Analysis 66: 8-1. \url{http://dx.doi.org/10.1016/j.csda.2013.03.024}
#'
#' Halperin, I. 1962. "The product of projection operators." Acta Scientiarum Mathematicarum (Szeged) 23(1-2): 96-99.
#'
#' @examples
#' library(fabricatr)
#' dat <- fabricate(
#'   N = 40,
#'   Y = rpois(N, lambda = 4),
#'   Z = rbinom(N, 1, prob = 0.4),
#'   D  = Z * rbinom(N, 1, prob = 0.8),
#'   X = rnorm(N),
#'   G = sample(letters[1:4], N, replace = TRUE)
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
#' # We can also specify fixed effects, that will be taken as exogenous regressors
#' # Speed gains with fixed effects are greatests with "stata" or "HC1" std.errors
#' tidy(iv_robust(Y ~ D | Z, data = dat, fixed_effects = ~ G, se_type = "HC1"))
#'
#' @export
iv_robust <- function(formula,
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
  datargs <- enquos(
    formula = formula,
    weights = weights,
    subset = subset,
    cluster = clusters,
    fixed_effects = fixed_effects
  )
  data <- enquo(data)
  model_data <- clean_model_data(data = data, datargs, estimator = "iv")

  if (ncol(model_data$instrument_matrix) < ncol(model_data$design_matrix)) {
    warning("More regressors than instruments")
  }

  fes <- is.character(model_data[["fixed_effects"]])
  if (fes) {
    yoriginal <- model_data[["outcome"]]
    Xoriginal <- model_data[["design_matrix"]]
    model_data <- demean_fes(model_data)
  } else {
    yoriginal <- NULL
    Xoriginal <- NULL
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
      fixed_effects = model_data$fixed_effects,
      ci = FALSE,
      se_type = "none",
      has_int = attr(model_data$terms, "intercept"),
      alpha = alpha,
      return_fit = TRUE,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      iv_stage = list(1)
    )

  # ------
  # Second stage
  # ------
  colnames(first_stage$fitted.values) <- colnames(model_data$design_matrix)
  if (!is.null(model_data$fixed_effects)) {
    attr(model_data$fixed_effects, "fe_rank") <- sum(model_data[["fe_levels"]]) + 1
  }

  second_stage <-
    lm_robust_fit(
      y = model_data$outcome,
      X = first_stage$fitted.values,
      yoriginal = yoriginal,
      Xoriginal = Xoriginal,
      weights = model_data$weights,
      cluster = model_data$cluster,
      fixed_effects = model_data$fixed_effects,
      ci = ci,
      se_type = se_type,
      has_int = attr(model_data$terms, "intercept"),
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      iv_stage = list(2, model_data$design_matrix)
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
