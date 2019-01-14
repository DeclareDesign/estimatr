#' Linear Hypothesis for Ordinary Least Squares with Robust Standard Errors
#'
#' @description This function is a wrapper for \code{\link{lm_robust}} that
#' is useful for estimating linear combination of coefficients.
#'
#' @param formula an object of class formula, as in \code{\link{lm}}
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
#' @param se_type The sort of standard error sought. If \code{clusters} is
#' not specified the options are "HC0", "HC1" (or "stata", the equivalent),
#'  "HC2" (default), "HC3", or
#' "classical". If \code{clusters} is specified the options are "CR0", "CR2" (default), or "stata". Can also specify "none", which may speed up estimation of the coefficients.
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
#' @param linearHypothesis A character string or a matrix specifying combination.
#' See \code{\link{linearHypothesis}} for more details.
#' @param lh_term A character string indicating the name of term to be assigned to the linear combination.
#' @details
#'
#' This function is simply a wrapper for \code{\link{lm_robust}} and for
#' \code{\link{linearHypothesis}}. It first runs \code{lm_robust} and
#' then passes \code{"lm_robust"} object as an argument to \code{linearHypothesis}.
#'
#' @return An object of class \code{"lm_robust"} (see \code{\link{lm_robust}}) and \code{"lh_robust"}.
#'
#'
#' An object of class \code{"lh_robust"} is a list containing the
#' following components:
#'   \item{coefficients}{the estimated linear combination of coefficients}
#'   \item{std.error}{the estimated standard errors of the linear combination}
#'   \item{statistic}{the t-statistic}
#'   \item{df}{the estimated degrees of freedom}
#'   \item{p.value}{the p-values from a two-sided t-test using \code{coefficients}, \code{std.error}, and \code{df}}
#'   \item{conf.low}{the lower bound of the \code{1 - alpha} percent confidence interval}
#'   \item{conf.high}{the upper bound of the \code{1 - alpha} percent confidence interval}
#'   \item{term}{a character vector of coefficient names}
#'   \item{alpha}{the significance level specified by the user}
#' @export
#'
lh_robust <- function(formula,
                      data,
                      weights,
                      subset,
                      clusters,
                      fixed_effects,
                      se_type = NULL,
                      ci = TRUE,
                      alpha = .05,
                      return_vcov = TRUE,
                      try_cholesky = FALSE,
                      linearHypothesis = NULL,
                      lh_term = NULL
                      ) {

  return_list <- list()


  datargs <- enquos(
    formula = formula,
    weights = weights,
    subset = subset,
    cluster = clusters,
    fixed_effects = fixed_effects
  )
  data <- enquo(data)
  model_data <- clean_model_data(data = data, datargs)

  fes <- is.character(model_data[["fixed_effects"]])
  if (fes) {
    yoriginal <- model_data[["outcome"]]
    Xoriginal <- model_data[["design_matrix"]]
    model_data <- demean_fes(model_data)
    attr(model_data$fixed_effects, "fe_rank") <- sum(model_data[["fe_levels"]]) + 1
  } else {
    Xoriginal <- NULL
    yoriginal <- NULL
  }

  return_list <-
    lm_robust_fit(
      y = model_data$outcome,
      X = model_data$design_matrix,
      yoriginal = yoriginal,
      Xoriginal = Xoriginal,
      weights = model_data$weights,
      cluster = model_data$cluster,
      fixed_effects = model_data$fixed_effects,
      ci = ci,
      se_type = se_type,
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      has_int = attr(model_data$terms, "intercept"),
      iv_stage = list(0)
    )

  return_lm_robust <- lm_return(
    return_list,
    model_data = model_data,
    formula = formula
  )


  model <- return_lm_robust

  out <- car::linearHypothesis(model, linearHypothesis,
                               level = 1-alpha)

  estimate  <- drop(attr(out, "value"))
  std.error <- sqrt(diag(attr(out, "vcov")))
  df <- model$df.residual
  tt <- estimate/std.error
  p.value <-  2 * pt(abs(tt), df)
  alpha_ <- alpha/2
  alpha_ <- c(alpha_, 1 - alpha_)
  ci <- estimate + std.error %o% qt(alpha_, df)
  if(is.null(lh_term)) lh_term <- linearHypothesis

  return_lh_robust <- list(   coefficients =  estimate,
                                std.error =  std.error,
                                statistic = abs(tt),
                                p.value = p.value,
                                alpha = alpha,
                                conf.low  =  ci[,1],
                                conf.high =  ci[,2],
                                df = df,
                                term = lh_term,
                                outcome =  return_lm_robust$outcome)



  attr( return_lh_robust , "class") <- "lh_robust"


  return_list <- list(lm_robust =  return_lm_robust,
                      lh_robust = return_lh_robust )

  attr( return_lh_robust , "class") <- "lh_wrapper"

  return(return_list)

}
