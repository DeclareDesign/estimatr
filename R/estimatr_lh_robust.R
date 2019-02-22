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
#' @details
#'
#' This function is a wrapper for \code{\link{lm_robust}} and for
#' \code{\link{car::linearHypothesis}}. It first runs \code{lm_robust} and
#' then passes \code{"lm_robust"} object as an argument to \code{car::linearHypothesis}.
#'
#' @return An object of class \code{"lh_robust"} comprised of two objects: one class \code{"lh_robust"},  and the other a tidy output for the class \code{linearHypothesis}.
#'
#'
#'  \code{"lh_robust"} contains two lists containing that subsequently contain the
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
#' @examples
#'
#' library(fabricatr)
#' dat <- fabricate(
#'   N = 40,
#'   y = rpois(N, lambda = 4),
#'   x = rnorm(N),
#'   z = rbinom(N, 1, prob = 0.4),
#'   clusterID = sample(1:4, 40, replace = TRUE)
#' )
#'
#' # Default variance estimator is HC2 robust standard errors
#' lhro <- lh_robust(y ~ x + z, data = dat, linearHypothesis = "z + 2x = 0")
#'
#' # Also recovers other sorts of standard erorrs just as specified in \code{\link{lm_robust}}
#' lh_robust(y ~ x + z, data = dat, linearHypothesis = "z + 2x = 0", se_types = classical", ")
#' lh_robust(y ~ x + z, data = dat, linearHypothesis = "z + 2x = 0", se_types = "HC0")
#' # Can specify clusters for cluster-robust inference
#' lh_robust(y ~ x + z,
#'           data = dat,
#'           linearHypothesis = "z + 2x = 0",
#'           clusters = clusterID )
#'
#' #  Can tidy() the data in to a data.frame
#' tidy(lmro)
#' # Can use summary() to get more statistics
#' summary(lmro)
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
                      linearHypothesis = NULL
                      ) {
    call <- match.call()
    arguments <- as.list(call)
    arguments[[1]] <- NULL
    arguments$linearHypothesis <- NULL
    model <- do.call(lm_robust, arguments)

  if(is.null(linearHypothesis)){
    warning("No linear hypothesis test performed")
    return(model)}

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

  return_lh <-  list( coefficients =  estimate,
                      std.error =  std.error,
                      statistic = abs(tt),
                      p.value = p.value,
                      alpha = alpha,
                      conf.low  =  ci[,1],
                      conf.high =  ci[,2],
                      df = df,
                      term = linearHypothesis,
                      outcome =  model$outcome)


  attr(return_lh, "linearHypothesis") <- out
  class( return_lh) <- "lh"


  return_list <- list(lm_robust =  model,
                      linearHypothesis =  return_lh )

  class(return_list)  <- "lh_robust"


  return(return_list)

}
