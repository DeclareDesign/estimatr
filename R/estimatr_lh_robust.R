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
#' effects that will be projected car_lht of the data, such as \code{~ blockID}. Do not
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
#' next passes \code{"lm_robust"} object as an argument to \code{car::linearHypothesis}.
#'
#' @return An object of class \code{"lh_robust"} containing the two following components:
#'
#'  \item{lm_robust} {an object as returned by \code{lm_robust()}}
#'  \item{lh} {A data frame with most of its columns pulled from \code{car::linearHypothesis}' output.}
#'
#' This function solely performs a \code{t-test} for the null hypothesis of no effects of the linear combination of coefficients specified by the user.
#' All other components are obtained by calling \code{linearHypothesis()} and \code{lm_robust()}
#'
#' The original output returned by  \code{car::linearHypothesis} is added as an attribute.
#'
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
#' lh_robust(y ~ x + z, data = dat, linearHypothesis = "z + 2x = 0", se_type = "classical")
#' lh_robust(y ~ x + z, data = dat, linearHypothesis = "z + 2x = 0", se_type =  "HC1")
#
# library(generics)
#' #  Can tidy() main output and subcomponents in to a data.frame
#' lhro <- lh_robust(y ~ x + z, data = dat, linearHypothesis = "z + 2x = 0")
#' tidy(lhro )
#' tidy(lhro$lm_robust)
#' tidy(lhro$lh)
#' # Can use summary() to get more statistics on the main output and subcomponents.
#' summary(lhro)
#' summary(lhro$lm_robust)
#' summary(lhro$lh)
#'
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
    model$call <- call

  if(is.null(linearHypothesis)){
    warning("No linear hypothesis test performed")
    return(model)}

  car_lht <- car::linearHypothesis(model, linearHypothesis,
                               level = 1-alpha)

  estimate  <- drop(attr(car_lht, "value"))
  std.error <- sqrt(diag(attr(car_lht, "vcov")))
  df <- model$df.residual
  tt <- estimate/std.error
  p.value <-  2 * pt(abs(tt), df)
  alpha_ <- alpha/2
  alpha_ <- c(alpha_, 1 - alpha_)
  ci <- estimate + std.error %o% qt(alpha_, df)

  return_lh <-  data.frame( coefficients =  estimate,
                      std.error =  std.error,
                      statistic = abs(tt),
                      p.value = p.value,
                      alpha = alpha,
                      conf.low  =  ci[,1],
                      conf.high =  ci[,2],
                      df = df,
                      term = linearHypothesis,
                      outcome =  model$outcome)


  attr(return_lh, "linearHypothesis") <-   car_lht
  class( return_lh) <- c("lh", "data.frame")


  return_list <- list(lm_robust =  model,
                      lh =  return_lh )

  class(return_list)  <- "lh_robust"


  return(return_list)

}
