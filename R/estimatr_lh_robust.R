#' Linear Hypothesis for Ordinary Least Squares with Robust Standard Errors
#'
#' @description This function fits a linear model with robust standard errors and performs linear hypothesis test.
#' @param ... Other arguments to be passed to  \code{\link{lm_robust}}
#' @param data A \code{data.frame}
#' @param linear_hypothesis A length 1 character string or a matrix specifying combination, to be passed to the hypothesis.matrix argument of car::linearHypothesis. Joint hypotheses are currently not handled by lh_robust.
#' See \code{\link[car]{linearHypothesis}} for more details.
#' @details
#'
#' This function is a wrapper for \code{\link{lm_robust}} and for
#' \code{\link[car]{linearHypothesis}}. It first runs \code{lm_robust} and
#' next passes \code{"lm_robust"} object as an argument to \code{linearHypothesis}.
#' Currently CR2 standard errors are not handled by lh_robust.
#'
#' @return An object of class \code{"lh_robust"} containing the two following components:
#'
#' \item{lm_robust}{an object as returned by \code{lm_robust}.}
#' \item{lh}{A data frame with most of its columns pulled from \code{linearHypothesis}' output.}
#'
#' The only analysis directly performed by \code{lh_robust} is a \code{t-test} for the null hypothesis of no effects of the linear combination of coefficients as specified by the user.
#' All other output components are either extracted from \code{linearHypothesis} or \code{lm_robust}.
#' Note that the estimate returned is the value of the LHS of an equation of the form f(X) = 0. Hyptheses "x - z = 1", "x +1= z + 2" and "x-z-1=0" will all return the value for "x-y-1"
#'
#' The original output returned by \code{linearHypothesis} is added as an attribute under the \code{"linear_hypothesis"} attribute.
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
#' lhro <- lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z + 2x = 0")
#'
#' # The linear hypothesis argument can be specified equivalently as:
#' lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z = 2x")
#' lh_robust(y ~ x + z, data = dat, linear_hypothesis = "2*x +1*z")
#' lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z + 2x = 0")
#'
#' # Also recovers other sorts of standard erorrs just as specified in \code{\link{lm_robust}}
#' lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z + 2x = 0", se_type = "classical")
#' lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z + 2x = 0", se_type =  "HC1")
#'
#' #  Can tidy() main output and subcomponents in to a data.frame
#' lhro <- lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z + 2x = 0")
#' tidy(lhro )
#' tidy(lhro$lm_robust)
#' tidy(lhro$lh)
#'
#' # Can use summary() to get more statistics on the main output and subcomponents.
#' summary(lhro)
#' summary(lhro$lm_robust)
#' summary(lhro$lh)
#'
#' @importFrom rlang quos eval_tidy
#'
#' @export
#'
lh_robust <- function(..., data, linear_hypothesis) {

  requireNamespace("car")

  # fit lm_robust model
  lm_robust_fit <- lm_robust(...,  data = data)

  alpha <- lm_robust_fit$alpha

  # Checks
  # This stop could also be limited to hypotheses involving more than one coefficient
  if(!is.null(lm_robust_fit$se_type) && lm_robust_fit$se_type == "CR2") {
    stop("lh_robust not available for CR2 standard errors")
  }

  if(lm_robust_fit$clustered  && is.null(lm_robust_fit$se_type)) {
    stop("lh_robust not available for CR2 standard errors; please specify CR0")
  }

  if(length(linear_hypothesis) > 1) {
    stop("lh_robust currently implements tests for hypotheses involving linear combinations of variables but not joint hypotheses (for instance X1 = X2, but not X1 = 0 and X2=0")
  }

  # calculate linear hypothesis
  car_lht <- car::linearHypothesis(
    lm_robust_fit, hypothesis.matrix = linear_hypothesis, level = 1 - alpha)

  estimate  <- drop(attr(car_lht, "value"))
  std.error <- sqrt(diag(attr(car_lht, "vcov")))

  if(length(estimate) > 1) {
    stop("lh_robust currently implements tests for hypotheses involving linear combinations of variables but not joint hypotheses (for instance X1 = X2, but not X1 = 0 and X2=0")
  }

  df <- lm_robust_fit$df[1]

  # appropriate when all elements of df are identical
  if(length(lm_robust_fit$df) >0 && var(lm_robust_fit$df > 0)) {
    warning("lh_robust inference may be inaccurate if degrees of freedom vary across coefficients")
  }

  statistic <- estimate / std.error
  p.value <-  2 * pt(abs(statistic), df, lower.tail = FALSE)
  ci <- estimate + std.error %o% qt(c(alpha / 2, 1 - alpha / 2), df)

  return_lh_robust <- data.frame(
    coefficients = estimate,
    std.error = std.error,
    statistic = statistic,
    p.value = p.value,
    alpha = alpha,
    conf.low  = ci[, 1],
    conf.high = ci[, 2],
    df = df,
    term = linear_hypothesis,
    outcome =  lm_robust_fit$outcome
  )

  attr(return_lh_robust, "linear_hypothesis") <- car_lht
  class(return_lh_robust) <- c("lh", "data.frame")

  return_lm_robust <- lm_robust_fit
  return_lm_robust[["call"]] <- match.call()

  return(structure(
    list(lm_robust = return_lm_robust, lh = return_lh_robust),
    class = "lh_robust"
  ))

}
