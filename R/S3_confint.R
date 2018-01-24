#' Confidence intervals for  \code{\link{lm_robust}} objects
#'
#' @param object an object of class \code{\link{lm_robust}}
#' @param parm a specification of which parameters are to be given confidence
#' intervals, either a vector of numbers or a vector of names. Returns all if missing.
#' @param level the significance level, defaults to 'alpha' specified in \code{\link{lm_robust}}
#' @param ... other arguments, unused
#'
#' @export
confint.lm_robust <-
  function(
           object,
           parm = NULL,
           level = NULL,
           ...) {
    cis <- get_ci_mat(object, level)

    if (!is.null(parm)) {
      cis <- cis[parm, , drop = F]
    }

    return(cis)
  }

#' Confidence intervals for \code{\link{difference_in_means}}objects
#'
#' @param object an object of class \code{\link{difference_in_means}}
#' @param parm unused
#' @param level the significance level, defaults to 'alpha' specified in \code{\link{difference_in_means}}
#' @param ... other arguments, unused
#'
#' @export
confint.difference_in_means <-
  function(
           object,
           parm = NULL,
           level = NULL,
           ...) {
    cis <- get_ci_mat(object, level)

    return(cis)
  }

#' Confidence intervals for \code{\link{horvitz_thompson}} objects
#'
#' @param object an object of class \code{\link{horvitz_thompson}}
#' @param parm unused
#' @param level the significance level, defaults to 'alpha' specified in \code{\link{horvitz_thompson}}
#' @param ... other arguments, unused
#'
#' @export
confint.horvitz_thompson <-
  function(
           object,
           parm = NULL,
           level = NULL,
           ...) {
    cis <- get_ci_mat(object, level)

    return(cis)
  }

## internal method that builds confidence intervals and labels the matrix to be returned
get_ci_mat <- function(object, level) {
  if (!is.null(level)) {
    ci_lower <- object$coefficients - qt(1 - level / 2, df = object$df) * object$se
    ci_upper <- object$coefficients + qt(1 - level / 2, df = object$df) * object$se
    cis <- cbind(ci_lower, ci_upper)
  } else {
    cis <- cbind(object$ci_lower, object$ci_upper)
    level <- object$alpha
  }

  dimnames(cis) <-
    list(
      object$coefficient_name,
      paste(level / 2 * c(100, -100) + c(0, 100), "%")
    )

  return(cis)
}
