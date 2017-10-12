#' Confidence intervals for  \link{\code{lm_robust}} objects
#'
#' @param obj an object of class 'lm_robust'
#' @param parm a specification of which parameters are to be given confidence
#' intervals, either a vector of numbers or a vector of names. Defaults to
#' 'coefficient_name' passed to  \link{\code{lm_robust}}. Returns all if missing.
#' @param level the significance level, defaults to 'alpha' specified in \link{\code{lm_robust}}
#'
#' @export
confint.lm_robust <-
  function(
    obj,
    parm = NULL,
    level = NULL
  ) {

    cis <- get_ci_mat(obj, level)

    if (!is.null(parm)) {
      cis <- cis[parm, , drop = F]
    } else if (!is.null(obj$which_covs)) {
      cis <- cis[obj$which_covs, , drop = F]
    }

    return(cis)

  }

#' Confidence intervals for \link{\code{difference_in_means}}objects
#'
#' @param obj an object of class 'difference_in_means'
#' @param level the significance level, defaults to 'alpha' specified in \link{\code{difference_in_means}}
#'
#' @export
confint.difference_in_means <-
  function(
    obj,
    level = NULL
  ) {

    cis <- get_ci_mat(obj, level)

    return(cis)

  }

#' Confidence intervals for \link{\code{horvitz_thompson}}objects
#'
#' @param obj an object of class 'horvitz_thompson'
#' @param level the significance level, defaults to 'alpha' specified in \link{\code{horvitz_thompson}}
#'
#' @export
confint.horvitz_thompson <-
  function(
    obj,
    level = NULL
  ) {

    cis <- get_ci_mat(obj, level)

    return(cis)

  }

## internal method that builds confidence intervals and labels the matrix to be returned
get_ci_mat <- function(obj, level) {
  if (!is.null(level)) {
    ci_lower <- obj$est - qt(1 - level / 2, df = obj$df) * obj$se
    ci_upper <- obj$est + qt(1 - level / 2, df = obj$df) * obj$se
    cis <- cbind(ci_lower, ci_upper)
  } else {
    cis <- cbind(obj$ci_lower, obj$ci_upper)
    level <- obj$alpha
  }

  dimnames(cis) <-
    list(
      obj$coefficient_name,
      paste(level / 2 * c(100, -100) + c(0, 100), '%')
    )

  return(cis)
}
