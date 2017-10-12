#' Extract coefficients from \link{\code{lm_robust}} object
#'
#' @param obj an object of class 'lm_robust'
#' @param parm a specification of which parameters are to be given confidence
#' intervals, either a vector of numbers or a vector of names. Defaults to
#' 'coefficient_name' passed to  \link{\code{lm_robust}}. Returns all if missing.
#' @param ... other arguments, unused
#'
#' @export
coef.lm_robust <-
  function(
    obj,
    parm = NULL,
    ...
  ) {

    coefs <- obj$est
    names(coefs) <- obj$coefficient_name

    if (!is.null(parm)) {
      coefs <- coefs[parm]
    } else if (!is.null(obj$which_covs)) {
      coefs <- coefs[obj$which_covs]
    }

    return(coefs)
  }

#' Extract estimate from \link{\code{difference_in_means}}
#'
#' @param obj an object of class \link{\code{difference_in_means}}
#' @param ... other arguments, unused
#'
#' @export
coef.difference_in_means <-
  function(
    obj,
    ...
  ) {

    coefs <- obj$est
    names(coefs) <- obj$coefficient_name

    return(coefs)
  }

#' Extract estimate from 'horvitz_thompson'
#'
#' @param obj an object of class 'horvitz_thompson'
#' @param ... other arguments, unused
#'
#' @export
coef.horvitz_thompson <-
  function(
    obj,
    ...
  ) {

    coefs <- obj$est
    names(coefs) <- obj$coefficient_name

    return(coefs)
  }
