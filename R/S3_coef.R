#' Extract coefficients from \code{\link{lm_robust}} object
#'
#' @param object an object of class 'lm_robust'
#' @param parm a specification of which parameters are to be given confidence
#' intervals, either a vector of numbers or a vector of names. Returns all if missing.
#' @param ... other arguments, unused
#'
#' @export
coef.lm_robust <-
  function(
    object,
    parm = NULL,
    ...
  ) {

    coefs <- object$est
    names(coefs) <- object$coefficient_name

    if (!is.null(parm)) {
      coefs <- coefs[parm]
    }

    return(coefs)
  }

#' Extract estimate from \code{\link{difference_in_means}}
#'
#' @param object an object of class \code{\link{difference_in_means}}
#' @param ... other arguments, unused
#'
#' @export
coef.difference_in_means <-
  function(
    object,
    ...
  ) {

    coefs <- object$est
    names(coefs) <- object$coefficient_name

    return(coefs)
  }

#' Extract estimate from \code{\link{horvitz_thompson}}
#'
#' @param object an object of class \code{\link{horvitz_thompson}}
#' @param ... other arguments, unused
#'
#' @export
coef.horvitz_thompson <-
  function(
    object,
    ...
  ) {

    coefs <- object$est
    names(coefs) <- object$coefficient_name

    return(coefs)
  }
