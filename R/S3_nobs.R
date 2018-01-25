#' Extract number of observations from \code{\link{lm_robust}} object
#'
#' @param object an object of class 'lm_robust'
#' @param ... other arguments, unused
#'
#' @export
nobs.lm_robust <-
  function(
    object,
    ...) {

    return(object$N)
  }

#' @export
nobs.summary.lm_robust <-
  function(
    object,
    ...) {

    return(object$N)
  }
