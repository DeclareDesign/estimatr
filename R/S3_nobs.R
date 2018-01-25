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
