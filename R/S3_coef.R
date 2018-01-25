#' @export
coef.lm_robust <-
  function(
           object,
           ...) {
    coefs <- object$coefficients

    return(coefs)
  }

#' @export
coef.difference_in_means <-
  function(
           object,
           ...) {
    coefs <- object$coefficients

    return(coefs)
  }

#' @export
coef.horvitz_thompson <-
  function(
           object,
           ...) {
    coefs <- object$coefficients

    return(coefs)
  }
