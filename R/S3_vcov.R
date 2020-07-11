#' @export
vcov.lm_robust <- function(object, complete = TRUE, ...) {
  vcov_simple(object, complete = complete)
}

#' @export
vcov.iv_robust <- vcov.lm_robust

#' @export
vcov.difference_in_means <- function(object, ...) {
  return(object$vcov)
}

#' @export
vcov.horvitz_thompson <- vcov.difference_in_means


# Helper function for extracting vcov when it is just an element in the object list
vcov_simple <- function(object, complete) {
  if (is.null(object$vcov)) {
    stop(
      "Object must have vcov matrix. Try setting `return_vcov = TRUE` in ",
      "the estimator function."
    )
  }
  if (complete && (object$rank < object$k)) {
    vc <- matrix(NA_real_, object$k, object$k,
                 dimnames = list(object$term, object$term))
    j <- which(!is.na(coef(object, complete = TRUE)))
    vc[j, j] <- object$vcov
    return(vc)
  } else {
    return(object$vcov)
  }
}
