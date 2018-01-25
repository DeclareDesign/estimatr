#' @export
vcov.lm_robust <-
  function(object, ...) {
    return(vcov_simple(object))
  }


#' @export
vcov.difference_in_means <-
  function(object, ...) {
    stop("vcov not supported for difference_in_means")
  }


#' @export
vcov.horvitz_thompson <-
  function(object, ...) {
    stop("vcov not supported for horvitz_thompson")
  }


# Helper function for extracting vcov when it is just an element in the object list
vcov_simple <-
  function(object) {
    if (is.null(object$vcov)) {
      stop(
        "Object must have vcov matrix. Try setting `return_vcov = TRUE` in ",
        "the estimator function."
      )
    } else {
      return(object$vcov)
    }
  }
