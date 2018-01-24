#' Return estimated variance covariance matrix of \code{\link{lm_robust}}
#'
#' @param object An object of class \code{\link{lm_robust}}
#' @param ... extra arguments, unused
#'
#' @export
vcov.lm_robust <-
  function(object, ...) {
    return(vcov_simple(object))
  }

#' Return estimated variance covariance matrix of \code{\link{difference_in_means}}
#'
#' @param object An object of class \code{\link{difference_in_means}}
#' @param ... extra arguments, unused
#'
#' @export
vcov.difference_in_means <-
  function(object, ...) {
    stop("vcov not supported for difference_in_means")
  }

#' Return estimated variance covariance matrix of \code{\link{horvitz_thompson}}
#'
#' @param object An object of class \code{\link{horvitz_thompson}}
#' @param ... extra arguments, unused
#'
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
