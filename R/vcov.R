#' Return estimated variance covariance matrix of lm_robust
#'
#' @param obj An object of class 'lm_robust'
#'
#' @export
vcov.lm_robust <-
  function(obj) {
    return(vcov_simple(obj))
  }

#' Return estimated variance covariance matrix of difference_in_means
#'
#' @param obj An object of class 'difference_in_means'
#'
#' @export
vcov.difference_in_means <-
  function(obj) {
    stop("vcov not supported for difference_in_means")
  }

#' Return estimated variance covariance matrix of horvitz_thompson
#'
#' @param obj An object of class 'horvitz_thompson'
#'
#' @export
vcov.difference_in_means <-
  function(obj) {
    stop("vcov not supported for horvitz_thompson")
  }


# Helper function for extracting vcov when it is just an element in the object list
vcov_simple <-
  function(obj) {
    if (is.null(obj$vcov)) {
      stop("Object cannot return vcov. Try setting return_vcov = TRUE.")
    } else {
      # todo: should we only keep the which_covs?
      return(obj$vcov)
    }
  }
