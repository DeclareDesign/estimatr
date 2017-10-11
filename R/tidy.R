# some code taken from the "broom" package
# https://github.com/tidyverse/broom

#' Tidy the result of a estimator into a data.frame
#'
#' @param object An object to be converted into a tidy data.frame
#' @param ... extra arguments
#'
#' @return a data.frame
#'
#' @export
tidy <- function(object, ...) UseMethod("tidy")


#' tidy on a NULL input
#'
#' tidy on a NULL input returns an empty data frame, which means it can be
#' combined with other data frames (treated as "empty")
#'
#' @param object A value NULL
#' @param ... extra arguments (not used)
#'
#' @return An empty data.frame
#'
#' @export
tidy.NULL <- function(object, ...) {
  data.frame()
}


#' Default tidying method
#'
#' In order to be consistent with the broom package, the default tidy
#' method just tries to cast the object as a data frame
#'
#' @param object an object to be tidied
#' @param ... extra arguments (not used)
#'
#' @export
tidy.default <- function(object, ...) {
  warning(paste("No method for tidying an S3 object of class",
                 class(object),
                 ", using as.data.frame"))
  as.data.frame(object)
}

#' Tidying lm_robust output to a data.frame
#'
#' @param object a lm_robust object to be tidied
#' @param ... extra arguments (not used)
#'
#' @return A data.frame with with coefficient names, estimates, standard
#' errors, confidence intervals, p-values, degrees of freedom for the covariates
#' specified by which_covs
#'
#' @export
tidy.lm_robust <- function(object, ...) {
  return_frame <- tidy_data_frame(object)
  if (!is.null(object$which_covs)) {
    return(return_frame[return_frame$coefficient_name %in% object$which_covs, ])
  } else {
    return(return_frame)
  }

}

#' Tidying difference_in_means output to a data.frame
#'
#' @param object a horvitz_thompson object to be tidied
#' @param ... extra arguments (not used)
#'
#' @return A data.frame with with coefficient names, estimates, standard
#' errors, confidence intervals, p-values, degrees of freedom
#'
#' @export
tidy.difference_in_means <- function(object, ...) {
  return_frame <- tidy_data_frame(object)
  return(return_frame)
}

#' Tidying horvitz_thompson output to a data.frame
#'
#' @param object a horvitz_thompson object to be tidied
#' @param ... extra arguments (not used)
#'
#' @return A data.frame with with coefficient names, estimates, standard
#' errors, confidence intervals, p-values, degrees of freedom
#'
#' @export
tidy.horvitz_thompson <- function(object, ...) {
  return_frame <- tidy_data_frame(object)
  return(return_frame)
}


tidy_data_frame <- function(object) {
  return_cols <-
    c("coefficient_name",
      "est",
      "se",
      "p",
      "ci_lower",
      "ci_upper",
      "df"
    )

  return_frame <- as.data.frame(object[return_cols])
}
