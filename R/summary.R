#' Summarizing \code{\link{lm_robust}} objects
#'
#' @param object an object of class \code{\link{lm_robust}}
#' @param ... arguments passed to \code{\link{tidy.lm_robust}}, unused
#'
#' @export
summary.lm_robust <-
  function(
    object,
    ...
  ) {

    tidy(object, ...)

  }


#' Summarizing \code{\link{difference_in_means}} objects
#'
#' @param object an object of class \code{\link{difference_in_means}}
#' @param ... arguments passed to \code{\link{tidy.difference_in_means}}, unused
#'
#' @export
summary.difference_in_means <-
  function(
    object,
    ...
  ) {

    tidy(object, ...)

  }


#' Summarizing \code{\link{horvitz_thompson}} objects
#'
#' @param object an object of class \code{\link{horvitz_thompson}}
#' @param ... arguments passed to \code{\link{tidy.horvitz_thompson}}, unused
#'
#' @export
summary.horvitz_thompson <-
  function(
    object,
    ...
  ) {

    tidy(object, ...)

  }
