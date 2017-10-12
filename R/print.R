#' Printing \link{\code{lm_robust}} objects
#'
#' @param obj an object of class \link{\code{lm_robust}}
#' @param ... arguments passed to \link{\code{tidy.lm_robust}}, unused
#'
#' @export
print.lm_robust <-
  function(
    obj,
    ...
  ) {

    tidy(obj, ...)

  }


#' Printing \link{\code{difference_in_means}} objects
#'
#' @param obj an object of class \link{\code{difference_in_means}}
#' @param ... arguments passed to \link{\code{tidy.difference_in_means}}, unused
#'
#' @export
print.difference_in_means <-
  function(
    obj,
    ...
  ) {

    tidy(obj, ...)

  }


#' Printing \link{\code{horvitz_thompson}} objects
#'
#' @param obj an object of class \link{\code{horvitz_thompson}}
#' @param ... arguments passed to \link{\code{tidy.horvitz_thompson}}, unused
#'
#' @export
print.horvitz_thompson <-
  function(
    obj,
    ...
  ) {

    tidy(obj, ...)

  }
