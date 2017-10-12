#' Summarizing 'lm_robust' objects
#'
#' @param obj an object of class 'lm_robust'
#' @param ... arguments passed to \link{\code{tidy.lm_robust}}, unused
#'
#' @export
summary.lm_robust <-
  function(
    obj,
    ...
  ) {

    tidy(obj, ...)

  }


#' Summarizing 'difference_in_means' objects
#'
#' @param obj an object of class 'difference_in_means'
#' @param ... arguments passed to \link{\code{tidy.difference_in_means}}, unused
#'
#' @export
summary.difference_in_means <-
  function(
    obj,
    ...
  ) {

    tidy(obj, ...)

  }


#' Summarizing 'horvitz_thompson' objects
#'
#' @param obj an object of class 'horvitz_thompson'
#' @param ... arguments passed to \link{\code{tidy.horvitz_thompson}}, unused
#'
#' @export
summary.horvitz_thompson <-
  function(
    obj,
    ...
  ) {

    tidy(obj, ...)

  }
