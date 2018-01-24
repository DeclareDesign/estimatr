#' Printing \code{\link{lm_robust}} objects
#'
#' @param x an object of class \code{\link{lm_robust}}
#' @param ... arguments passed to \code{\link{tidy.lm_robust}}, unused
#'
#' @export
print.lm_robust <-
  function(
           x,
           ...) {
    print(tidy(x, ...))
  }


#' Printing \code{\link{difference_in_means}} objects
#'
#' @param x an object of class \code{\link{difference_in_means}}
#' @param ... arguments passed to \code{\link{tidy.difference_in_means}}, unused
#'
#' @export
print.difference_in_means <-
  function(
           x,
           ...) {
    print(paste0("Design: ", x$design))
    print(tidy(x, ...))
  }


#' Printing \code{\link{horvitz_thompson}} objects
#'
#' @param x an object of class \code{\link{horvitz_thompson}}
#' @param ... arguments passed to \code{\link{tidy.horvitz_thompson}}, unused
#'
#' @export
print.horvitz_thompson <-
  function(
           x,
           ...) {
    print(tidy(x, ...))
  }
