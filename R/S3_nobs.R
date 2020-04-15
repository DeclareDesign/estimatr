#' @export
nobs.lm_robust <- function(object, ...) object$N

#' @export
nobs.iv_robust <- function(object, ...) object$N

#' @export
nobs.difference_in_means <- function(object, ...) object$N

#' @export
nobs.horvitz_thompson <- function(object, ...) object$N

#' @export
nobs.summary.lm_robust <- nobs.lm_robust
