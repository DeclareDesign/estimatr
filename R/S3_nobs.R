#' @export
nobs.lm_robust <- function(object, ...) object$nobs

#' @export
nobs.lh_robust <- function(object, ...) object$nobs

#' @export
nobs.iv_robust <- function(object, ...) object$nobs

#' @export
nobs.summary.lm_robust <- nobs.lm_robust

#' @export
nobs.horvitz_thompson <- function(object, ...) object$nobs

#' @export
nobs.difference_in_means <- function(object, ...) object$nobs

