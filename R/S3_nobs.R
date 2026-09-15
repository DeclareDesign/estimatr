#' @importFrom stats nobs
#' @export
nobs.lm_robust <- function(object, ...) object$nobs

#' @export
# An lh_robust fit keeps its observations on the lm_robust inside it; reading
# the top level returned NULL, in 1.0.6 as well.
nobs.lh_robust <- function(object, ...) object$lm_robust$nobs

#' @export
nobs.iv_robust <- function(object, ...) object$nobs

#' @export
nobs.summary.lm_robust <- nobs.lm_robust

#' @export
nobs.difference_in_means <- function(object, ...) object$nobs

#' @export
nobs.horvitz_thompson <- function(object, ...) object$nobs
