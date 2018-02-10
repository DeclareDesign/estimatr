#' @export
nobs.lm_robust <- function(object, ...) object$N

#' @export
nobs.summary.lm_robust <- nobs.lm_robust
