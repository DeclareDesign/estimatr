coef_estimatr <- function(object, ...) {
  object$estimate
}

#' @export
coef.lm_robust <- coef_estimatr

#' @export
coef.iv_robust <- coef_estimatr

#' @export
coef.lm_lin <- coef_estimatr

#' @export
coef.horvitz_thompson <- coef_estimatr

#' @export
coef.difference_in_means <- coef_estimatr
