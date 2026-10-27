# lm's meanings. A fit stores its residuals on the data's scale, so the
# weighted residual sum of squares is sum(w * e^2).
#' @importFrom stats deviance
#' @export
deviance.lm_robust <- function(object, ...) {
  e <- as.matrix(object[["residuals"]])
  if (ncol(e) > 1) {
    stop("`deviance()` is not available for a fit with several outcomes.",
         call. = FALSE)
  }
  w <- object[["weights"]]
  sum(if (is.null(w)) e^2 else w * e^2)
}

#' @export
deviance.iv_robust <- deviance.lm_robust

#' @importFrom stats sigma
#' @export
sigma.lm_robust <- function(object, ...) {
  sqrt(deviance(object) / object[["df.residual"]])
}

#' @export
sigma.iv_robust <- sigma.lm_robust

# The row names of the data rows the fit used; `full = FALSE` drops those with
# zero weight, as case.names.lm() does.
#' @importFrom stats case.names
#' @export
case.names.lm_robust <- function(object, full = FALSE, ...) {
  dn <- rownames(model.frame(object))
  w <- object[["weights"]]
  if (full || is.null(w)) dn else dn[w != 0]
}

#' @export
case.names.iv_robust <- case.names.lm_robust
