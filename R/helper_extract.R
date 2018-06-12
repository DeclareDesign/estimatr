# This code modified from
# https://github.com/leifeld/texreg/blob/master/R/extract.R (no LICENSE)
#' Extract model data for \pkg{texreg} package
#' @rdname extract.lm_robust
#'
#' @description Prepares a \code{"lm_robust"} or \code{"iv_robust"} object for the \pkg{texreg}
#' package. This is largely a clone of the \code{extract.lm}
#' method.
#'
#' @param model an object of class \code{\link{lm_robust}} or \code{"iv_robust"}
#' @param include.ci logical. Defaults to TRUE
#' @param include.rsquared logical. Defaults to TRUE
#' @param include.adjrs logical. Defaults to TRUE
#' @param include.nobs logical. Defaults to TRUE
#' @param include.fstatistic logical. Defaults to TRUE
#' @param include.rmse logical. Defaults to TRUE
#' @param ... unused
#'
extract.robust_default <- function(model,
                              include.ci = TRUE,
                              include.rsquared = TRUE,
                              include.adjrs = TRUE,
                              include.nobs = TRUE,
                              include.fstatistic = FALSE,
                              include.rmse = TRUE,
                              ...) {
  s <- tidy(model)

  names <- s[["term"]]
  co <- s[["estimate"]]
  se <- s[["std.error"]]
  pval <- s[["p.value"]]
  cilow <- numeric()
  ciupper <- numeric()
  if (include.ci) {
    cilow <- s[["conf.low"]]
    ciupper <- s[["conf.high"]]
  }

  rs <- model$r.squared # extract R-squared
  adj <- model$adj.r.squared # extract adjusted R-squared
  n <- nobs(model) # extract number of observations

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rsquared) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs) {
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.fstatistic) {
    fstat <- model[["fstatistic"]][[1]]
    gof <- c(gof, fstat)
    gof.names <- c(gof.names, "F statistic")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rmse && !is.null(model[["res_var"]])) {
    rmse <- sqrt(model[["res_var"]])
    gof <- c(gof, rmse)
    gof.names <- c(gof.names, "RMSE")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  tr <- texreg::createTexreg(
    coef.names = names,
    coef = co,
    se = se,
    pvalues = pval,
    ci.low = cilow,
    ci.up = ciupper,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

#' @rdname extract.lm_robust
#'
#' @export
extract.lm_robust <- extract.robust_default

#' @rdname extract.lm_robust
#'
#' @export
extract.iv_robust <- extract.robust_default

