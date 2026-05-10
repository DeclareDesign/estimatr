#' estimatrZero: Fast Estimators for Design-Based Inference
#'
#' A lean rewrite of estimatr. Provides `lm_robust`, `lm_lin`, `iv_robust`,
#' `lh_robust`, and `difference_in_means` with HC2/CR2 standard errors.
#' Fixed effects and Horvitz-Thompson are intentionally omitted.
#' Weighted R-squared bug fixed.
#'
#' @docType package
#' @name estimatrZero-package
#' @useDynLib estimatrZero, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats coef var setNames delete.response model.frame model.matrix na.pass pt qt pf pnorm qnorm pchisq weighted.mean update reformulate formula printCoefmat .checkMFClasses .getXlevels
"_PACKAGE"
