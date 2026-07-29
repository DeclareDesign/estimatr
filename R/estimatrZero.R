#' estimatrZero: Fast Estimators for Design-Based Inference
#'
#' A lean rewrite of estimatr targeting the DeclareDesign workflow. Provides
#' `lm_robust`, `lm_lin`, `iv_robust`, `lh_robust`, `difference_in_means`, and
#' `horvitz_thompson`, with feols-style fixed effects absorption and
#' design-aware Horvitz-Thompson variance. See `vignette("estimatrZero")` for
#' what changes relative to estimatr and what does not.
#'
#' @docType package
#' @name estimatrZero-package
#' @useDynLib estimatrZero, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats coef var setNames delete.response model.frame model.matrix na.pass pt qt pf pnorm qnorm pchisq weighted.mean update reformulate formula printCoefmat .checkMFClasses .getXlevels
"_PACKAGE"
