#' estimatr: Fast Estimators for Design-Based Inference
#'
#' Fast procedures for a small set of commonly-used, design-appropriate
#' estimators with robust standard errors and confidence intervals. Provides
#' `lm_robust`, `lm_lin`, `iv_robust`, `lh_robust`, `difference_in_means`, and
#' `horvitz_thompson`, with feols-style fixed effects absorption and
#' design-aware Horvitz-Thompson variance. See `vignette("estimatr2.0")` for
#' what changes in version 2.0 and what does not.
#'
#' @docType package
#' @name estimatr-package
#' @useDynLib estimatr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats coef var setNames delete.response model.frame model.matrix na.pass pt qt pf pnorm qnorm pchisq weighted.mean update reformulate formula printCoefmat .checkMFClasses .getXlevels
"_PACKAGE"
