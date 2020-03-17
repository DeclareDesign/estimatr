#' estimatr
#'
#' @description Fast procedures for small set of commonly-used, design-appropriate estimators with robust standard errors and confidence intervals. Includes estimators for linear regression, instrumental variables regression, difference-in-means, Horvitz-Thompson estimation, and regression improving precision of experimental estimates by interacting treatment with centered pre-treatment covariates introduced by Lin (2013) <doi:10.1214/12-AOAS583>.
#'
#' @docType package
#' @useDynLib estimatr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats sd var model.matrix.default pt qt var weighted.mean lm
#' vcov model.frame.default model.response complete.cases terms reformulate
#' update model.extract setNames delete.response .checkMFClasses model.frame
#' model.matrix na.pass nobs coef pf .getXlevels df.residual fitted.values
#' formula model.matrix.lm resid weights lm.fit na.omit pchisq printCoefmat
#' residuals
#' @importFrom methods setGeneric setMethod isGeneric className
#' @importFrom Formula as.Formula
#' @importFrom rlang enquos enquo eval_tidy quo_get_expr quo_set_expr quo_is_missing sym quo
#' @name estimatr
NULL
