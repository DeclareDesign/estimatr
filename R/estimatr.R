#' estimatr
#'
#' @docType package
#' @useDynLib estimatr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats sd var model.matrix.default pt qt var weighted.mean lm
#' vcov model.frame.default model.response complete.cases terms reformulate
#' update model.extract setNames delete.response .checkMFClasses model.frame
#' model.matrix na.pass nobs coef pf
#' @importFrom methods setGeneric setMethod isGeneric className
#' @name estimatr
NULL
