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
#' @importFrom Formula as.Formula
#' @importFrom rlang enquos enquo eval_tidy quo_get_expr quo_set_expr quo_is_missing sym quo
#' @name estimatr
NULL
