#' @importFrom Formula Formula
#' @importFrom stats getCall
#' @export
update.iv_robust <- function(object, formula., ..., evaluate = TRUE) {
  if (is.null(call <- getCall(object)))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.))
    call$formula <- formula(update(Formula(formula(object)), formula.))
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate)
    eval(call, parent.frame())
  else call
}
