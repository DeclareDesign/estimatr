#' Extra logging on na.omit handler
#'
#' @param object a data.frame
#' @param ... unused
#'
#' @return a normal \code{omit} object, with the extra attribute \code{why_omit},
#' which contains the leftmost column containing an NA for each row that was dropped, by
#' column name, if any were dropped.
#'
#' @seealso \code{\link{na.omit}}
na.omit_detailed.data.frame <- function(object){

  why_omit <- naomitwhy(object, function(x) is.na(x))

  if(is.integer(why_omit)) {
    object <- if(length(dim(object))) object[-why_omit, , drop=FALSE] else object[-why_omit]
    attr(object, "na.action") <- why_omit
  }
  object

}
