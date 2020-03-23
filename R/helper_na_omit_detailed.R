#' Extra logging on na.omit handler
#'
#' @param object a data.frame
#'
#' @return a normal \code{omit} object, with the extra attribute \code{why_omit},
#' which contains the leftmost column containing an NA for each row that was dropped, by
#' column name, if any were dropped.
#'
#' @seealso \code{\link{na.omit}}
na.omit_detailed.data.frame <- function(object){

  naomitwhy(object, function(x, w) x[w, , drop=FALSE])

}

