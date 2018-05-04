
#' Prepare model fits for stargazer
#'
#' @param ... a list of lm_robust objects
#' @param stat either "se" (the default) or "p"
#'
#' @return a list of vectors of extracted statistics for stargazers
#' @export
#'
#' @examples
#'
#'
#'
starprep <- function(..., stat = "se"){
  fitlist = list(...)
  if(stat == "se") {
    out <- lapply(fitlist, function(x) x$std.error)
  } else if (stat == "p") {
    out <- lapply(fitlist, function(x){
      ps <- x$p.value
      ps["(Intercept)"] <- 1
      ps
    } )
  }
  return(out)
}
