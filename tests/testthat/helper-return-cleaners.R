# This fn removes calls from function returns to make testing easier
rmcall <- function(obj) {
  if (!is.null(obj[["call"]])) {
    obj[["call"]] <- NULL
  }
  return(obj)
}
