# This fn removes calls from function returns to make testing easier
rmcall <- function(obj) {
  if (!is.null(obj[["call"]])) {
    obj[["call"]] <- NULL
  }
  return(obj)
}

# Casts conditions as character objects for equality purposes
condchr <- function(obj) {
  obj[["condition2"]] <- as.character(obj[["condition2"]])
  obj[["condition1"]] <- as.character(obj[["condition1"]])

  obj
}
