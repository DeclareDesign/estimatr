# This fn removes calls from function returns to make testing easier
rmcall <- function(obj) {
  if (!is.null(obj[["call"]])) {
    obj[["call"]] <- NULL
  }
  return(obj)
}

# This fn casts tibbles as data.frames for equivalency tests
expect_equivalent_tbl <- function(obj1, obj2) {
  expect_equivalent(as.data.frame(obj1), as.data.frame(obj2))
}

expect_equal_tbl <- function(obj1, obj2) {
  expect_equal(as.data.frame(obj1), as.data.frame(obj2))
}

expect_identical_tbl <- function(obj1, obj2) {
  expect_identical(as.data.frame(obj1), as.data.frame(obj2))
}
