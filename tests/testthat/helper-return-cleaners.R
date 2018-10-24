# This fn removes calls from function returns to make testing easier
rmcall <- function(obj) structure(obj, call = NULL)

# This fn casts tibbles as data.frames for equivalency tests
# TODO implement this everywhere
expect_equivalent_tbl <- function(obj1, obj2) {
  expect_equivalent(as.data.frame(obj1), as.data.frame(obj2))
}
