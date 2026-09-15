na_omit_detailed <- function(object) {
  naomitwhy(object, function(x, w) x[w, , drop = FALSE])
}
