# This code modified from
# https://github.com/atahk/bucky/blob/master/R/zzz.R (GPL 3.0)
.onLoad <- function(libname, pkgname) {
  if (suppressWarnings(requireNamespace("texreg", quietly=TRUE))) {
    setGeneric("extract", function(model, ...) standardGeneric("extract"),
               package = "texreg")
    setMethod("extract",
              signature = className("lm_robust", pkgname),
              definition = extract.lm_robust)
  }
  invisible()
}






















