# Some of this is code modified from
# https://github.com/atahk/bucky/blob/master/R/zzz.R (GPL 3.0)
.onLoad <- function(libname, pkgname) {
  if (suppressWarnings(requireNamespace("texreg", quietly = TRUE))) {
    setGeneric("extract", function(model, ...) standardGeneric("extract"),
      package = "texreg"
    )
    setMethod("extract",
      signature = className("lm_robust", pkgname),
      definition = extract.lm_robust
    )
    setMethod("extract",
              signature = className("iv_robust", pkgname),
              definition = extract.iv_robust
    )
  }
  if(requireNamespace("emmeans", quietly = TRUE)) {
    emmeans::.emm_register("lm_robust", pkgname)
  }
  invisible()
}

#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
  if (isNamespaceLoaded("broom") && packageVersion("broom") <= "0.5.0") {
    packageStartupMessage(
      "Warning: the `broom` package version 0.5.0 or earlier is loaded.\n",
      "Please upgrade `broom` or `tidy` methods may not work as expected."
    )
  }
  invisible()
}
