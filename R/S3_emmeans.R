### Support for emmeans package
# Note: recover_data and emm_basis methods are registered dynamically in zzz.R

#' @importFrom utils getS3method

recover_data.lm_robust <- function(object, ...) {
  # `envir` is not optional. `recover_data` is emmeans' generic, so without it
  # getS3method() searches the caller's path and fails outright whenever
  # emmeans has been loaded rather than attached, which is what
  # `emmeans::emmeans(...)` does.
  data <- getS3method("recover_data", "lm", envir = asNamespace("emmeans"))(object, ...)
  if (object$rank < object$k)
    attr(data, "pass.it.on") <- TRUE
  data
}

emm_basis.lm_robust <- function(object, trms, xlev, grid, ...) {
  bhat <- coef(object)
  n.mult <- ifelse(is.matrix(bhat), ncol(bhat), 1)
  m <- suppressWarnings(model.frame(trms, grid, na.action = na.pass, xlev = xlev))
  X <- model.matrix(trms, m, contrasts.arg = object$contrasts)
  V <- emmeans::.my.vcov(object, ...)

  if (!anyNA(bhat))
    nbasis <- estimability::all.estble
  else {
    desmat <- model.matrix(trms, data = attr(object, "data"))
    nbasis <- estimability::nonest.basis(desmat)
  }
  misc <- list()
  if (n.mult > 1) {
    eye <- diag(n.mult)
    X <- kronecker(eye, X)
    nbasis <- kronecker(eye, nbasis)
    if (is.null(colnames(bhat)))
      colnames(bhat) <- seq_len(n.mult)
    misc$ylevs <- list(rep.meas = colnames(bhat))
    bhat <- as.numeric(bhat)
  }
  dfargs <- list(df = object$df.residual)
  dffun <- function(k, dfargs) dfargs$df
  list(X = X, bhat = bhat, nbasis = nbasis, V = V,
       dffun = dffun, dfargs = dfargs, misc = misc)
}
