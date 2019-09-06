### Support for emmeans package
#
# Note: the recover_data and emm_basis methods are registered dynamically
#   (see zzz.R). So these functions are not exported

#' @importFrom utils getS3method

recover_data.lm_robust <- function(object, ...) {
  data <- getS3method("recover_data", "lm")(object, ...)
  if (object$rank < object$k)  # rank-deficient. Need to pass dataset to emm_basis
    attr(data, "pass.it.on") <- TRUE
  data
}

emm_basis.lm_robust <- function(object, trms, xlev, grid, ...) {
  # coef() works right for lm but coef.aov tosses out NAs
  bhat <- coef(object)
  n.mult <- ifelse(is.matrix(bhat), ncol(bhat), 1)  # columns in mult response
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
  if (n.mult > 1) { # multivariate case. Need to expand some matrices
    eye <- diag(n.mult)
    X <- kronecker(eye, X)
    nbasis <- kronecker(eye, nbasis)
    if(is.null(colnames(bhat)))
      colnames(bhat) <- seq_len(n.mult)
    misc$ylevs <- list(rep.meas = colnames(bhat))
    bhat <- as.numeric(bhat) # stretch coefs into a vector
  }
  dfargs <- list(df = object$df.residual)
  dffun <- function(k, dfargs) dfargs$df
  list(X = X, bhat = bhat, nbasis = nbasis, V = V,
       dffun = dffun, dfargs = dfargs, misc = misc)
}
