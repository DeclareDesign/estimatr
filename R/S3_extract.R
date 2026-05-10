extract.robust_default <- function(model,
                              include.ci = TRUE,
                              include.rsquared = TRUE,
                              include.adjrs = TRUE,
                              include.nobs = TRUE,
                              include.fstatistic = FALSE,
                              include.rmse = TRUE,
                              include.nclusts = TRUE,
                              ...) {
  s <- tidy(model)

  names <- s[["term"]]
  co <- s[["estimate"]]
  se <- s[["std.error"]]
  pval <- s[["p.value"]]
  cilow <- numeric()
  ciupper <- numeric()
  if (include.ci) {
    cilow <- s[["conf.low"]]
    ciupper <- s[["conf.high"]]
  }

  rs <- model$r.squared
  adj <- model$adj.r.squared
  n <- nobs(model)

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rsquared) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs) {
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj.\\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.fstatistic) {
    fstat <- model[["fstatistic"]][[1]]
    gof <- c(gof, fstat)
    gof.names <- c(gof.names, "F statistic")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rmse && !is.null(model[["res_var"]])) {
    rmse <- sqrt(model[["res_var"]])
    gof <- c(gof, rmse)
    gof.names <- c(gof.names, "RMSE")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nclusts && model[["clustered"]]) {
    gof <- c(gof, model[["nclusters"]])
    gof.names <- c(gof.names, "N Clusters")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  tr <- texreg::createTexreg(
    coef.names = names,
    coef = co,
    se = se,
    pvalues = pval,
    ci.low = cilow,
    ci.up = ciupper,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  return(tr)
}

extract.lm_robust <- extract.robust_default
extract.iv_robust <- extract.robust_default
