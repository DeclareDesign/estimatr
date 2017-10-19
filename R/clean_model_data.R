# Internal method to check for missingness on auxiliary variables and warn
find_warn_missing <- function(x, type) {
  x_missing <- is.na(x)

  if (any(x_missing)) {
    warning(
      sprintf(
        "Some observations have missingness in the %s variable. These observations have been dropped.",
        type
      )
    )
  }

  return(x_missing)
}


# Internal method to process data
clean_model_data <- function(formula,
                             data,
                             subset,
                             weights,
                             cluster,
                             where) {

  mf <- match.call()
  mf <- mf[c(1, match(names(formals(sys.function())), names(mf),0L))] #drop formals left missing
  mf[[1]] <- quote(stats::model.frame)
  mf[["where"]] <- NULL # drop the where argument
  mf[["na.action"]] <- quote(estimatr:::na.omit_detailed.data.frame) #TODO fix :::via roxygen

  # Weights and clusters may be quoted...
  # TODO helper function
  if(hasName(mf, "weights") && is.character(mf[['weights']]))
    mf[["weights"]] <- as.symbol(mf[["weights"]])

  # Clusters...
  if(hasName(mf, "cluster") && is.character(mf[['cluster']]))
    mf[["cluster"]] <- as.symbol(mf[["cluster"]])

  mf <- eval(mf, where)

  local({
    na.action <- attr(mf, "na.action")
    why_omit  <- attr(na.action, "why_omit")

    # Todo generalize to all extra components
    if(hasName(why_omit, "(cluster)")){
      warning(
        "Some observations have missingness in the cluster variable but not in the outcome or covariates. These observations have been dropped."
      )
    }

    if(hasName(why_omit, "(weights)")){
      warning(
        "Some observations have missingness in the weights variable but not in the outcome or covariates. These observations have been dropped."
      )
    }
  })


  ret <- list(
    outcome=model.response(mf),
    design_matrix=model.matrix.default(formula, data = mf)
  )

  if(!missing(weights)){
    ret[["weights"]] <- model.extract(mf, "weights")
  }

  if(!missing(cluster)){
    ret[["cluster"]] <- model.extract(mf, "cluster")
    if(is.character(ret[["cluster"]]))
      ret[["cluster"]] <- as.factor(ret[["cluster"]])
  }

  return(ret)
}
