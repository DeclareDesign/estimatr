# Internal method to process data
clean_model_data <- function(formula,
                             data,
                             subset,
                             weights,
                             block,
                             condition_pr,
                             cluster,
                             where) {
  mf <- match.call()
  mf <- mf[c(1, match(names(formals(sys.function())), names(mf), 0L))] # drop formals left missing
  mf[[1]] <- quote(stats::model.frame)
  mf[["where"]] <- NULL # drop the where argument
  mf[["na.action"]] <- quote(estimatr::na.omit_detailed.data.frame)

  # Weights and clusters may be quoted...
  # TODO helper function
  if (!is.null(mf$weights) && is.character(mf[["weights"]])) {
    mf[["weights"]] <- as.symbol(mf[["weights"]])
  }

  # Clusters...
  if (!is.null(mf$cluster) && is.character(mf[["cluster"]])) {
    mf[["cluster"]] <- as.symbol(mf[["cluster"]])
  }

  # Blocks...
  if (!is.null(mf$block) && is.character(mf[["block"]])) {
    mf[["block"]] <- as.symbol(mf[["block"]])
  }

  # condition_prs...
  if (!is.null(mf$condition_pr) && is.character(mf[["condition_pr"]])) {
    mf[["condition_pr"]] <- as.symbol(mf[["condition_pr"]])
  }

  mf <- eval(mf, where)

  local({
    na.action <- attr(mf, "na.action")
    why_omit <- attr(na.action, "why_omit")

    # Todo generalize to all extra components
    if (!is.null(why_omit[["(cluster)"]])) {
      warning(
        "Some observations have missingness in the cluster variable but not ",
        "in the outcome or covariates. These observations have been dropped."
      )
    }

    if (!is.null(why_omit[["(condition_pr)"]])) {
      warning(
        "Some observations have missingness in the condition_pr variable but ",
        "not in the outcome or covariates. These observations have been dropped."
      )
    }

    if (!is.null(why_omit[["(block)"]])) {
      warning(
        "Some observations have missingness in the block variable but not in ",
        "the outcome or covariates. These observations have been dropped."
      )
    }

    if (!is.null(why_omit[["(weights)"]])) {
      warning(
        "Some observations have missingness in the weights variable but not in ",
        "the outcome or covariates. These observations have been dropped."
      )
    }
  })

  # TODO when using . it adds weights and clusters to model!
  ret <- list(
    outcome = model.response(mf, type = "numeric"),
    design_matrix = model.matrix.default(terms(mf), data = mf)
  )

  # Keep the original treatment vector for DiM and HT
  # They will never have a model frame larger than 6 covars
  # so we can add a check that prevents slowing down large
  # lm_robust calls
  if (ncol(mf) < 6) {
    ret[["original_treatment"]] <- mf[, colnames(mf) == all.vars(terms(mf)[[3]])[1]]
  }

  if (!missing(weights)) {
    ret[["weights"]] <- model.extract(mf, "weights")
    if (any(ret[["weights"]] < 0)) {
      stop("`weights` must not be negative")
    }
  }

  if (!missing(cluster)) {
    ret[["cluster"]] <- model.extract(mf, "cluster")
    if (!(class(ret[["cluster"]]) %in% c("factor", "integer"))) {
      ret[["cluster"]] <- as.factor(ret[["cluster"]])
    }
  }

  if (!missing(block)) {
    ret[["block"]] <- model.extract(mf, "block")
  }

  if (!missing(condition_pr)) {
    ret[["condition_pr"]] <- model.extract(mf, "condition_pr")

    if (any(ret[["condition_pr"]] <= 0 | ret[["condition_pr"]] > 1)) {
      stop(
        "`condition_prs` must be a vector of positive values no greater than 1"
      )
    }
  }

  ret[["terms"]] <- attr(mf, "terms")

  return(ret)
}
