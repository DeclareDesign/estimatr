# Internal method to process data
#' @importFrom rlang f_rhs
clean_model_data <- function(data, datargs) {

  # if data exists, evaluate it
  if (quo_is_missing(data)) {
    data <- quo(NULL)
  }

  mfargs <- Filter(Negate(quo_is_missing), datargs)
  for (da in setdiff(names(mfargs), names(formals(stats::model.frame)))) {
    mfargs[[da]] <- f_rhs(mfargs[[da]]) # replace quo(w) with quote(w)
    if(is.character(mfargs[[da]])) mfargs[[da]] <- sym(mfargs[[da]])
  }

  mfargs[["formula"]] <- quo(Formula::as.Formula(!!mfargs[["formula"]]))
  mfargs[["na.action"]] <- quo(estimatr::na.omit_detailed.data.frame)
  mfargs[["drop.unused.levels"]] <- TRUE

  # Get model frame
  mf <- eval_tidy(quo((stats::model.frame)(data=!!data, !!!mfargs)))

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


  if (!is.null(attr(terms(mf), "Formula_without_dot"))) {
    formula <- attr(terms(mf), "Formula_without_dot")
  } else {
    formula <- eval_tidy(mfargs[["formula"]]) # unwrap quosure => a formula
  }

  ret <- list(
    outcome = model.response(mf, type = "numeric"),
    design_matrix = model.matrix(terms(formula, rhs = 1), data = mf)
  )

  if (any(grepl("\\|", formula[[3]]))) {
    ret[["instrument_matrix"]] <- model.matrix(terms(formula, rhs = 2), data = mf)
  }

  # Keep the original treatment vector for DiM and HT
  # They will never have a model frame larger than 6 covars
  # so we can add a check that prevents slowing down large
  # lm_robust calls
  if (ncol(mf) < 6) {
    ret[["original_treatment"]] <- mf[, colnames(mf) == all.vars(terms(mf)[[3]])[1]]
  }

  if (!is.null(mfargs$weights)) {
    ret[["weights"]] <- model.extract(mf, "weights")
    if (any(ret[["weights"]] < 0)) {
      stop("`weights` must not be negative")
    }
  }

  if (!is.null(mfargs$cluster)) {
    ret[["cluster"]] <- model.extract(mf, "cluster")
    if (!(class(ret[["cluster"]]) %in% c("factor", "integer"))) {
      ret[["cluster"]] <- as.factor(ret[["cluster"]])
    }
  }

  if (!is.null(mfargs$block)) {
    ret[["block"]] <- model.extract(mf, "block")
  }

  if (!is.null(mfargs$condition_pr)) {
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
