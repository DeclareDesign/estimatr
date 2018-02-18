# Internal method to process data
clean_model_data <- function(data, datargs) {

  for (da in names(datargs)) {
    if (quo_is_missing(datargs[[da]])) {
      datargs[[da]] <- NULL
    } else {
      if (is.character(quo_get_expr(datargs[[da]]))) {
        datargs[[da]] <- quo_set_expr(
          datargs[[da]],
          sym(quo_get_expr(datargs[[da]]))
        )
      }
    }
  }

  # Get expressions
  mf <- lapply(datargs, quo_get_expr)
  # if data exists, add it
  if (!quo_is_missing(data)) {
    mf[["data"]] <- expr(eval_tidy(data))
  }

  mf[["na.action"]] <- quote(estimatr::na.omit_detailed.data.frame)

  # Get model frame
  mf <- do.call(stats::model.frame, mf)

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
    formula <- quo(attr(terms(mf), "Formula_without_dot"))
  }

  ret <- list(
    outcome = model.response(mf, type = "numeric"),
    design_matrix = model.matrix(terms(eval_tidy(formula), rhs = 1), data = mf)
  )

  if (any(grepl("\\|", formula[[3]]))) {
    ret[["instrument_matrix"]] <- model.matrix(terms(eval_tidy(formula), rhs = 2), data = mf)
  }

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
