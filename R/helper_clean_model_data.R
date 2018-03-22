# library(estimatr)
# f <- function(w) {
#   dat <- data.frame(x = rnorm(10), y = rnorm(10))
#   lm_robust(y ~ x, data = dat, w = w)
# }
# f(NULL)
# f(1:10)


# Internal method to process data
#' @importFrom rlang f_rhs
clean_model_data <- function(data, datargs, estimator = "") {

  # if data exists, evaluate it
  data <- if (quo_is_missing(data)) NULL else eval_tidy(data)

  if (getOption("estimatr.debug.clean_model_data", FALSE)) browser()

  mfargs <- Filter(Negate(quo_is_missing), datargs)

  m_formula <- eval_tidy(mfargs[["formula"]])
  m_formula_env <- environment(m_formula)

  args_ignored <- c("subset", "se_type")
  # For each ... that would go to model.fram .default, early eval, save to formula env, and point to it
  # subset is also non-standard eval
  to_process <- setdiff(
    names(mfargs),
    setdiff(names(formals(stats::model.frame.default)), args_ignored)
  )

  for (da in to_process) {
    name <- sprintf(".__%s%%%d__", da, sample.int(.Machine$integer.max, 1))
    m_formula_env[[name]] <- eval_tidy(mfargs[[da]], data = data)
    mfargs[[da]] <- sym(name)
  }

  mfargs[["formula"]] <- Formula::as.Formula(m_formula)

  # Get model frame
  mf <- eval_tidy(quo((stats::model.frame)(
    !!!mfargs,
    data = data,
    na.action = na.omit_detailed.data.frame,
    drop.unused.levels = TRUE
  )))

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

  if (estimator == "iv") {
    if (length(formula)[2] != 2) {
      stop(
        "Must specify a `formula` with both regressors and instruments. For ",
        "example, `formula = y ~ x1 + x2 | x1 + z2` where x1 and x2 are the ",
        "regressors and z1 and z2 are the instruments.\n\nSee ?iv_robust."
      )
    }
    ret[["instrument_matrix"]] <- model.matrix(terms(formula, rhs = 2), data = mf)
    ret[["terms_regressors"]] <- terms(formula, rhs = 1)
  } else if (estimator %in% c("ht", "dim")) {
    ret[["original_treatment"]] <- mf[, colnames(mf) == all.vars(terms(mf)[[3]])[1]]
  }

  ret[["weights"]] <- model.extract(mf, "weights")
  if (any(ret[["weights"]] < 0)) {
    stop("`weights` must not be negative")
  }

  ret[["cluster"]] <- model.extract(mf, "cluster")
  if (!(class(ret[["cluster"]]) %in% c("factor", "integer")) && !is.null(ret[["cluster"]])) {
    ret[["cluster"]] <- as.factor(ret[["cluster"]])
  }

  ret[["block"]] <- model.extract(mf, "block")

  ret[["condition_pr"]] <- model.extract(mf, "condition_pr")

  if (any(ret[["condition_pr"]] <= 0 | ret[["condition_pr"]] > 1)) {
    stop(
      "`condition_prs` must be a vector of positive values no greater than 1"
    )
  }

  ret[["terms"]] <- attr(mf, "terms")

  return(ret)
}
