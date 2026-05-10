#' @importFrom rlang f_rhs %||% quo_is_missing eval_tidy sym
#' @importFrom stats terms model.matrix model.response model.extract
clean_model_data <- function(data, datargs, estimator = "") {

  data <- if (quo_is_missing(data)) NULL else eval_tidy(data)

  mfargs <- Filter(Negate(quo_is_missing), datargs)

  m_formula <- eval_tidy(mfargs[["formula"]])
  m_formula_env <- environment(m_formula)

  stopifnot("`formula` argument must be a formula" = inherits(m_formula, "formula"))

  mfargs <- as.list(mfargs)

  args_ignored <- "se_type"
  to_process <- setdiff(
    names(mfargs),
    c(
      setdiff(names(formals(stats::model.frame.default)), "subset"),
      args_ignored
    )
  )

  for (da in to_process) {
    name <- sprintf(".__%s%%%d__", da, sample.int(.Machine$integer.max, 1))
    m_formula_env[[name]] <- eval_tidy(mfargs[[da]], data = data)
    mfargs[[da]] <- sym(name)
  }

  mfargs[["formula"]] <- Formula::as.Formula(m_formula)

  mf <- eval_tidy(rlang::quo((stats::model.frame)(
    !!!mfargs,
    data = data,
    na.action = na_omit_detailed,
    drop.unused.levels = TRUE
  )))

  local({
    na.action <- attr(mf, "na.action")
    why_omit <- attr(na.action, "why_omit")

    missing_warning <- c(
      "Some observations have missingness in the %s variable(s) but not in ",
      "the outcome or covariates. These observations have been dropped."
    )

    to_check_if_missing <- c("cluster", "block", "weights")

    for (x in to_check_if_missing) {
      if (!is.null(why_omit[[sprintf("(%s)", x)]])) {
        warning(sprintf(missing_warning, x))
      }
    }
  })

  if (!is.null(attr(terms(mf), "Formula_without_dot"))) {
    formula <- attr(terms(mf), "Formula_without_dot")
  } else {
    formula <- eval_tidy(mfargs[["formula"]])
  }

  ret <- list(
    outcome = stats::model.response(mf, type = "numeric"),
    design_matrix = stats::model.matrix(terms(formula, rhs = 1), data = mf),
    formula = formula
  )

  if (estimator == "iv") {
    if (length(formula)[2] != 2) {
      stop(
        "Must specify a `formula` with both regressors and instruments. For ",
        "example, `formula = y ~ x1 + x2 | x1 + z2`.\n\nSee ?iv_robust."
      )
    }
    ret[["instrument_matrix"]] <- stats::model.matrix(terms(formula, rhs = 2), data = mf)
    ret[["terms_regressors"]] <- terms(formula, rhs = 1)
  } else if (estimator == "dim") {
    ret[["original_treatment"]] <- mf[, colnames(mf) == all.vars(terms(mf)[[3]])[1]]
  }

  ret[["weights"]] <- stats::model.extract(mf, "weights")
  if (any(ret[["weights"]] < 0)) {
    stop("`weights` must not be negative")
  }

  ret[["cluster"]] <- stats::model.extract(mf, "cluster")
  if (!(class(ret[["cluster"]]) %in% c("factor", "integer")) &&
      !is.null(ret[["cluster"]])) {
    ret[["cluster"]] <- as.factor(ret[["cluster"]])
  }

  ret[["block"]] <- stats::model.extract(mf, "block")

  ret[["terms"]] <- attr(mf, "terms")
  dcs <- attr(ret[["terms"]], "dataClasses")
  drop_vars <- c("(block)", "(cluster)")
  attr(ret[["terms"]], "dataClasses") <- dcs[setdiff(names(dcs), drop_vars)]
  ret[["xlevels"]] <- .getXlevels(ret[["terms"]], mf)

  return(ret)
}
