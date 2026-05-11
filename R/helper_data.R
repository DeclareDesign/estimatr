#' @importFrom rlang f_rhs %||% quo_is_missing eval_tidy sym
#' @importFrom stats terms model.matrix model.response model.extract
clean_model_data <- function(data, datargs, estimator = "") {

  data <- if (quo_is_missing(data)) NULL else eval_tidy(data)

  mfargs <- Filter(Negate(quo_is_missing), datargs)

  m_formula <- eval_tidy(mfargs[["formula"]])
  m_formula_env <- environment(m_formula)

  stopifnot("`formula` argument must be a formula" = inherits(m_formula, "formula"))

  mfargs <- as.list(mfargs)

  args_ignored <- c("se_type", "fixed_effects")
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

  # fixed_effects: evaluate and store as factor matrix in formula env so
  # model.frame can attach it as an ancillary variable
  if ("fixed_effects" %in% names(mfargs)) {
    name <- sprintf(".__fixed_effects%%%d__", sample.int(.Machine$integer.max, 1))
    fe_formula <- eval_tidy(mfargs[["fixed_effects"]], data = data)
    m_formula_env[[name]] <- sapply(
      stats::model.frame.default(fe_formula, data = data, na.action = NULL),
      FUN = as.factor
    )
    mfargs[["fixed_effects"]] <- sym(name)
  }

  # IV needs Formula (for the | separator); everything else uses plain formula
  # and avoids the ~80µs overhead of Formula::as.Formula + model.frame.Formula
  if (estimator == "iv") {
    mfargs[["formula"]] <- Formula::as.Formula(m_formula)
  } else {
    mfargs[["formula"]] <- m_formula
  }

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

    to_check_if_missing <- c("cluster", "block", "weights", "fixed_effects")

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
  if (!is.null(ret[["cluster"]]) &&
      !is.factor(ret[["cluster"]]) &&
      !is.integer(ret[["cluster"]])) {
    ret[["cluster"]] <- as.factor(ret[["cluster"]])
  }

  ret[["block"]] <- stats::model.extract(mf, "block")

  ret[["fixed_effects"]] <- stats::model.extract(mf, "fixed_effects")
  if (is.character(ret[["fixed_effects"]])) {
    ret[["fixed_effects"]] <- as.matrix(ret[["fixed_effects"]])
  }

  ret[["terms"]] <- attr(mf, "terms")
  dcs <- attr(ret[["terms"]], "dataClasses")
  drop_vars <- c("(block)", "(cluster)", "(fixed_effects)")
  attr(ret[["terms"]], "dataClasses") <- dcs[setdiff(names(dcs), drop_vars)]
  ret[["xlevels"]] <- .getXlevels(ret[["terms"]], mf)

  return(ret)
}

# Demean outcome and design matrix by fixed effects (Frisch-Waugh-Lovell).
# FE demeaning via alternating projections in C++ (demean_cpp in lm_robust_helper.cpp).
# One-way FE converges in exactly 1 iteration; multi-way iterates to eps = 1e-8.
# Returns model_data with demeaned outcome and design matrix (intercept absorbed).
demean_fes <- function(model_data) {
  fe_df     <- as.data.frame(model_data[["fixed_effects"]], stringsAsFactors = TRUE)
  fe_levels <- vapply(fe_df, nlevels, 0L)
  fe_codes  <- lapply(fe_df, as.integer)
  # At this point model_data[["weights"]] is the RAW weight vector (not yet
  # sqrt-transformed — that happens in lm_robust_fit / prep_lm_data).
  # WLS group means use raw weights as the metric.
  w       <- model_data[["weights"]] %||% numeric(0L)
  has_int <- attr(model_data$terms, "intercept")

  model_data[["yoriginal"]] <- as.matrix(model_data[["outcome"]])
  y_dm <- demean_cpp(as.matrix(model_data[["outcome"]]), fe_codes, w)
  colnames(y_dm) <- colnames(model_data[["outcome"]])
  model_data[["outcome"]] <- y_dm

  X         <- model_data[["design_matrix"]]
  keep_cols <- if (has_int) colnames(X) != "(Intercept)" else rep(TRUE, ncol(X))
  X_sub     <- X[, keep_cols, drop = FALSE]
  model_data[["design_matrix"]] <- if (ncol(X_sub) > 0L) {
    X_dm <- demean_cpp(X_sub, fe_codes, w)
    colnames(X_dm) <- colnames(X_sub)
    X_dm
  } else {
    X_sub
  }

  model_data[["fe_levels"]] <- fe_levels
  model_data
}

# Demean an arbitrary matrix by the same FE structure as model_data.
# Used in iv_robust to demean the instrument matrix.
demean_matrix_by_fes <- function(mat, model_data) {
  fe_df    <- as.data.frame(model_data[["fixed_effects"]], stringsAsFactors = TRUE)
  fe_codes <- lapply(fe_df, as.integer)
  w        <- model_data[["weights"]] %||% numeric(0L)
  has_int  <- !is.null(colnames(mat)) && "(Intercept)" %in% colnames(mat)

  orig_colnames <- colnames(mat)
  mat <- demean_cpp(mat, fe_codes, w)
  colnames(mat) <- orig_colnames
  if (has_int) mat <- mat[, colnames(mat) != "(Intercept)", drop = FALSE]
  mat
}
