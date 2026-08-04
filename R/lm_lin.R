#' Linear Regression with Lin (2013) Covariate Adjustment
#'
#' @param formula an object of class formula with only the treatment on the RHS
#' @param covariates a right-sided formula with pre-treatment covariates
#' @param data A `data.frame`
#' @param weights the bare (unquoted) name of the weights variable
#' @param subset An optional bare (unquoted) expression specifying a subset
#' @param clusters An optional bare (unquoted) name of the cluster variable
#' @param se_type The sort of standard error (see [lm_robust()])
#' @param ci logical. Whether to compute p-values and confidence intervals.
#' @param alpha The significance level, 0.05 by default.
#' @param return_vcov logical. Whether to return the vcov matrix.
#' @param try_cholesky logical. Whether to try Cholesky decomposition.
#'
#' @return An object of class `"lm_robust"`.
#'
#' @references Lin, Winston. 2013. "Agnostic Notes on Regression Adjustments to
#'   Experimental Data: Reexamining Freedman's Critique." The Annals of Applied
#'   Statistics 7(1): 295-318. \doi{10.1214/12-AOAS583}.
#'
#' @export
lm_lin <- function(formula,
                   covariates,
                   data,
                   weights,
                   subset,
                   clusters,
                   se_type = NULL,
                   ci = TRUE,
                   alpha = .05,
                   return_vcov = TRUE,
                   try_cholesky = FALSE) {

  if (length(all.vars(rlang::f_rhs(formula))) > 1) {
    stop(
      "The `formula` argument must only have the treatment variable on the ",
      "right-hand side. Covariates should go in the `covariates` argument."
    )
  }

  if (!inherits(covariates, "formula")) {
    stop(
      "The `covariates` argument must be specified as a formula:\n",
      "You passed an object of class ", class(covariates)
    )
  }

  cov_terms <- terms(covariates)

  if (attr(cov_terms, "response") != 0) {
    stop(
      "Must not specify a response variable in `covariates` formula.\n",
      "`covariates` must be a right-sided formula, such as '~ x1 + x2 + x3'"
    )
  }

  if (length(attr(cov_terms, "order")) == 0) {
    stop(
      "`covariates` must have a variable on the right-hand side, not 0 or 1"
    )
  }

  full_formula <- update(
    formula,
    reformulate(c(".", labels(cov_terms)))
  )

  datargs <- rlang::enquos(
    formula = full_formula,
    weights = weights,
    subset = subset,
    cluster = clusters
  )
  data <- rlang::enquo(data)
  model_data <- clean_model_data(data = data, datargs)

  outcome <- as.matrix(model_data$outcome)
  n <- nrow(outcome)
  design_matrix <- model_data$design_matrix
  weights <- model_data$weights
  cluster <- model_data$cluster

  has_intercept <- attr(terms(formula), "intercept")
  treat_col <- which(attr(design_matrix, "assign") == 1)
  treatment <- design_matrix[, treat_col, drop = FALSE]
  design_mat_treatment <- colnames(design_matrix)[treat_col]

  # Kept on the fit so predict() expands a numeric multi-valued treatment
  # against the levels seen here rather than whatever happens to appear in
  # `newdata`. NULL for a binary treatment and for a factor, where the design
  # matrix already carries one column per level.
  treatment_vals <- NULL

  if (any(!(treatment %in% c(0, 1)))) {
    vals <- sort(unique(treatment))
    if (has_intercept) vals <- vals[-1]

    names(vals) <- paste0(colnames(design_matrix)[treat_col], vals)

    treatment <-
      outer(
        drop(treatment),
        vals,
        function(x, y) as.numeric(x == y)
      )
    treatment_vals <- vals
  }

  demeaned_covars <-
    design_matrix[
      ,
      setdiff(colnames(design_matrix), c(design_mat_treatment, "(Intercept)")),
      drop = FALSE
    ]

  if (is.numeric(weights)) {
    center <- apply(demeaned_covars, 2, weighted.mean, weights)
  } else {
    center <- colMeans(demeaned_covars)
  }

  demeaned_covars <- sweep(demeaned_covars, 2, center)

  colnames(demeaned_covars) <- lin_covar_names(colnames(demeaned_covars))

  n_treat_cols <- ncol(treatment)
  n_covars <- ncol(demeaned_covars)

  n_int_covar_cols <- n_covars * (n_treat_cols)
  interacted_covars <- matrix(0, nrow = n, ncol = n_int_covar_cols)
  interacted_covars_names <- character(n_int_covar_cols)
  for (i in 1:n_covars) {
    covar_name <- colnames(demeaned_covars)[i]

    cols <- (i - 1) * n_treat_cols + (1:n_treat_cols)
    interacted_covars[, cols] <- treatment * demeaned_covars[, i]
    interacted_covars_names[cols] <- paste0(colnames(treatment), ":", covar_name)
  }
  colnames(interacted_covars) <- interacted_covars_names

  if (has_intercept) {
    X <- cbind(
      matrix(1, nrow = n, ncol = 1, dimnames = list(NULL, "(Intercept)")),
      treatment,
      demeaned_covars,
      interacted_covars
    )
  } else {
    if (n_treat_cols == 1) {
      X <- cbind(
        treatment,
        demeaned_covars,
        interacted_covars
      )
    } else {
      X <- cbind(
        treatment,
        interacted_covars
      )
    }
  }

  return_list <-
    lm_robust_fit(
      y = outcome,
      X = X,
      weights = weights,
      cluster = cluster,
      ci = ci,
      se_type = se_type,
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      has_int = has_intercept,
      iv_stage = list(0)
    )

  return_list <- lm_return(
    return_list,
    model_data = model_data,
    formula = formula
  )

  # `center` already carries the covariates' original names, which is what
  # predict() indexes `newdata` by. (There was a setNames() call here whose
  # result was discarded, so it had never done anything.)
  return_list[["scaled_center"]] <- center
  return_list[["treatment_vals"]] <- treatment_vals

  return_list[["call"]] <- match.call()

  return(return_list)
}
