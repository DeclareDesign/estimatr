
#' Linear regression with the Lin (2013) covariate adjustment
#'
#' @param formula An object of class "formula", such as Y ~ Z. Should only have the outcome and the treatment.
#' @param covariates A one-sided formula with all of the covariates on the right hand side, such as ~ x1 + x2 + x3.
#' @param data A data.frame.
#' @param weights the bare (unquoted) names of the weights variable in the supplied data.
#' @param subset An optional bare (unquoted) expression specifying a subset of observations to be used.
#' @param clusters An optional bare (unquoted) name of the variable that corresponds to the clusters in the data.
#' @param se_type The sort of standard error sought. Without clustering: "HC0", "HC1", "HC2" (default), "HC3", or "classical". With clustering: "BM" (default), "stata".
#' @param ci A boolean for whether to compute and return pvalues and confidence intervals, TRUE by default.
#' @param alpha The significance level, 0.05 by default.
#' @param coefficient_name a character or character vector that indicates which coefficients should be reported. If left unspecified, returns all coefficients.
#' @param return_vcov a boolean for whether to return the vcov matrix for later usage, TRUE by default.
#' @param try_cholesky a boolean for whether to try using a cholesky decomposition to solve LS instead of a QR decomposition, FALSE by default. See 'details'.

#' @export
#'
lm_lin <- function(formula,
                   covariates,
                   data,
                   weights,
                   subset,
                   clusters,
                   se_type = NULL,
                   ci = TRUE,
                   alpha = .05,
                   coefficient_name = NULL,
                   return_vcov = TRUE,
                   try_cholesky = FALSE) {

  # Check formula
  if (length(all.vars(formula[[3]])) > 1) {
    stop(
      "The formula should only include one variable on the right-hand side: the treatment variable."
    )
  }

  cov_terms <- terms(covariates)
  # Check covariates is right hand sided fn
  if (attr(cov_terms, "response") != 0) {
    stop(
      "The covariate formula should be a right-hand sided equation, such as '~ x1 + x2 + x3'"
    )
  }
  cov_names <- all.vars(covariates)

  # Get all variables for the design matrix
  full_formula <- update(formula, reformulate(c(".", cov_names), "."))

  where <- parent.frame()
  model_data <- eval(substitute(
    clean_model_data(
      formula = full_formula,
      data = data,
      subset = subset,
      cluster = clusters,
      weights = weights,
      where = where
    )
  ))

  outcome <- model_data$outcome
  n <- length(outcome)
  design_matrix <- model_data$design_matrix
  weights <- model_data$weights
  cluster <- model_data$cluster

  # Get treatment columns
  has_intercept <- attr(terms(formula), "intercept")
  treat_col <- which(attr(design_matrix, "assign") == 1)
  treatment <- design_matrix[, treat_col, drop = FALSE]
  design_mat_treatment <- colnames(design_matrix)[treat_col]

  # Check case where treatment is not factor and is not binary
  if (any(!(treatment %in% c(0, 1)))) {
    if (ncol(treatment) > 1) {
      stop(
        "Treatment variable must be binary or a factor for the Lin estimator."
      )
    } else {
      # create dummies for non-factor treatment variable

      # Drop out first group if there is an intercept
      vals <- sort(unique(treatment))[if(has_intercept) {-1} else {TRUE}]
      n_treats <- length(vals)
      # TODO warn if too many values?

      treatment_mat <- matrix(
        NA,
        nrow = n,
        ncol = n_treats,
        dimnames = list(NULL,
                        paste0(colnames(design_matrix)[treat_col], vals))
      )

      for (i in 1:n_treats) {
        treatment_mat[, i] <- as.numeric(treatment == vals[i])
      }

      treatment <- treatment_mat

    }

  }

  # center all covariates
  demeaned_covars <-
    scale(
      design_matrix[
        ,
        setdiff(colnames(design_matrix), c(design_mat_treatment, "(Intercept)")),
        drop = FALSE
      ],
      center = TRUE,
      scale = FALSE
    )

  original_covar_names <- colnames(demeaned_covars)

  # Change name of centered covariates to end in bar
  colnames(demeaned_covars) <- paste0(colnames(demeaned_covars), "_bar")

  n_treat_cols <- ncol(treatment)
  n_covars <- ncol(demeaned_covars)

  # Interacted
  #n_int_covar_cols <- n_covars * (n_treat_cols + has_intercept)
  n_int_covar_cols <- n_covars * (n_treat_cols)
  interacted_covars <- matrix(0, nrow = n, ncol = n_int_covar_cols)
  interacted_covars_names <- character(n_int_covar_cols)
  for (i in 1:n_covars) {

    covar_name <- colnames(demeaned_covars)[i]

    cols <- (i - 1) * n_treat_cols + (1:n_treat_cols)
    interacted_covars[, cols] <-  apply(treatment, 2, `*`, demeaned_covars[, i])
    interacted_covars_names[cols] <- paste0(colnames(treatment), ":", covar_name)
  }
  colnames(interacted_covars) <- interacted_covars_names
  #print(interacted_covars)

  if (has_intercept) {
    # Have to manually create intercept if treatment wasn't a factor
    X <- cbind(matrix(1, nrow = n, ncol = 1, dimnames = list(NULL, c("(Intercept)"))),
               treatment,
               demeaned_covars,
               interacted_covars)
  } else {
    # If no intercept, but treatment is only one column, need to add base terms for covariates
    if (n_treat_cols == 1) {
      X <- cbind(treatment,
                 demeaned_covars,
                 interacted_covars)
    } else {
      X <- cbind(treatment,
                 interacted_covars)
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
      coefficient_name = coefficient_name,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky
    )

  return_list <- lm_return(return_list,
                           model_data = model_data,
                           formula = formula)

  return_list[["scaled_center"]] <- attr(demeaned_covars, "scaled:center")
  setNames(return_list[["scaled_center"]], original_covar_names)

  return(return_list)
}
