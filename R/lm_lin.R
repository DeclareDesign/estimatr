
#' Linear regression with the Lin (2013) covariate adjustment
#'
#' @param formula An object of class "formula", such as Y ~ Z. Should only have the outcome and the treatment.
#'
#' @param data A data.frame.
#' @param covariates A one-sided formula with all of the covariates on the right hand side, such as ~ x1 + x2 + x3.
#' @param weights the bare (unquoted) names of the weights variable in the supplied data.
#' @param subset An optional bare (unquoted) expression specifying a subset of observations to be used.
#' @param cluster_variable_name An optional bare (unquoted) name of the variable that corresponds to the clusters in the data.
#' @param se_type The sort of standard error sought. Without clustering: "HC0", "HC1", "HC2" (default), "HC3", or "classical". With clustering: "BM" (default), "stata".
#' @param ci A boolean for whether to compute and return pvalues and confidence intervals, TRUE by default.
#' @param alpha The significance level, 0.05 by default.
#' @param coefficient_name a character or character vector that indicates which coefficients should be reported. If left unspecified, returns all coefficients.
#' @param return_vcov a boolean for whether to return the vcov matrix for later usage, TRUE by default.
#'
#' @export
#'
lm_lin <- function(formula,
                   data,
                   covariates,
                   weights,
                   subset,
                   cluster_variable_name,
                   se_type = NULL,
                   ci = TRUE,
                   alpha = .05,
                   coefficient_name = NULL,
                   return_vcov = TRUE) {

  ## Check formula
  if (length(all.vars(formula[[3]])) > 1) {
    stop(
      "The formula should only include one variable on the right-hand side: the treatment variable."
    )
  }

  cov_terms <- terms(covariates)
  ## check covariates is right hand sided fn
  if (attr(cov_terms, "response") != 0) {
    stop(
      "The covariate formula should be a right-hand sided equation, such as '~ x1 + x2 + x3'"
    )
  }
  cov_names <- all.vars(covariates)

  # Get all variables for the design matrix
  full_formula <- update(formula, reformulate(c('.', cov_names), "."))

  model_data <- eval.parent(substitute(
    clean_model_data(
      formula = full_formula,
      data = data,
      subset = subset,
      cluster = cluster_variable_name,
      weights = weights
    )
  ))

  outcome <- model_data$outcome
  design_matrix <- model_data$design_matrix
  weights <- model_data$weights
  cluster <- model_data$cluster

  # If Z is a factor, can't use variable name
  # So get first column non intercept column
  treat_col <- which(attr(design_matrix, "assign") == 1)
  treat_name <- colnames(design_matrix)[treat_col]
  treatment <- design_matrix[, treat_col]

  if (any(!(treatment %in% c(0, 1)))) {
    stop(
      "Treatment variable must be binary for the Lin estimator."
    )
  }

  # center all covariates
  demeaned_covars <-
    scale(
      design_matrix[
        ,
        setdiff(colnames(design_matrix), c(treat_name, '(Intercept)')),
        drop = F
      ],
      center = TRUE,
      scale = FALSE
    )

  # Change name of centered covariates to end in bar
  colnames(demeaned_covars) <- paste0(colnames(demeaned_covars), '_bar')

  # Interacted
  interacted_covars <- treatment * demeaned_covars
  colnames(interacted_covars) <- paste0(treat_name, ':', colnames(demeaned_covars))

  # Interact with treatment
  X <- cbind(design_matrix[, attr(design_matrix, "assign") <= 1, drop = F],
             demeaned_covars,
             interacted_covars)

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
      return_vcov = return_vcov
    )

  return(return_list)
}
