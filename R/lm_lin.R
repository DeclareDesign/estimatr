
#' A wrapper for \code{lm_robust} that applies the Lin 2013 covariate adjustment
#'
#' @param formula An object of class "formula", such as Y ~ Z. Should only have the outcome and the treatment.
#'
#' @param data A data.frame.
#' @param covariates A one-sided formula with all of the covariates on the right hand side, such as ~ x1 + x2 + x3.
#' @param ... All other arguments that are passed through to \code{\link{lm_robust}}
#'
#' @export
#'
lm_lin <- function(formula,
                   data,
                   covariates,
                   weights = NULL,
                   subset = NULL,
                   cluster_variable_name = NULL,
                   se_type = NULL,
                   ci = TRUE,
                   alpha = .05,
                   coefficient_name = NULL) {

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

  model_data <-
    clean_model_data(formula = full_formula,
                     data = data,
                     condition_call = substitute(subset),
                     cluster_variable_name = substitute(cluster_variable_name),
                     weights = substitute(weights))

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

  return_frame <- lm_fit(y = outcome,
                         design_matrix = X,
                         weights = weights,
                         cluster = cluster,
                         ci = ci,
                         se_type = se_type,
                         alpha = alpha,
                         coefficient_name = coefficient_name)

  return(return_frame)
}
