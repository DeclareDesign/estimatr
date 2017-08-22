
#' A Wrapper for Lin 2013 Covariate Adjustment for OLS with Robust Standard Errors
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
                   subset = NULL,
                   ...) {

  ## Check formulae
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

  full_formula <- update(formula, reformulate(c('.', cov_names), "."))

  ## Get design matrix
  condition_call <- substitute(subset)

  if (!is.null(condition_call)) {
    r <- eval(condition_call, data)
    data <- data[r,]
  }

  design_matrix <- model.matrix.default(full_formula, data = data)
  # If Z is a factor, can't use variable name, so get first column not named
  # (Intercept); and it should be the first or second column
  treat_col <- which(colnames(design_matrix)[1:2] != '(Intercept)')[1]
  treat_name <- colnames(design_matrix)[treat_col]
  treatment <- design_matrix[, treat_col]

  if (any(!(treatment %in% c(0, 1)))) {
    stop(
      "Treatment variable must be binary for the Lin estimator."
    )
  }

  # check speed
  demeaned_covars <-
    apply(
      design_matrix[, setdiff(colnames(design_matrix), c(treat_name, '(Intercept)'))],
      2,
      function(x) { xbar <- x - mean(x) }
    )

  colnames(demeaned_covars) <- paste0(colnames(demeaned_covars), '_bar')

  lin_formula <-
    update(
      formula,
      reformulate(paste0(treat_name, ' + ',
                         paste0(treat_name, '*', colnames(demeaned_covars))))
    )

  new_data <- data
  new_data[, treat_name] <- treatment

  return(lm_robust(formula = lin_formula,
                   data = cbind(new_data, demeaned_covars),
                   ...))
}
