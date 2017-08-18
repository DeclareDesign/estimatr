
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
      "The covariate formula should be a right hand sided equation, such as '~ x1 + x2 + x3'"
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
  treatment <- design_matrix[, all.vars(formula[[3]])]

  # check speed
  demeaned_interacted_covars <-
    lapply(
      setdiff(colnames(design_matrix), c(all.vars(formula), '(Intercept)')),
      function(x) {
        xbar <- design_matrix[, x] - mean(design_matrix[, x])
        matrix(
          c(xbar, xbar * treatment),
          ncol = 2,
          dimnames = list(
            NULL,
            c(paste0(x, '_bar'), paste0(all.vars(formula[[3]]), ':', x, '_bar'))
          )
        )
      }
    )

  demeaned_interacted_mat <- do.call("cbind", demeaned_interacted_covars)

  lin_formula <-
    update(
      formula,
      reformulate(c('.', colnames(demeaned_interacted_covars)))
    )

  return(lm_robust(formula = lin_formula,
                   data = cbind(data, demeaned_interacted_mat),
                   subset = subset,
                   ...))
}
