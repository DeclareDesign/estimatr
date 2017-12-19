
# todo: add return doc
#' Ordinary Least Squares with Robust Standard Errors
#'
#' @param formula an object of class formula, as in \code{\link{lm}}.
#'
#' @param data A data.frame.
#' @param weights the bare (unquoted) names of the weights variable in the supplied data.
#' @param subset An optional bare (unquoted) expression specifying a subset of observations to be used.
#' @param clusters An optional bare (unquoted) name of the variable that corresponds to the clusters in the data.
#' @param se_type The sort of standard error sought. Without clustering: "HC0", "HC1", "HC2" (default), "HC3", or "classical". With clustering: "CR2" (default), "stata".
#' @param ci A boolean for whether to compute and return pvalues and confidence intervals, TRUE by default.
#' @param alpha The significance level, 0.05 by default.
#' @param coefficient_name a character or character vector that indicates which coefficients should be reported. If left unspecified, returns all coefficients.
#' @param return_vcov a boolean for whether to return the vcov matrix for later usage, TRUE by default.
#' @param trychol a boolean for whether to try using a cholesky decomposition to solve LS instead of a QR decomposition, FALSE by default. See 'details'.
#'
#' @export
#'
lm_robust <- function(formula,
                      data,
                      weights,
                      subset,
                      clusters,
                      se_type = NULL,
                      ci = TRUE,
                      alpha = .05,
                      coefficient_name = NULL,
                      return_vcov = TRUE,
                      trychol = FALSE) {

  where <- parent.frame()
  model_data <- eval(substitute(
    clean_model_data(
      formula = formula,
      data = data,
      subset = subset,
      cluster = clusters,
      weights = weights,
      where = where
    )
  ))

  return_list <-
    lm_robust_fit(
      y = model_data$outcome,
      X = model_data$design_matrix,
      weights = model_data$weights,
      cluster = model_data$cluster,
      ci = ci,
      se_type = se_type,
      alpha = alpha,
      coefficient_name = coefficient_name,
      return_vcov = return_vcov,
      trychol = trychol
    )

  return_list <- lm_return(return_list,
                           model_data = model_data,
                           formula = formula)

  return(return_list)
}
