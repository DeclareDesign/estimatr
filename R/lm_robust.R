

#' Ordinary Least Squares with Robust Standard Errors
#'
#' @param formula an object of class formula, as in \code{\link{lm}}.
#'
#' @param data A data.frame.
#' @param alpha The significance level, 0.05 by default.
#' @param se_type The sort of standard error sought -- "HCO", "HC1", "HC2", "HC3", or "classical". "HC2" by default
#' @param coefficient_name a character or character vector that indicates which coefficients should be reported. Defaults to "Z".
#'
#' @export
#'
lm_robust_se <- function(formula,
                         data,
                         alpha = .05,
                         se_type = "HC2",
                         coefficient_name = "Z") {

  design_matrix <- model.matrix.default(formula, data = data)
  variable_names <- colnames(design_matrix)

  outcome <- data[, all.vars(formula[[2]])]
  fit <- lm_robust_helper(y = outcome, X = design_matrix, type = se_type)

  N <- nrow(design_matrix)
  k <- ncol(design_matrix)
  df <- N - k

  coef <- fit$beta_hat
  se <- sqrt(diag(fit$Vcov_hat))

  p <- 2 * pt(abs(coef / se), df = df, lower.tail = FALSE)
  ci_lower <- coef - qt(1 - alpha / 2, df = df) * se
  ci_upper <- coef + qt(1 - alpha / 2, df = df) * se

  return_frame <-
  data.frame(
    variable_names = variable_names,
    est = coef,
    se = se,
    p = p,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )

    which_coefs <- return_frame$variable_names %in% coefficient_name

  # if ever we can figure out all the use cases in the test....
  # which_coefs <- return_frame$variable_names %in% deparse(substitute(coefficient_name))

  return(return_frame[which_coefs, ])

}
