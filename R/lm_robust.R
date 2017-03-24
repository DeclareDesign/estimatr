

#' @export
lm_robust_se <- function(formula,
                         data = data,
                         alpha = .05,
                         se_type = "HC2") {

  design_matrix <- model.matrix.default(formula, data = data)
  variable_names <- colnames(design_matrix)

  outcome <- data[, all.vars(formula[[2]])]
  fit <- lm_robust_helper(y = outcome, X = design_matrix, type = type)

  N <- nrow(design_matrix)
  k <- ncol(design_matrix)
  df <- N - k

  coef <- fit$beta_hat
  se <- sqrt(diag(fit$Vcov_hat))

  p <- 2 * pt(abs(coef), df = df, lower.tail = FALSE)
  ci_lower <- coef - qt(1 - alpha / 2, df = df) * se
  ci_upper <- coef + qt(1 - alpha / 2, df = df) * se

  data.frame(
    variable_names = variable_names,
    est = coef,
    se = se,
    p = p,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )

}
