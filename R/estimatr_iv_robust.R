#' @export
iv_robust <- function(formula,
                      data,
                      weights,
                      subset,
                      clusters,
                      se_type = NULL,
                      ci = TRUE,
                      alpha = .05,
                      return_vcov = TRUE,
                      try_cholesky = FALSE) {

  where <- parent.frame()
  model_data <- eval(substitute(
    clean_model_data(
      formula = Formula::as.Formula(formula),
      data = data,
      subset = subset,
      cluster = clusters,
      weights = weights,
      where = where
    )
  ))

  # -----------
  # First stage
  # -----------
  print(model_data)
  first_return <-
    lm_robust_fit(
      y = model_data$design_matrix,
      X = model_data$instrument_matrix,
      weights = model_data$weights,
      cluster = model_data$cluster,
      ci = FALSE,
      se_type = "none",
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      has_int = attr(model_data$terms, "intercept")
    )

  # print(first_return)
  # ------
  # Second stage
  # ------
  second_return <-
    lm_robust_fit(
      y = model_data$outcome,
      X = first_return$fit,
      weights = model_data$weights,
      cluster = model_data$cluster,
      ci = FALSE,
      se_type = se_type,
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      has_int = attr(model_data$terms, "intercept"),
      Xfirst = model_data$design_matrix
    )

  return_list <- lm_return(
    second_return,
    model_data = model_data,
    formula = formula
  )

  return_list[["call"]] <- match.call()

  return(return_list)
}
