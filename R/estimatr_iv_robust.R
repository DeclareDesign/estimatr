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

  datargs <- enquos(
    formula = formula,
    weights = weights,
    subset = subset,
    cluster = clusters
  )
  data <- enquo(data)
  model_data <- clean_model_data(data = data, datargs)

  # -----------
  # First stage
  # -----------

  first_stage <-
    lm_robust_fit(
      y = model_data$design_matrix,
      X = model_data$instrument_matrix,
      weights = model_data$weights,
      cluster = model_data$cluster,
      ci = FALSE,
      se_type = "none",
      has_int = attr(model_data$terms, "intercept"),
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky
    )

  # ------
  # Second stage
  # ------
  colnames(first_stage$fitted.values) <- colnames(model_data$design_matrix)

  second_stage <-
    lm_robust_fit(
      y = model_data$outcome,
      X = first_stage$fitted.values,
      weights = model_data$weights,
      cluster = model_data$cluster,
      ci = TRUE,
      se_type = se_type,
      has_int = attr(model_data$terms, "intercept"),
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      X_first_stage = model_data$design_matrix
    )

  return_list <- lm_return(
    second_stage,
    model_data = model_data,
    formula = formula
  )

  return_list[["call"]] <- match.call()

  return(return_list)
}
