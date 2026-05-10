#' Two-Stage Least Squares Instrumental Variables Regression
#'
#' @param formula an object of class formula with regressors and instruments,
#'   e.g. `y ~ x1 + x2 | z1 + z2`.
#' @param data A `data.frame`
#' @param weights the bare (unquoted) name of the weights variable
#' @param subset An optional bare (unquoted) expression specifying a subset
#' @param clusters An optional bare (unquoted) name of the cluster variable
#' @param se_type The sort of standard error (see [lm_robust()])
#' @param ci logical. Whether to compute p-values and confidence intervals.
#' @param alpha The significance level, 0.05 by default.
#' @param diagnostics logical. Whether to compute IV diagnostic statistics.
#' @param return_vcov logical. Whether to return the vcov matrix.
#' @param try_cholesky logical. Whether to try Cholesky decomposition.
#'
#' @return An object of class `"iv_robust"`.
#'
#' @export
iv_robust <- function(formula,
                      data,
                      weights,
                      subset,
                      clusters,
                      se_type = NULL,
                      ci = TRUE,
                      alpha = .05,
                      diagnostics = FALSE,
                      return_vcov = TRUE,
                      try_cholesky = FALSE) {
  datargs <- rlang::enquos(
    formula = formula,
    weights = weights,
    subset = subset,
    cluster = clusters
  )
  data <- rlang::enquo(data)
  model_data <- clean_model_data(data = data, datargs, estimator = "iv")

  if (ncol(model_data$instrument_matrix) < ncol(model_data$design_matrix)) {
    warning("More regressors than instruments")
  }

  # -----------
  # First stage
  # -----------

  has_int <- attr(model_data$terms, "intercept")
  first_stage <-
    lm_robust_fit(
      y = model_data$design_matrix,
      X = model_data$instrument_matrix,
      weights = model_data$weights,
      cluster = model_data$cluster,
      ci = FALSE,
      se_type = "none",
      has_int = has_int,
      alpha = alpha,
      return_fit = TRUE,
      return_vcov = FALSE,
      try_cholesky = try_cholesky,
      iv_stage = list(1)
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
      ci = ci,
      se_type = se_type,
      has_int = attr(model_data$terms, "intercept"),
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      iv_stage = list(2, model_data$design_matrix)
    )


  return_list <- lm_return(
    second_stage,
    model_data = model_data,
    formula = model_data$formula
  )

  se_type <- return_list[["se_type"]]

  # ------
  # diagnostics
  # ------
  if (diagnostics) {

    instruments <- setdiff(
      colnames(model_data$instrument_matrix),
      colnames(model_data$design_matrix)
    )
    endog <- setdiff(
      colnames(model_data$design_matrix),
      colnames(model_data$instrument_matrix)
    )

    first_stage_fits <- first_stage[["fitted.values"]][, endog, drop = FALSE]
    colnames(first_stage_fits) <- paste0("fit_", colnames(first_stage_fits))

    first_stage_residuals <- model_data$design_matrix - first_stage[["fitted.values"]]
    colnames(first_stage_residuals) <- paste0("resid_", colnames(first_stage_residuals))

    wu_hausman_ftest_val <- wu_hausman_reg_ftest(model_data, first_stage_residuals, se_type)

    extra_instruments <- length(instruments) - length(endog)

    if (extra_instruments && is.null(model_data$weights)) {
      ss_residuals <- model_data$outcome - second_stage[["fitted.values"]]

      if (se_type == "classical") {
        overid_chisq_val <- sargan_chisq(model_data, ss_residuals)
      } else {
        overid_chisq_val <- wooldridge_score_chisq(
          model_data = model_data,
          endog = endog,
          instruments = instruments,
          ss_residuals = ss_residuals,
          first_stage_fits = first_stage_fits,
          m = extra_instruments
        )
      }

      overid_chisqtest_val <- c(
        overid_chisq_val,
        extra_instruments,
        pchisq(overid_chisq_val, extra_instruments, lower.tail = FALSE)
      )

    } else {
      overid_chisqtest_val <- c(NA_real_, 0, NA_real_)
    }
    names(overid_chisqtest_val) <- c("value", "df", "p.value")

    first_stage_ftest_val <- first_stage_ftest(model_data, endog, instruments, se_type)

    return_list[["diagnostic_first_stage_fstatistic"]] <- first_stage_ftest_val
    return_list[["diagnostic_endogeneity_test"]] <- wu_hausman_ftest_val
    return_list[["diagnostic_overid_test"]] <- overid_chisqtest_val
  }
  return_list[["call"]] <- match.call()

  return_list[["terms_regressors"]] <- model_data[["terms_regressors"]]
  return_list[["formula"]] <- formula(formula)
  class(return_list) <- "iv_robust"

  return(return_list)
}

get_dendf <- function(lm_fit) {
  if (is.numeric(lm_fit[["nclusters"]])) {
    lm_fit[["nclusters"]] - 1
  } else {
    lm_fit[["df.residual"]]
  }
}

first_stage_ftest <- function(model_data, endog, instruments, se_type) {

  lm_instruments <- lm_robust_fit(
    y = model_data$design_matrix[, endog, drop = FALSE],
    X = model_data$instrument_matrix,
    weights = model_data$weights,
    cluster = model_data$cluster,
    se_type = se_type,
    has_int = 0 %in% attr(model_data$instrument_matrix, "assign"),
    return_fit = TRUE,
    return_vcov = TRUE,
    ci = FALSE
  )
  coef_inst <- as.matrix(lm_instruments[["coefficients"]])

  if (all(colnames(model_data$instrument_matrix) %in% instruments)) {
    indices <- seq_len(nrow(coef_inst))
  } else {
    indices <- which(colnames(model_data$instrument_matrix) %in% instruments)
  }
  nomdf <- length(indices)

  fstat <- compute_fstat(
    coef_matrix = coef_inst,
    coef_indices = indices,
    vcov_fit = lm_instruments[["vcov"]],
    rank = lm_instruments[["rank"]],
    nomdf = nomdf
  )
  dendf <- get_dendf(lm_instruments)

  first_stage_ftest_val <- c(
    setNames(fstat, paste0(endog, ":value")),
    numdf = nomdf,
    dendf = dendf
  )

  first_stage_ftest_val
}

wu_hausman_reg_ftest <- function(model_data, first_stage_residuals, se_type) {

  aug_design <- cbind(model_data$design_matrix, first_stage_residuals)

  wu_hausman_lm <- lm_robust_fit(
    y = model_data$outcome,
    X = aug_design,
    weights = model_data$weights,
    cluster = model_data$cluster,
    se_type = se_type,
    has_int = attr(model_data$terms, "intercept"),
    return_fit = FALSE,
    return_vcov = TRUE,
    ci = FALSE
  )

  n_resid_cols <- ncol(first_stage_residuals)
  endogeneity_indices <- seq(
    ncol(model_data$design_matrix) + 1,
    ncol(aug_design)
  )

  coef_wu_hausman <- as.matrix(wu_hausman_lm[["coefficients"]])
  fstat <- compute_fstat(
    coef_matrix = coef_wu_hausman,
    coef_indices = endogeneity_indices,
    vcov_fit = wu_hausman_lm[["vcov"]],
    rank = wu_hausman_lm[["rank"]],
    nomdf = n_resid_cols
  )
  dendf <- get_dendf(wu_hausman_lm)

  wu_hausman_ftest_val <- c(
    value = fstat,
    numdf = n_resid_cols,
    dendf = dendf,
    p.value = pf(fstat, n_resid_cols, dendf, lower.tail = FALSE)
  )

  wu_hausman_ftest_val
}

sargan_chisq <- function(model_data, ss_residuals) {
  ss_resid_lm <- lm_robust_fit(
    y = ss_residuals,
    X = model_data$instrument_matrix,
    weights = NULL,
    cluster = NULL,
    se_type = "classical",
    has_int = attr(model_data$terms, "intercept"),
    return_fit = FALSE,
    return_vcov = FALSE,
    ci = FALSE
  )

  nrow(model_data$instrument_matrix) * ss_resid_lm[["r.squared"]]
}

wooldridge_score_chisq <- function(model_data,
                                    endog,
                                    instruments,
                                    ss_residuals,
                                    first_stage_fits,
                                    m) {

  aug_instrument_mat <- cbind(
    model_data$instrument_matrix,
    first_stage_fits[, endog[!(endog %in% colnames(model_data$instrument_matrix))], drop = FALSE]
  )

  wooldridge_lm <- lm_robust_fit(
    y = ss_residuals,
    X = aug_instrument_mat,
    weights = NULL,
    cluster = NULL,
    se_type = "classical",
    has_int = attr(model_data$terms, "intercept"),
    return_fit = FALSE,
    return_vcov = FALSE,
    ci = FALSE
  )

  nrow(aug_instrument_mat) * wooldridge_lm[["r.squared"]]
}
