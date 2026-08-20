#' Two-Stage Least Squares Instrumental Variables Regression
#'
#' @param formula an object of class formula with regressors and instruments,
#'   e.g. `y ~ x1 + x2 | z1 + z2`.
#' @param data A `data.frame`
#' @param weights the bare (unquoted) name of the weights variable
#' @param subset An optional bare (unquoted) expression specifying a subset
#' @param clusters An optional bare (unquoted) name of the cluster variable
#' @param fixed_effects An optional one-sided formula of fixed effects to absorb,
#'   such as `~ blockID`. Uses FWL demeaning (see [lm_robust()] for details and
#'   SE type restrictions). Diagnostics are not available with `fixed_effects`.
#' @param se_type The standard error type. `"HC2"`, `"HC3"` and `"CR2"` are
#'   **not available** with `fixed_effects` here. [lm_robust()] lifts that
#'   restriction for a single FE factor via a leverage identity, but the
#'   identity has not been derived for the two-stage case, where the second
#'   stage runs on fitted rather than observed regressors, so it is not assumed.
#'   Defaults: `"HC2"` (no clusters, no FE), `"CR2"` (clusters, no FE),
#'   `"stata"` (no clusters, with FE), `"CR0"` (clusters, with FE).
#' @param ci logical. Whether to compute p-values and confidence intervals.
#' @param alpha The significance level, 0.05 by default.
#' @param diagnostics logical. Whether to compute IV diagnostic statistics.
#' @param return_vcov logical. Whether to return the vcov matrix.
#' @param try_cholesky logical. Whether to try Cholesky decomposition.
#'
#' @return An object of class `"iv_robust"`.
#'
#' @importFrom stats na.omit
#' @examples
#' set.seed(25)
#' n <- 200
#' dat <- data.frame(z = rbinom(n, 1, 0.5), cl = rep(1:20, each = 10))
#' dat$x <- dat$z * rbinom(n, 1, 0.7)
#' dat$y <- dat$x + rnorm(n)
#'
#' # Endogenous regressor on the left of the bar, instrument on the right
#' fit <- iv_robust(y ~ x | z, data = dat)
#' tidy(fit)
#'
#' # The same variance menu as lm_robust()
#' iv_robust(y ~ x | z, data = dat, se_type = "classical")
#' iv_robust(y ~ x | z, data = dat, clusters = cl)
#'
#' # Weak-instrument, endogeneity and overidentification tests
#' summary(iv_robust(y ~ x | z, data = dat, diagnostics = TRUE))
#'
#' @export
iv_robust <- function(formula,
                      data,
                      weights,
                      subset,
                      clusters,
                      fixed_effects,
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
    cluster = clusters,
    fixed_effects = fixed_effects
  )
  data <- rlang::enquo(data)
  model_data <- clean_model_data(data = data, datargs, estimator = "iv")

  if (ncol(model_data$instrument_matrix) < ncol(model_data$design_matrix)) {
    warning("More regressors than instruments")
  }

  has_fe  <- is.character(model_data[["fixed_effects"]])
  fe_rank <- 0L
  yoriginal <- NULL

  if (has_fe) {
    if (diagnostics) {
      warning("Diagnostics are not available with `fixed_effects`. Skipping.")
      diagnostics <- FALSE
    }
    yoriginal  <- as.matrix(model_data[["outcome"]])
    model_data <- demean_fes(model_data)
    model_data[["instrument_matrix"]] <- demean_matrix_by_fes(
      model_data[["instrument_matrix"]], model_data
    )
    fe_rank <- sum(model_data[["fe_levels"]]) - length(model_data[["fe_levels"]]) + 1L
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
      iv_stage = list(1),
      fe_rank = fe_rank
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
      iv_stage = list(2, model_data$design_matrix),
      fe_rank = fe_rank
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
  if (has_fe) {
    for (nm in c("r.squared", "adj.r.squared", "tss", "fstatistic")) {
      if (!is.null(return_list[[nm]])) {
        return_list[[paste0("proj_", nm)]] <- return_list[[nm]]
        return_list[[nm]] <- NULL
      }
    }
    residuals_proj <- drop(return_list[["residuals"]])
    return_list[["fitted.values"]] <- drop(yoriginal) - residuals_proj

    n_obs <- nrow(yoriginal)
    tss_full <- sum((yoriginal - mean(yoriginal))^2)
    r2_full  <- 1 - sum(residuals_proj^2) / tss_full
    return_list[["r.squared"]]     <- r2_full
    return_list[["adj.r.squared"]] <- 1 - (1 - r2_full) * (n_obs - 1L) / return_list[["df.residual"]]
    return_list[["tss"]]           <- tss_full
    return_list[["felevels"]]      <- model_data[["fe_level_names"]]
    return_list[["fixed_effects"]] <- absorbed_group_effects(
      return_list[["fitted.values"]], return_list[["coefficients"]], model_data
    )
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
    firststage_nomdf <- lm_instruments[["rank"]]
    firststage_fstat_value <- lm_instruments[["fstatistic"]][seq_len(length(endog))]
  } else {
    lm_noinstruments <- lm_robust_fit(
      y = model_data$design_matrix[, endog, drop = FALSE],
      X = model_data$instrument_matrix[
        ,
        !(colnames(model_data$instrument_matrix) %in% instruments),
        drop = FALSE
      ],
      weights = model_data$weights,
      cluster = model_data$cluster,
      se_type = "none",
      has_int = FALSE,
      ci = FALSE,
      return_fit = TRUE,
      return_vcov = FALSE
    )

    coef_noinst <- as.matrix(lm_noinstruments[["coefficients"]])
    inst_indices <- which(!(rownames(coef_inst) %in% rownames(coef_noinst)))
    firststage_nomdf <- lm_instruments[["rank"]] - lm_noinstruments[["rank"]]
    firststage_fstat_value <- compute_fstat(
      coef_matrix = coef_inst,
      coef_indices = inst_indices,
      vcov_fit = lm_instruments[["vcov"]],
      rank = lm_instruments[["rank"]],
      nomdf = firststage_nomdf
    )
  }

  fstat_names <- if (ncol(coef_inst) > 1) {
    paste0(colnames(coef_inst), ":value")
  } else {
    "value"
  }

  dendf <- get_dendf(lm_instruments)

  c(
    setNames(firststage_fstat_value, fstat_names),
    nomdf = firststage_nomdf,
    dendf = dendf,
    setNames(
      vapply(
        firststage_fstat_value,
        function(x) pf(x, firststage_nomdf, dendf, lower.tail = FALSE),
        numeric(1)
      ),
      gsub("value", "p.value", fstat_names)
    )
  )
}

wu_hausman_reg_ftest <- function(model_data, first_stage_residuals, se_type) {

  has_int <- 0 %in% attr(model_data$design_matrix, "assign")

  lm_noresids <- lm_robust_fit(
    y = model_data$outcome,
    X = model_data$design_matrix,
    weights = model_data$weights,
    cluster = model_data$cluster,
    se_type = "none",
    has_int = has_int,
    ci = FALSE,
    return_fit = TRUE,
    return_vcov = FALSE
  )

  lm_resids <- lm_robust_fit(
    y = model_data$outcome,
    X = cbind(model_data$design_matrix, first_stage_residuals),
    weights = model_data$weights,
    cluster = model_data$cluster,
    se_type = se_type,
    has_int = has_int,
    ci = FALSE,
    return_fit = TRUE,
    return_vcov = TRUE
  )

  coef_noresids <- na.omit(lm_noresids[["coefficients"]])
  coef_resids   <- na.omit(lm_resids[["coefficients"]])
  ovar <- which(!(names(coef_resids) %in% names(coef_noresids)))
  wu_hausman_nomdf <- lm_resids[["rank"]] - lm_noresids[["rank"]]

  wu_hausman_fstat <- compute_fstat(
    coef_matrix = as.matrix(coef_resids),
    coef_indices = ovar,
    vcov_fit = lm_resids[["vcov"]],
    rank = lm_resids[["rank"]],
    nomdf = wu_hausman_nomdf
  )

  dendf <- get_dendf(lm_resids)

  c(
    value   = wu_hausman_fstat,
    numdf   = wu_hausman_nomdf,
    dendf   = dendf,
    p.value = pf(wu_hausman_fstat, wu_hausman_nomdf, dendf, lower.tail = FALSE)
  )
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
