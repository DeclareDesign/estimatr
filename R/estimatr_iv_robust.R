#' Two-Stage Least Squares Instrumental Variables Regression
#'
#' @description This formula estimates an instrumental variables regression
#' using two-stage least squares with a variety of options for robust
#' standard errors
#'
#' @param formula an object of class formula of the regression and the instruments.
#' For example, the formula \code{y ~ x1 + x2 | z1 + z2} specifies \code{x1} and \code{x2}
#' as endogenous regressors and \code{z1} and \code{z2} as their respective instruments.
#' @param data A \code{data.frame}
#' @param weights the bare (unquoted) names of the weights variable in the
#' supplied data.
#' @param subset An optional bare (unquoted) expression specifying a subset
#' of observations to be used.
#' @param clusters An optional bare (unquoted) name of the variable that
#' corresponds to the clusters in the data.
#' @param fixed_effects An optional right-sided formula containing the fixed
#' effects that will be projected out of the data, such as \code{~ blockID}. Do not
#' pass multiple-fixed effects with intersecting groups. Speed gains are greatest for
#' variables with large numbers of groups and when using "HC1" or "stata" standard errors.
#' See 'Details'.
#' @param se_type The sort of standard error sought. If \code{clusters} is
#' not specified the options are "HC0", "HC1" (or "stata", the equivalent),
#'  "HC2" (default), "HC3", or
#' "classical". If \code{clusters} is specified the options are "CR0", "CR2" (default), or "stata". Can also specify "none", which may speed up estimation of the coefficients.
#' @param ci logical. Whether to compute and return p-values and confidence
#' intervals, TRUE by default.
#' @param alpha The significance level, 0.05 by default.
#' @param diagnostics logical. Whether to compute and return instrumental variable diagnostic statistics and tests.
#' @param return_vcov logical. Whether to return the variance-covariance
#' matrix for later usage, TRUE by default.
#' @param try_cholesky logical. Whether to try using a Cholesky
#' decomposition to solve least squares instead of a QR decomposition,
#' FALSE by default. Using a Cholesky decomposition may result in speed gains, but should only
#' be used if users are sure their model is full-rank (i.e., there is no
#' perfect multi-collinearity)
#'
#' @details
#'
#' This function performs two-stage least squares estimation to fit
#' instrumental variables regression. The syntax is similar to that in
#' \code{ivreg} from the \code{AER} package. Regressors and instruments
#' should be specified in a two-part formula, such as
#' \code{y ~ x1 + x2 | z1 + z2 + z3}, where \code{x1} and \code{x2} are
#' regressors and \code{z1}, \code{z2}, and \code{z3} are instruments. Unlike
#' \code{ivreg}, you must explicitly specify all exogenous regressors on
#' both sides of the bar.
#'
#' The default variance estimators are the same as in \code{\link{lm_robust}}.
#' Without clusters, we default to \code{HC2} standard errors, and with clusters
#' we default to \code{CR2} standard errors. 2SLS variance estimates are
#' computed using the same estimators as in \code{\link{lm_robust}}, however the
#' design matrix used are the second-stage regressors, which includes the estimated
#' endogenous regressors, and the residuals used are the difference
#' between the outcome and a fit produced by the second-stage coefficients and the
#' first-stage (endogenous) regressors. More notes on this can be found at
#' \href{https://declaredesign.org/r/estimatr/articles/mathematical-notes.html}{the mathematical appendix}.
#'
#' If \code{fixed_effects} are specified, both the outcome, regressors, and instruments
#' are centered using the method of alternating projections (Halperin 1962; Gaure 2013). Specifying
#' fixed effects in this way will result in large speed gains with standard error
#' estimators that do not need to invert the matrix of fixed effects. This means using
#' "classical", "HC0", "HC1", "CR0", or "stata" standard errors will be faster than other
#' standard error estimators. Be wary when specifying fixed effects that may result
#' in perfect fits for some observations or if there are intersecting groups across
#' multiple fixed effect variables (e.g. if you specify both "year" and "country" fixed effects
#' with an unbalanced panel where one year you only have data for one country).
#'
#' If \code{diagnostics} are requested, we compute and return three sets of diagnostics.
#' First, we return tests for weak instruments using first-stage F-statistics (\code{diagnostic_first_stage_fstatistic}). Specifically,
#' the F-statistics reported compare the model regressing each endogeneous variable on both the
#' included exogenous variables and the instruments to a model where each endogenous variable is
#' regressed only on the included exogenous variables (without the instruments). A significant F-test
#' for weak instruments provides evidence against the null hypothesis that the instruments are weak.
#' Second, we return tests for the endogeneity of the endogenous variables, often called the Wu-Hausman
#' test  (\code{diagnostic_endogeneity_test}). We implement the regression test from Hausman (1978), which allows for robust variance estimation.
#' A significant endogeneity test provides evidence against the null that all the variables are exogenous.
#' Third, we return a test for the correlation between the instruments and the error term  (\code{diagnostic_overid_test}). We implement
#' the Wooldridge (1995) robust score test, which is identical to Sargan's (1958) test with classical
#' standard errors. This test is only reported if the model is overidentified (i.e. the number of
#' instruments is greater than the number of endogenous regressors), and if no weights are specified.
#'
#' @return An object of class \code{"iv_robust"}.
#'
#' The post-estimation commands functions \code{summary} and \code{\link{tidy}}
#' return results in a \code{data.frame}. To get useful data out of the return,
#' you can use these data frames, you can use the resulting list directly, or
#' you can use the generic accessor functions \code{coef}, \code{vcov},
#' \code{confint}, and \code{predict}.
#'
#' An object of class \code{"iv_robust"} is a list containing at least the
#' following components:
#'   \item{coefficients}{the estimated coefficients}
#'   \item{std.error}{the estimated standard errors}
#'   \item{statistic}{the t-statistic}
#'   \item{df}{the estimated degrees of freedom}
#'   \item{p.value}{the p-values from a two-sided t-test using \code{coefficients}, \code{std.error}, and \code{df}}
#'   \item{conf.low}{the lower bound of the \code{1 - alpha} percent confidence interval}
#'   \item{conf.high}{the upper bound of the \code{1 - alpha} percent confidence interval}
#'   \item{term}{a character vector of coefficient names}
#'   \item{alpha}{the significance level specified by the user}
#'   \item{se_type}{the standard error type specified by the user}
#'   \item{res_var}{the residual variance}
#'   \item{nobs}{the number of observations used}
#'   \item{k}{the number of columns in the design matrix (includes linearly dependent columns!)}
#'   \item{rank}{the rank of the fitted model}
#'   \item{vcov}{the fitted variance covariance matrix}
#'   \item{r.squared}{the \eqn{R^2} of the second stage regression}
#'   \item{adj.r.squared}{the \eqn{R^2} of the second stage regression, but penalized for having more parameters, \code{rank}}
#'   \item{fstatistic}{a vector with the value of the second stage F-statistic with the numerator and denominator degrees of freedom}
#'   \item{firststage_fstatistic}{a vector with the value of the first stage F-statistic with the numerator and denominator degrees of freedom, useful for a test for weak instruments}
#'   \item{weighted}{whether or not weights were applied}
#'   \item{call}{the original function call}
#'   \item{fitted.values}{the matrix of predicted means}
#' We also return \code{terms} with the second stage terms and \code{terms_regressors} with the first stage terms, both of which used by \code{predict}. If \code{fixed_effects} are specified, then we return \code{proj_fstatistic}, \code{proj_r.squared}, and \code{proj_adj.r.squared}, which are model fit statistics that are computed on the projected model (after demeaning the fixed effects).
#'
#' We also return various diagnostics when \code{`diagnostics` == TRUE}. These are stored in \code{diagnostic_first_stage_fstatistic}, \code{diagnostic_endogeneity_test}, and \code{diagnostic_overid_test}. They have the test statistic, relevant degrees of freedom, and p.value in a named vector. See 'Details' for more. These are printed in a formatted table when the model object is passed to \code{summary()}.
#'
#' @references
#'
#' Gaure, Simon. 2013. "OLS with multiple high dimensional category variables." Computational Statistics & Data Analysis 66: 8-1. \doi{10.1016/j.csda.2013.03.024}
#'
#' Halperin, I. 1962. "The product of projection operators." Acta Scientiarum Mathematicarum (Szeged) 23(1-2): 96-99.
#'
#' @examples
#' library(fabricatr)
#' dat <- fabricate(
#'   N = 40,
#'   Y = rpois(N, lambda = 4),
#'   Z = rbinom(N, 1, prob = 0.4),
#'   D  = Z * rbinom(N, 1, prob = 0.8),
#'   X = rnorm(N),
#'   G = sample(letters[1:4], N, replace = TRUE)
#' )
#'
#' # Instrument for treatment `D` with encouragement `Z`
#' tidy(iv_robust(Y ~ D + X | Z + X, data = dat))
#'
#' # Instrument with Stata's `ivregress 2sls , small rob` HC1 variance
#' tidy(iv_robust(Y ~ D | Z, data = dat, se_type = "stata"))
#'
#' # With clusters, we use CR2 errors by default
#' dat$cl <- rep(letters[1:5], length.out = nrow(dat))
#' tidy(iv_robust(Y ~ D | Z, data = dat, clusters = cl))
#'
#' # Again, easy to replicate Stata (again with `small` correction in Stata)
#' tidy(iv_robust(Y ~ D | Z, data = dat, clusters = cl, se_type = "stata"))
#'
#' # We can also specify fixed effects, that will be taken as exogenous regressors
#' # Speed gains with fixed effects are greatests with "stata" or "HC1" std.errors
#' tidy(iv_robust(Y ~ D | Z, data = dat, fixed_effects = ~ G, se_type = "HC1"))
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
  datargs <- enquos(
    formula = formula,
    weights = weights,
    subset = subset,
    cluster = clusters,
    fixed_effects = fixed_effects
  )
  data <- enquo(data)
  model_data <- clean_model_data(data = data, datargs, estimator = "iv")

  if (ncol(model_data$instrument_matrix) < ncol(model_data$design_matrix)) {
    warning("More regressors than instruments")
  }

  fes <- is.character(model_data[["fixed_effects"]])
  if (fes) {
    if (diagnostics) {
      warning("Will not return `diagnostics` if `fixed_effects` are used.")
      diagnostics <- FALSE
    }
    yoriginal <- model_data[["outcome"]]
    Xoriginal <- model_data[["design_matrix"]]
    model_data <- demean_fes(model_data)
  } else {
    yoriginal <- NULL
    Xoriginal <- NULL
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
      fixed_effects = model_data$fixed_effects,
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
  if (!is.null(model_data$fixed_effects)) {
    attr(model_data$fixed_effects, "fe_rank") <- sum(model_data[["fe_levels"]]) + 1
  }

  second_stage <-
    lm_robust_fit(
      y = model_data$outcome,
      X = first_stage$fitted.values,
      yoriginal = yoriginal,
      Xoriginal = Xoriginal,
      weights = model_data$weights,
      cluster = model_data$cluster,
      fixed_effects = model_data$fixed_effects,
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

    # values for use by multiple diagnostic tests
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

    # Wu-Hausman f-test for endogeneity
    wu_hausman_ftest_val <- wu_hausman_reg_ftest(model_data, first_stage_residuals, se_type)

    # Overidentification test (only computed if n(instruments) > n(endog regressors)
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

    # Weak instrument test (first stage f-test)
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

# IV diagnostic test functions
# helper to get denominator degress of freedom
get_dendf <- function(lm_fit) {
  if (is.numeric(lm_fit[["nclusters"]])) {
    lm_fit[["nclusters"]] - 1
  } else {
    lm_fit[["df.residual"]]
  }
}

# Weak first-stage ftest
first_stage_ftest <- function(model_data, endog, instruments, se_type) {

  lm_instruments <- lm_robust_fit(
    y = model_data$design_matrix[, endog, drop = FALSE],
    X = model_data$instrument_matrix,
    weights = model_data$weights,
    cluster = model_data$cluster,
    fixed_effects = model_data$fixed_effects,
    se_type = se_type,
    has_int = 0 %in% attr(model_data$instrument_matrix, "assign"),
    return_fit = TRUE,
    return_vcov = TRUE,
    ci = FALSE
  )
  coef_inst <- as.matrix(lm_instruments[["coefficients"]])

  if (all(colnames(model_data$instrument_matrix) %in% instruments)) {
    # if all instruments (including intercept!) are only instruments
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
      fixed_effects = model_data$fixed_effects,
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
        function(x) {
          pf(x, firststage_nomdf, dendf, lower.tail = FALSE)
        },
        numeric(1)
      ),
      gsub("value", "p.value", fstat_names)
    )
  )
}

# Wu-Hausman f-test for endogeneity
wu_hausman_reg_ftest <- function(model_data, first_stage_residuals, se_type) {

  lm_noresids <- lm_robust_fit(
    y = model_data$outcome,
    X = model_data$design_matrix,
    weights = model_data$weights,
    cluster = model_data$cluster,
    fixed_effects = model_data$fixed_effects,
    se_type = "none",
    has_int = 0 %in% attr(model_data$design_matrix, "assign"),
    ci = FALSE,
    return_fit = TRUE,
    return_vcov = FALSE
  )

  lm_resids <- lm_robust_fit(
    y = model_data$outcome,
    X = cbind(model_data$design_matrix, first_stage_residuals),
    weights = model_data$weights,
    cluster = model_data$cluster,
    fixed_effects = model_data$fixed_effects,
    se_type = se_type,
    has_int = 0 %in% attr(model_data$design_matrix, "assign"),
    ci = FALSE,
    return_fit = TRUE,
    return_vcov = TRUE
  )

  coef_noresids <- na.omit(lm_noresids[["coefficients"]])
  coef_resids <- na.omit(lm_resids[["coefficients"]])
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
    "value" = wu_hausman_fstat,
    "numdf" = wu_hausman_nomdf,
    "dendf" = dendf,
    "p.value" = pf(wu_hausman_fstat, wu_hausman_nomdf, dendf, lower.tail = FALSE)
  )
}

# Overidentification tests
# Sargan (classical ses)
sargan_chisq <- function(model_data, ss_residuals) {
  lmr <- lm_robust_fit(
    y = as.matrix(ss_residuals),
    X = model_data$instrument_matrix,
    weights = model_data$weights,
    cluster = model_data$cluster,
    fixed_effects = model_data$fixed_effects,
    se_type = "classical",
    has_int = 0 %in% attr(model_data, "assign"),
    ci = FALSE,
    return_fit = FALSE,
    return_vcov = FALSE
  )

  lmr[["r.squared"]] * lmr[["nobs"]]
}

# wooldridge score test (robust SEs)
wooldridge_score_chisq <- function(model_data, endog, instruments, ss_residuals, first_stage_fits, m) {

  # Using notation following stata ivregress postestimation help
  qhat_fit <- lm_robust_fit(
    y = model_data$instrument_matrix[, instruments[seq_len(m)], drop = FALSE],
    X = cbind(
      model_data$design_matrix[, -which(colnames(model_data$design_matrix) %in% endog), drop = FALSE],
      first_stage_fits
    ),
    weights = model_data$weights,
    cluster = model_data$cluster,
    fixed_effects = model_data$fixed_effects,
    se_type = "none",
    has_int = TRUE,
    ci = FALSE,
    return_fit = TRUE,
    return_vcov = FALSE
  )

  kmat <- as.matrix(
      model_data$instrument_matrix[, instruments[seq_len(m)], drop = FALSE] -
      qhat_fit[["fitted.values"]]
    ) *  as.vector(ss_residuals)

  if (!is.null(model_data[["weights"]])) {
    kmat_fit <- lm.fit(kmat * model_data[["weights"]], as.matrix(rep(1, times = length(ss_residuals))))
  } else {
    kmat_fit <- lm.fit(kmat, as.matrix(rep(1, times = length(ss_residuals))))
  }

  return(length(ss_residuals) - sum(residuals(kmat_fit)))
}

build_ivreg_diagnostics_mat <- function(iv_robust_out, stata = FALSE) {
  weakinst <- iv_robust_out[["diagnostic_first_stage_fstatistic"]]
  wu_hausman <- iv_robust_out[["diagnostic_endogeneity_test"]]
  overid <- iv_robust_out[["diagnostic_overid_test"]]
  n_weak_inst_fstats <- (length(weakinst) - 2) / 2

  diag_mat <- rbind(
    matrix(
      c(
        weakinst[seq_len(n_weak_inst_fstats)],
        rep(weakinst["nomdf"], n_weak_inst_fstats),
        rep(weakinst["dendf"], n_weak_inst_fstats),
        weakinst[n_weak_inst_fstats + 2 + seq_len(n_weak_inst_fstats)]
      ),
      nrow = n_weak_inst_fstats
    ),
    wu_hausman,
    c(overid[1], if (stata & overid[2] == 0) NA else overid[2], NA, overid[3])
  )[, c(2, 3, 1, 4)]

  weak_names <- "Weak instruments"
  if (n_weak_inst_fstats > 1) {
    weak_names <- paste0(
      weak_names,
      " (",
      gsub("\\:*value", "", names(weakinst[seq_len(n_weak_inst_fstats)])),
      ")"
    )
  }

  rownames(diag_mat) <- c(
    weak_names,
    "Wu-Hausman",
    "Overidentifying"
  )

  diag_mat
}
