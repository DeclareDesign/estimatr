#' Two-Stage Least Squares Instrumental Variables Regression
#'
#' Fits a two-stage least squares instrumental variables regression and
#' returns heteroskedasticity-robust or cluster-robust standard errors, with
#' optional weak-instrument, Wu-Hausman, and overidentification diagnostics.
#'
#' @param formula (required) An object of class formula with regressors and
#'   instruments, in either of two forms. The two-part `y ~ x1 + d | x1 + z1 +
#'   z2` lists the regressors left of the bar and every instrument right of it,
#'   so an exogenous regressor appears on both sides. The three-part `y ~ x1 |
#'   d | z1 + z2`, as in `ivreg::ivreg()` and `fixest::feols()`, lists the
#'   exogenous regressors, the endogenous regressors, and the excluded
#'   instruments once each, and fits the same model; as in `ivreg::ivreg()`,
#'   `formula(fit)` returns its two-part form with the endogenous regressors
#'   first, `y ~ d + x1 | z1 + z2 + x1`, and the coefficients follow that
#'   order.
#' @param data (optional) A `data.frame`
#' @param weights (optional) The bare (unquoted) name of the weights variable
#' @param subset (optional) A bare (unquoted) expression specifying a subset
#' @param clusters (optional) A bare (unquoted) name of the cluster variable
#' @param fixed_effects (optional) A one-sided formula of fixed effects to absorb,
#'   such as `~ blockID`. Uses FWL demeaning (see [lm_robust()] for details and
#'   SE type restrictions). The diagnostics are those of the same model with
#'   the fixed effects as dummies among both the regressors and the instruments.
#' @param se_type (optional) The standard error type. `"HC2"` and `"HC3"` work with
#'   `fixed_effects` at any number of factors: the second stage runs on fitted
#'   regressors, but those are demeaned by the same fixed effects, so the
#'   leverage decomposition [lm_robust()] describes applies unchanged. `"CR2"`
#'   and `"CR3"` with `fixed_effects` expand the dummies, as CR2 did in
#'   estimatr 1.0.6.
#'   Defaults: `"HC2"` (no clusters, with or without FE), `"CR2"` (clusters, no
#'   FE), `"CR0"` (clusters, with FE). `"stata"` and `"stata_nested"` count
#'   K as [lm_robust()] describes.
#' @param ci (optional) Logical. Whether to compute p-values and confidence intervals.
#' @param alpha (optional) The significance level, 0.05 by default.
#' @param diagnostics (optional) Logical. Whether to compute IV diagnostic
#'   statistics: the first-stage F test of the excluded instruments for each
#'   endogenous regressor, a regression-based Wu-Hausman test of endogeneity,
#'   and, when the model is overidentified, a test of the overidentifying
#'   restrictions. That test is Sargan's with `se_type = "classical"` and
#'   Wooldridge's (1995) robust score test otherwise, with the score's variance
#'   summed within clusters when `clusters` is given. With `weights`, each
#'   test is the one on the model with every row multiplied by the square root
#'   of its weight. The first-stage F and Wu-Hausman tests are Wald tests
#'   under the fit's own `se_type`, so a classical one is valid exactly when
#'   the classical weighted standard errors are and a robust one exactly when
#'   the robust ones are. The robust score test uses the score's HC0 or CR0
#'   sandwich under every robust `se_type`, as Wooldridge (1995) and Stata
#'   define it; the HC1, HC2, HC3, CR2, and `"stata"` refinements correct a
#'   coefficient covariance and have no counterpart in a score test. The
#'   overidentification test is `NA`, with a warning, for a clustered fit with
#'   no more clusters than restrictions. All three reproduce Stata's
#'   `estat firststage`, `estat endogenous`, and `estat overid` on every row of
#'   the test suite's Stata fixture that Stata answers, except the robust score
#'   test after aweights under `forceweights`, where Stata computes the
#'   frequency-weight statistic instead.
#' @param return_vcov (optional) Logical. Whether to return the vcov matrix.
#' @param df_type (optional) The degrees of freedom (see [lm_robust()])
#' @param try_cholesky (optional) Logical. Whether to solve by Cholesky
#'   decomposition of `X'X` rather than by the default pivoted QR. `FALSE` by
#'   default, and worth turning on in most applied settings: about 1.4 times
#'   faster at n = 100,000 with two regressors, and 1.7 times faster at
#'   n = 200,000 with 60 regressors, where it is 0.15s against 0.25s. The
#'   saving is per fit, so it is worth most in a simulation that fits the same
#'   design thousands of times.
#'
#'   Rank deficiency is caught on either path. Redundant columns come back as
#'   `NA` exactly as they do from [lm()] whichever path ran, and a design that
#'   is rank deficient falls back to the QR.
#'
#'   Whether it is safe turns on one question, whether two regressors are
#'   nearly the same variable. Forming `X'X` squares the condition number, so
#'   the Cholesky path has about twice the rounding error of the QR, and
#'   only near-collinearity makes that visible. Differences of scale do not,
#'   because the columns are normalized before either decomposition, so a
#'   covariate in dollars beside one in years costs nothing. For a treatment
#'   indicator, a few covariates, block or cluster dummies, the centered
#'   interactions [lm_lin()] builds, or a factorial, the two paths agree to at
#'   least 10 significant digits, which is why [difference_in_means()] sets it
#'   to `TRUE` internally. Agreement falls to about 3 digits as the scaled
#'   condition index reaches `1e6`, and the QR fallback takes over above
#'   roughly `1e8`. Nothing interpretable lives in that range: a design at
#'   `1e6` returns a coefficient of 4.8e4 with a standard error of 4.6e4 on a
#'   regressor whose true effect is zero. To check a design directly, scale the
#'   columns first, since the unscaled condition number of a design in mixed
#'   units is large for a reason that does not affect the fit:
#'   `kappa(sweep(X, 2, sqrt(colSums(X^2)), "/"), exact = TRUE)`.
#'
#' @return An object of class `"iv_robust"`, a list holding the estimate table in `coefficients`, `std.error`, `df`, `statistic`,
#'   `p.value`, `conf.low`, `conf.high`, `term`, and `outcome`; the fit in
#'   `fitted.values`, `residuals`, `vcov`, `nobs`, `k`, `rank`, `df.residual`,
#'   and `res_var`; the summary statistics `r.squared`, `adj.r.squared`,
#'   `tss`, and `fstatistic`; and `se_type`, `weighted`, `clustered`, `fes`,
#'   `alpha`, `terms`, `xlevels`, and `call`.
#'
#'   `residuals` are the structural residuals, `y - X beta`, rather than the
#'   second-stage ones. `ei.iv`, `terms_regressors`, and `formula` record the
#'   two-stage structure. `first_stage` is a list of [lm_robust()] fits, one
#'   per endogenous regressor and named by it, each the regression of that
#'   regressor on every instrument with the fit's `se_type`, clusters, weights,
#'   fixed effects, and `df_type`, so `tidy()`, `summary()`, and table packages
#'   take them as they are (see the examples). Each equals the [lm_robust()]
#'   call recorded in its `call`, which is `NULL` for the dummy column of an
#'   endogenous factor. With `diagnostics = TRUE` the object also holds
#'   `diagnostic_first_stage_fstatistic`, `diagnostic_endogeneity_test`, and
#'   `diagnostic_overid_test`.
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
#' # Exogenous regressors, endogenous regressors, and instruments, each once
#' dat$w <- rnorm(n)
#' iv_robust(y ~ w | x | z, data = dat)
#'
#' # The same variance menu as lm_robust()
#' iv_robust(y ~ x | z, data = dat, se_type = "classical")
#' iv_robust(y ~ x | z, data = dat, clusters = cl)
#'
#' # The first stage, one lm_robust fit per endogenous regressor
#' tidy(fit$first_stage$x)
#'
#' # Weak-instrument, endogeneity, and overidentification tests
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
                      try_cholesky = FALSE,
                      df_type = NULL) {
  datargs <- rlang::enquos(
    formula = formula,
    weights = weights,
    subset = subset,
    cluster = clusters,
    fixed_effects = fixed_effects
  )
  data <- rlang::enquo(data)
  model_data <- clean_model_data(data = data, datargs, estimator = "iv")

  has_fe  <- !is.null(model_data[["fixed_effects"]])
  fe_rank <- 0L
  fe_rank_k <- NULL
  fe_lev <- NULL
  yoriginal <- NULL

  # The Wu-Hausman and over-identification tests are tests on one outcome's
  # structural residuals, and the diagnostic fields hold one statistic each.
  # With several outcomes the score test used to die on a dims mismatch
  # (1.0.6 too); skip with the same notice fixed effects give.
  if (diagnostics && ncol(as.matrix(model_data[["outcome"]])) > 1L) {
    warning("Diagnostics are not available with multiple outcomes. Skipping.")
    diagnostics <- FALSE
  }

  if (has_fe) {
    yoriginal  <- as.matrix(model_data[["outcome"]])
    zoriginal <- model_data[["instrument_matrix"]]
    zoriginal <- zoriginal[, colnames(zoriginal) != "(Intercept)", drop = FALSE]
    model_data <- demean_fes(model_data)
    model_data[["instrument_matrix"]] <- demean_matrix_by_fes(
      model_data[["instrument_matrix"]], model_data
    )
    fe_rank <- sum(model_data[["fe_levels"]]) - length(model_data[["fe_levels"]]) + 1L

    # The nominal count above overstates the rank whenever one factor is partly
    # spanned by the others -- a nested factor, or a disconnected design --
    # which inflates the rank correction and shrinks the residual degrees of
    # freedom. fe_leverage() returns the exact rank from the same
    # eigendecomposition that gives the leverage, so it is used for every
    # se_type; only HC2 and HC3 also need the vector.
    fe_proj <- fe_leverage(model_data[["fe_codes"]], model_data[["weights"]],
                           leverage = needs_fe_leverage(se_type, !is.null(model_data[["cluster"]])))
    fe_rank <- fe_proj[["rank"]]
    fe_lev <- fe_proj[["leverage"]]
    if (identical(se_type, "stata_nested") && !is.null(model_data[["cluster"]])) {
      fe_rank_k <- fe_rank_unnested(model_data[["fe_codes"]], model_data[["cluster"]],
                                    model_data[["weights"]], fe_rank)
    }

    # There is nothing to instrument once the fixed effects have absorbed the
    # whole design. 1.0.6 reached the solver and died on "length of 'dimnames'
    # [2] not equal to array extent"; lm_robust() answers this case, so say
    # which function to call rather than repeating the cryptic failure.
    if (ncol(model_data[["design_matrix"]]) == 0L) {
      stop(
        "The model has no regressors left once `fixed_effects` are absorbed, ",
        "so there is nothing to instrument.\nUse `lm_robust()` with the same ",
        "`fixed_effects` if the absorbed fit is what you want."
      )
    }
  }

  # The absorbed fixed effects as every fit on the demeaned data needs them:
  # the second stage, the first-stage fits, and the regressions behind the
  # diagnostics. By FWL each of those is the fit with the dummies, given the
  # dummies' rank, leverage, and, for CR2 and CR3, the dummies themselves.
  fe_args <- list(
    fe_rank = fe_rank,
    fe_leverage = fe_lev,
    femat = if (has_fe && needs_fe_dummies(se_type, df_type, !is.null(model_data[["cluster"]])))
      fe_dummy_matrix(model_data)
      else NULL,
    fe_rank_k = if (is.null(fe_rank_k)) fe_rank else fe_rank_k
  )

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

  fitted_values <- zero_underidentified(model_data, first_stage)

  second_stage <- do.call(lm_robust_fit, c(list(
    y = model_data$outcome,
    X = fitted_values,
    weights = model_data$weights,
    cluster = model_data$cluster,
    ci = ci,
    se_type = se_type,
    has_int = attr(model_data$terms, "intercept"),
    alpha = alpha,
    return_vcov = return_vcov,
    try_cholesky = try_cholesky,
    iv_stage = list(2, model_data$design_matrix),
    df_type = df_type
  ), fe_args))

  return_list <- lm_return(
    second_stage,
    model_data = model_data,
    formula = model_data$formula
  )

  se_type <- return_list[["se_type"]]

  # ------
  # First-stage fits
  # ------
  # One lm_robust fit per endogenous regressor, on every instrument, with the
  # second stage's variance, clusters, weights, and absorbed fixed effects, so
  # that each is the fit lm_robust() returns for `x ~ <instruments>` called
  # with the same arguments (#335).
  roles <- instrument_roles(model_data, first_stage)
  iv_call <- match.call()
  instruments_rhs <- formula(model_data[["formula"]], lhs = 0, rhs = 2)[[2]]
  first_stage_data <- model_data
  first_stage_data[["design_matrix"]] <- model_data[["instrument_matrix"]]
  first_stage_data[["terms"]] <- terms(model_data[["formula"]], rhs = 2)
  if (has_fe) first_stage_data[["Xoriginal"]] <- zoriginal
  return_list[["first_stage"]] <- setNames(lapply(roles[["endog"]], function(endog) {
    first_stage_data[["outcome"]] <- model_data[["design_matrix"]][, endog, drop = FALSE]
    fit <- do.call(lm_robust_fit, c(list(
      y = first_stage_data[["outcome"]],
      X = first_stage_data[["design_matrix"]],
      weights = model_data[["weights"]],
      cluster = model_data[["cluster"]],
      ci = ci,
      se_type = se_type,
      has_int = attr(first_stage_data[["terms"]], "intercept"),
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      df_type = df_type
    ), fe_args))
    fit <- lm_return(fit, model_data = first_stage_data, formula = NULL)
    if (has_fe) {
      fit <- fe_full_scale(fit, model_data[["Xoriginal"]][, endog], first_stage_data)
    }
    # The call lm_robust() would be given for this fit, where the endogenous
    # regressor is a term of its own. A factor's dummy column is not, and
    # there is no such call to record.
    if (endog %in% attr(model_data[["terms_regressors"]], "term.labels")) {
      fit_call <- iv_call[names(iv_call) != "diagnostics"]
      fit_call[[1L]] <- quote(lm_robust)
      fit_call[["formula"]] <- call("~", str2lang(endog), instruments_rhs)
      fit[["call"]] <- fit_call
    }
    class(fit) <- "lm_robust"
    fit
  }), roles[["endog"]])

  # ------
  # diagnostics
  # ------
  if (diagnostics) {

    endog <- roles[["endog"]]

    first_stage_fits <- first_stage[["fitted.values"]][, endog, drop = FALSE]
    colnames(first_stage_fits) <- paste0("fit_", colnames(first_stage_fits))

    # Only the endogenous regressors. An exogenous regressor's first-stage
    # residual is zero up to rounding, and carrying it left the Wu-Hausman
    # auxiliary regression relying on rank detection to drop a column of noise.
    # Both stats::lm() and a pivoted QR at tol = 1e-7 keep such a column, which
    # adds a spurious degree of freedom to the test's numerator.
    first_stage_residuals <-
      (model_data$design_matrix - first_stage[["fitted.values"]])[, endog, drop = FALSE]
    colnames(first_stage_residuals) <- paste0("resid_", endog)

    extra_instruments <- ncol(roles[["excluded"]]) - length(endog)

    # With weights, every diagnostic is the textbook test on the model with each
    # row multiplied by sqrt(w). The first-stage F and Wu-Hausman tests are
    # Wald tests under the fit's own se_type, so a classical one is valid
    # exactly when the classical weighted standard errors are, which is when
    # the weights are inverse error variances, and a robust one exactly when
    # the robust standard errors are. The robust score test's S is the score's
    # HC0 or CR0 sandwich under every robust se_type, as Wooldridge and Stata
    # define it: the leverage and small-sample refinements correct a
    # coefficient covariance and have no counterpart in a score test.
    wu_hausman_ftest_val <- wu_hausman_reg_ftest(model_data, first_stage_residuals, se_type, fe_args)

    if (extra_instruments) {
      ss_residuals <- model_data$outcome - second_stage[["fitted.values"]]

      if (se_type == "classical") {
        overid_chisq_val <- sargan_chisq(model_data, ss_residuals)
      } else {
        overid_chisq_val <- wooldridge_score_chisq(
          model_data = model_data,
          exogenous = roles[["exogenous"]],
          excess = roles[["excluded"]][, seq_len(extra_instruments), drop = FALSE],
          ss_residuals = ss_residuals,
          first_stage_fits = first_stage_fits
        )
      }

      overid_chisqtest_val <- c(
        overid_chisq_val,
        extra_instruments,
        pchisq(overid_chisq_val, extra_instruments, lower.tail = FALSE)
      )

    } else {
      overid_chisqtest_val <- c(NA_real_, extra_instruments, NA_real_)
    }
    names(overid_chisqtest_val) <- c("value", "df", "p.value")

    first_stage_ftest_val <- first_stage_ftest(model_data, roles, se_type, fe_args)

    return_list[["diagnostic_first_stage_fstatistic"]] <- first_stage_ftest_val
    return_list[["diagnostic_endogeneity_test"]] <- wu_hausman_ftest_val
    return_list[["diagnostic_overid_test"]] <- overid_chisqtest_val
  }
  if (has_fe) {
    return_list <- fe_full_scale(return_list, yoriginal, model_data)
  }

  return_list[["call"]] <- match.call()

  return_list[["terms_regressors"]] <- model_data[["terms_regressors"]]
  return_list[["formula"]] <- iv_two_part_formula(formula(formula))
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

# Which regressors are endogenous, and a basis for the instruments they
# exclude, decided by the column space rather than by name.
#
# A regressor is exogenous when the instruments reproduce it: its first-stage
# residual has norm at most 1e-7 of its own, the per-column test rank detection
# uses. Deciding by name made every diagnostic depend on how the instrument
# formula was spelled. With x1 in the span of the instruments but not listed
# among them, x1 counted as endogenous, the Wu-Hausman auxiliary regression
# carried a column of rounding error, and the test returned 0.144 on 2 degrees
# of freedom where the same model with x1 listed returns 0.109 on 1, as
# AER::ivreg and estimatr 1.0.6 both do for either spelling.
#
# The excluded instruments are the instrument columns not named among the
# regressors, residualized on the exogenous regressors and reduced to a basis
# by a pivoted QR at the same tolerance. The weak-instrument F and the
# overidentification tests are invariant to which basis of that space is used,
# so where the names already line up no number moves beyond rounding.
#
# Both rank decisions are made in the metric the first stage was fit in, rows
# scaled by sqrt(w). The fitted values are on the scale of the data, so where
# zero weights make a dummy collinear (#305) the fit reproduces it on every
# row that counts and misses it on the rows that do not. Measured unweighted,
# that miss called an exogenous dummy endogenous, gave it a first-stage F of
# its own, and left the overidentification test on -1 degrees of freedom.
instrument_roles <- function(model_data, first_stage) {
  X <- model_data$design_matrix
  Z <- model_data$instrument_matrix
  sw <- if (is.null(model_data$weights)) 1 else sqrt(model_data$weights)

  resid_norm <- sqrt(colSums((sw * (X - first_stage[["fitted.values"]]))^2))
  reproduced <- resid_norm <= 1e-7 * sqrt(colSums((sw * X)^2))
  exogenous <- colnames(X)[reproduced]

  candidates <- Z[, !(colnames(Z) %in% colnames(X)), drop = FALSE]
  weighted_candidates <- sw * candidates
  if (length(exogenous)) {
    candidates <- qr.resid(qr(X[, exogenous, drop = FALSE]), candidates)
    weighted_candidates <-
      qr.resid(qr(sw * X[, exogenous, drop = FALSE]), weighted_candidates)
  }
  basis <- qr(weighted_candidates, tol = 1e-7)

  list(
    endog = colnames(X)[!reproduced],
    exogenous = exogenous,
    excluded = candidates[, basis$pivot[seq_len(basis$rank)], drop = FALSE]
  )
}


# An underidentified model has no 2SLS estimate for the regressors the
# instruments cannot reproduce, and fitting one anyway returned a
# clean-looking object. The second stage regresses on the first-stage fitted
# values, which span no more than the instruments do, so rank detection
# dropped whichever column lm_solver()'s pivot ranked last. On
# `mpg ~ hp + cyl | am` that was the intercept, and hp and cyl came back with
# standard errors that nothing identifies; on `y ~ x | z` with a flat first
# stage the intercept went too, and `x` came back carrying the intercept's
# value under the endogenous regressor's name, which is the number a reader
# takes for the LATE.
#
# Which column goes is therefore an answer here rather than a naming question,
# and it is not left to the pivot. The exogenous regressors and the intercept
# are kept, since the first stage reproduces them exactly, and the rest of the
# drop set comes from a dqrdc2 QR of the fitted values with those columns
# moved to the front. dqrdc2 keeps the earliest columns it can, so the
# endogenous regressors go in reverse formula order, as stats::lm() and
# AER::ivreg() drop them. The reordering is what does the work: `y ~ x + w |
# z + w` puts the endogenous regressor ahead of the exogenous one, which is
# what a drop rule reading position will honour.
#
# The reordering is not optional, even though lm_solver() is itself order
# preserving. lm_solver() reads the design in formula order, and on
# `y ~ x + w | w` the fitted values of x and w are collinear with w written
# second, so a position rule on the unordered design keeps x and drops w, the
# exogenous regressor the first stage reproduces exactly. Moving the exogenous
# columns to the front is what makes x, the regressor the instruments cannot
# reproduce, the column that goes.
#
# The dropped columns are zeroed rather than removed. A zero column has zero
# residual norm, so lm_solver() drops it wherever it sits and returns its
# coefficient as NA in place: the fit on the rest is unchanged,
# drop_collinear() takes the same columns off the design matrix, and
# lm_return() reports the NAs in the words it already uses for a collinear
# regressor. Only endogenous columns are ever zeroed, so collinearity among
# the regressors themselves is still resolved by lm_robust()'s own rule.
zero_underidentified <- function(model_data, first_stage) {
  fitted_values <- first_stage[["fitted.values"]]
  if (qr(fitted_values)$rank >= qr(model_data$design_matrix)$rank) {
    return(fitted_values)
  }

  roles <- instrument_roles(model_data, first_stage)
  ordered <- c(roles[["exogenous"]], roles[["endog"]])
  ordered_qr <- qr(fitted_values[, ordered, drop = FALSE])
  kept <- ordered[ordered_qr$pivot[seq_len(ordered_qr$rank)]]
  fitted_values[, setdiff(roles[["endog"]], kept)] <- 0

  fitted_values
}

first_stage_ftest <- function(model_data, roles, se_type, fe_args = list()) {

  endog <- roles[["endog"]]
  n_exogenous <- length(roles[["exogenous"]])

  lm_instruments <- do.call(lm_robust_fit, c(list(
    y = model_data$design_matrix[, endog, drop = FALSE],
    X = cbind(
      model_data$design_matrix[, roles[["exogenous"]], drop = FALSE],
      roles[["excluded"]]
    ),
    weights = model_data$weights,
    cluster = model_data$cluster,
    se_type = se_type,
    has_int = 0 %in% attr(model_data$instrument_matrix, "assign"),
    return_fit = TRUE,
    return_vcov = TRUE,
    ci = FALSE
  ), fe_args))
  coef_inst <- as.matrix(lm_instruments[["coefficients"]])

  if (n_exogenous == 0L) {
    firststage_nomdf <- lm_instruments[["rank"]]
    firststage_fstat_value <- lm_instruments[["fstatistic"]][seq_len(length(endog))]
  } else {
    # The coefficients are padded with NA for any column the first stage
    # dropped (a dummy that weights have made collinear, say), but the
    # variance matrix is built on the columns kept. Positions in the one are
    # not positions in the other, so test only the excluded instruments that
    # survived and index the variance by their place among the survivors.
    excluded_idx <- n_exogenous + seq_len(ncol(roles[["excluded"]]))
    kept <- which(!is.na(coef_inst[, 1L]))
    tested <- intersect(excluded_idx, kept)
    firststage_nomdf <- length(tested)
    firststage_fstat_value <- if (firststage_nomdf == 0L) {
      rep(NA_real_, ncol(coef_inst))
    } else {
      compute_fstat(
        coef_matrix = coef_inst[kept, , drop = FALSE],
        coef_indices = match(tested, kept),
        vcov_fit = lm_instruments[["vcov"]],
        rank = lm_instruments[["rank"]],
        nomdf = firststage_nomdf
      )
    }
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

wu_hausman_reg_ftest <- function(model_data, first_stage_residuals, se_type, fe_args = list()) {

  has_int <- 0 %in% attr(model_data$design_matrix, "assign")

  lm_noresids <- do.call(lm_robust_fit, c(list(
    y = model_data$outcome,
    X = model_data$design_matrix,
    weights = model_data$weights,
    cluster = model_data$cluster,
    se_type = "none",
    has_int = has_int,
    ci = FALSE,
    return_fit = TRUE,
    return_vcov = FALSE
  ), fe_args))

  lm_resids <- do.call(lm_robust_fit, c(list(
    y = model_data$outcome,
    X = cbind(model_data$design_matrix, first_stage_residuals),
    weights = model_data$weights,
    cluster = model_data$cluster,
    se_type = se_type,
    has_int = has_int,
    ci = FALSE,
    return_fit = TRUE,
    return_vcov = TRUE
  ), fe_args))

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
  weights <- model_data[["weights"]]
  ss_resid_lm <- lm_robust_fit(
    y = ss_residuals,
    X = model_data$instrument_matrix,
    weights = weights,
    cluster = NULL,
    se_type = "classical",
    has_int = attr(model_data$terms, "intercept"),
    return_fit = FALSE,
    return_vcov = FALSE,
    ci = FALSE
  )

  # A zero-weight row is not an observation, as everywhere else in the package.
  n_obs <- if (is.null(weights)) nrow(model_data$instrument_matrix) else sum(weights > 0)
  n_obs * ss_resid_lm[["r.squared"]]
}

# Wooldridge's robust score test for over-identifying restrictions. The
# rewrite had once replaced the construction with `n * R^2` of the second-stage
# residuals on [Z, xhat], which is algebraically the Sargan statistic rather
# than the robust score test, and indexed `first_stage_fits` by the endogenous
# variable names after renaming its columns `fit_<name>`, so every
# over-identified fit with a non-classical `se_type` errored with "subscript out
# of bounds" instead of returning either number.
#
# `excess` is any m columns of a basis for the excluded instruments. Wooldridge
# shows the statistic does not depend on which m are chosen.
#
# The score contributions k are the excess instruments residualized on
# [exogenous, xhat], times the structural residual. The statistic is s' S^-1 s,
# with s the sum of k and S its outer product summed within clusters, or over
# observations when there are none. Unclustered, that is n minus the residual
# sum of squares from regressing a vector of ones on k, the form Wooldridge
# (1995) gives and Stata's `estat overid` reports. estimatr 1.0.6 used that form
# for clustered fits as well, which ignores the clustering.
wooldridge_score_chisq <- function(model_data,
                                    exogenous,
                                    excess,
                                    ss_residuals,
                                    first_stage_fits) {

  qhat_fit <- lm_robust_fit(
    y = excess,
    X = cbind(
      model_data$design_matrix[, exogenous, drop = FALSE],
      first_stage_fits
    ),
    weights = model_data$weights,
    cluster = model_data$cluster,
    se_type = "none",
    has_int = TRUE,
    ci = FALSE,
    return_fit = TRUE,
    return_vcov = FALSE
  )

  kmat <- as.matrix(excess - qhat_fit[["fitted.values"]]) * as.vector(ss_residuals)
  # On the weighted fit's transformed data each contribution is sqrt(w) r times
  # sqrt(w) u, so s sums w k and S sums w^2 k k', the weighting of the HC meat.
  if (!is.null(model_data[["weights"]])) {
    kmat <- kmat * model_data[["weights"]]
  }
  cluster <- model_data[["cluster"]]
  kmat_sums <- if (is.null(cluster)) kmat else rowsum(kmat, cluster)
  meat_qr <- qr(crossprod(kmat_sums))

  # With no more clusters than restrictions S is singular, and s' S^+ s would
  # equal the number of clusters whatever the data. A cluster whose weights are
  # all zero is not a cluster: it contributes a zero row to `kmat_sums`, and
  # counting rows let three clusters with one weighted out pass this guard and
  # return exactly 2 on 2 df, with no warning, whatever the data.
  weights <- model_data[["weights"]]
  observed <- if (is.null(weights)) rep(TRUE, nrow(kmat)) else weights > 0
  n_groups <- if (is.null(cluster)) sum(observed) else length(unique(cluster[observed]))
  if (n_groups <= ncol(kmat) || meat_qr$rank < ncol(kmat)) {
    warning(
      "`diagnostic_overid_test` is NA: the variance of the score for the ",
      "overidentifying restrictions is singular, as it is whenever there are ",
      "no more clusters than overidentifying restrictions."
    )
    return(NA_real_)
  }

  score <- colSums(kmat)
  drop(crossprod(score, qr.coef(meat_qr, score)))
}
