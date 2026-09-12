#' Linear Regression with Lin (2013) Covariate Adjustment
#'
#' Estimates an average treatment effect with covariate adjustment following
#' Lin (2013): every covariate is centered, interacted with treatment, and
#' entered alongside it. Centering is what makes the treatment coefficient
#' the effect estimate, and the interactions avoid the bias Freedman (2008)
#' identified in ordinary covariate-adjusted regression.
#'
#' @param formula (required) An object of class formula with only the treatment on the RHS
#' @param covariates (required) A right-sided formula with pre-treatment covariates
#' @param data (optional) A `data.frame`
#' @param weights (optional) The bare (unquoted) name of the weights variable
#' @param subset (optional) A bare (unquoted) expression specifying a subset
#' @param clusters (optional) A bare (unquoted) name of the cluster variable
#' @param se_type (optional) The sort of standard error (see [lm_robust()])
#' @param ci (optional) Logical. Whether to compute p-values and confidence intervals.
#' @param alpha (optional) The significance level, 0.05 by default.
#' @param return_vcov (optional) Logical. Whether to return the vcov matrix.
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
#'   the Cholesky path carries about twice the rounding error of the QR, and
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
#' @return An object of class `"lm_robust"`, as returned by [lm_robust()],
#'   with two additions: `scaled_center`, the covariate means used for
#'   centering (taken after any function in the formula is evaluated), and
#'   `treatment_levels`. The treatment row of `coefficients` is the estimate
#'   of the average treatment effect.
#'
#' @references Lin, Winston. 2013. "Agnostic Notes on Regression Adjustments to
#'   Experimental Data: Reexamining Freedman's Critique." The Annals of Applied
#'   Statistics 7(1): 295-318. \doi{10.1214/12-AOAS583}.
#'
#' @examples
#' set.seed(20)
#' dat <- data.frame(
#'   x  = rnorm(40, mean = 2.3),
#'   x2 = rpois(40, lambda = 2),
#'   x3 = runif(40),
#'   z  = rep(0:1, 20),
#'   cl = rep(1:20, each = 2)
#' )
#' dat$y <- rnorm(40) + dat$x + 0.35 * dat$z
#'
#' # lm_robust's interface plus one argument
#' fit <- lm_lin(y ~ z, covariates = ~ x, data = dat)
#' tidy(fit)
#'
#' # Several covariates
#' lm_lin(y ~ z, covariates = ~ x + x2, data = dat)
#'
#' # Covariates are centered after any function in the formula is evaluated
#' fit2 <- lm_lin(y ~ z, covariates = ~ x + log(x3), data = dat)
#' fit2$scaled_center["log(x3)"]
#' mean(log(dat$x3))
#'
#' # Clusters, and multi-valued treatments whether or not they are factors
#' lm_lin(y ~ z, covariates = ~ x, data = dat, clusters = cl)
#' dat$z3 <- rep(1:3, length.out = 40)
#' lm_lin(y ~ z3, covariates = ~ x, data = dat)
#' lm_lin(y ~ factor(z3), covariates = ~ x, data = dat)
#'
#' # Dropping the intercept gives the mean outcome under each condition
#' lm_lin(y ~ z3 - 1, covariates = ~ x, data = dat)
#'
#' @export
lm_lin <- function(formula,
                   covariates,
                   data,
                   weights,
                   subset,
                   clusters,
                   se_type = NULL,
                   ci = TRUE,
                   alpha = .05,
                   return_vcov = TRUE,
                   try_cholesky = FALSE) {

  if (length(all.vars(rlang::f_rhs(formula))) > 1) {
    stop(
      "The `formula` argument must only have the treatment variable on the ",
      "right-hand side. Covariates should go in the `covariates` argument."
    )
  }

  if (!inherits(covariates, "formula")) {
    stop(
      "The `covariates` argument must be specified as a formula:\n",
      "You passed an object of class ", class(covariates)
    )
  }

  cov_terms <- terms(covariates)

  if (attr(cov_terms, "response") != 0) {
    stop(
      "Must not specify a response variable in `covariates` formula.\n",
      "`covariates` must be a right-sided formula, such as '~ x1 + x2 + x3'"
    )
  }

  if (length(attr(cov_terms, "order")) == 0) {
    stop(
      "`covariates` must have a variable on the right-hand side, not 0 or 1"
    )
  }

  full_formula <- update(
    formula,
    reformulate(c(".", labels(cov_terms)))
  )

  datargs <- rlang::enquos(
    formula = full_formula,
    weights = weights,
    subset = subset,
    cluster = clusters
  )
  data <- rlang::enquo(data)
  model_data <- clean_model_data(data = data, datargs)

  outcome <- as.matrix(model_data$outcome)
  n <- nrow(outcome)
  design_matrix <- model_data$design_matrix
  weights <- model_data$weights
  cluster <- model_data$cluster

  has_intercept <- attr(terms(formula), "intercept")
  treat_col <- which(attr(design_matrix, "assign") == 1)
  treatment <- design_matrix[, treat_col, drop = FALSE]
  design_mat_treatment <- colnames(design_matrix)[treat_col]

  # Kept on the fit so predict() expands a numeric multi-valued treatment
  # against the levels seen here rather than whatever happens to appear in
  # `newdata`. NULL for a binary treatment and for a factor, where the design
  # matrix already carries one column per level.
  treatment_vals <- NULL

  # Every level of the treatment, baseline included, which is what the
  # `treatment_levels` field has always carried. Distinct from
  # `treatment_vals` below, which holds only the non-baseline levels of a
  # numeric multi-valued treatment and is what predict() expands against.
  treatment_levels <- if (length(treat_col) > 1L) {
    sub(paste0("^", design_mat_treatment[1L]), "",
        design_mat_treatment)[seq_along(treat_col)]
  } else {
    sort(unique(drop(treatment)))
  }

  # The second clause is 1.0.6's and had been dropped. Without an intercept
  # there is no baseline to absorb the control group, so a 0/1 treatment has to
  # expand into both indicators: `lm_lin(y ~ z - 1)` returned `z, x_c, z:x_c`
  # and lost the control-group intercept entirely, where 1.0.6 returned
  # `z0, z1, z0:x_c, z1:x_c`.
  if (any(!(treatment %in% c(0, 1))) || (!has_intercept && ncol(treatment) == 1L)) {
    vals <- sort(unique(treatment))
    if (has_intercept) vals <- vals[-1]

    names(vals) <- paste0(colnames(design_matrix)[treat_col], vals)

    treatment <-
      outer(
        drop(treatment),
        vals,
        function(x, y) as.numeric(x == y)
      )
    treatment_vals <- vals
  }

  demeaned_covars <-
    design_matrix[
      ,
      setdiff(colnames(design_matrix), c(design_mat_treatment, "(Intercept)")),
      drop = FALSE
    ]

  if (is.numeric(weights)) {
    center <- apply(demeaned_covars, 2, weighted.mean, weights)
  } else {
    center <- colMeans(demeaned_covars)
  }

  demeaned_covars <- sweep(demeaned_covars, 2, center)

  colnames(demeaned_covars) <- lin_covar_names(colnames(demeaned_covars))

  n_treat_cols <- ncol(treatment)
  n_covars <- ncol(demeaned_covars)

  n_int_covar_cols <- n_covars * (n_treat_cols)
  interacted_covars <- matrix(0, nrow = n, ncol = n_int_covar_cols)
  interacted_covars_names <- character(n_int_covar_cols)
  for (i in 1:n_covars) {
    covar_name <- colnames(demeaned_covars)[i]

    cols <- (i - 1) * n_treat_cols + (1:n_treat_cols)
    interacted_covars[, cols] <- treatment * demeaned_covars[, i]
    interacted_covars_names[cols] <- paste0(colnames(treatment), ":", covar_name)
  }
  colnames(interacted_covars) <- interacted_covars_names

  if (has_intercept) {
    X <- cbind(
      matrix(1, nrow = n, ncol = 1, dimnames = list(NULL, "(Intercept)")),
      treatment,
      demeaned_covars,
      interacted_covars
    )
  } else {
    if (n_treat_cols == 1) {
      X <- cbind(
        treatment,
        demeaned_covars,
        interacted_covars
      )
    } else {
      X <- cbind(
        treatment,
        interacted_covars
      )
    }
  }

  return_list <-
    lm_robust_fit(
      y = outcome,
      X = X,
      weights = weights,
      cluster = cluster,
      ci = ci,
      se_type = se_type,
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      has_int = has_intercept,
      iv_stage = list(0)
    )

  return_list <- lm_return(
    return_list,
    model_data = model_data,
    formula = formula
  )

  # `center` already carries the covariates' original names, which is what
  # predict() indexes `newdata` by. (There was a setNames() call here whose
  # result was discarded, so it had never done anything.)
  return_list[["scaled_center"]] <- center
  return_list[["treatment_vals"]] <- treatment_vals
  return_list[["treatment_levels"]] <- treatment_levels

  return_list[["call"]] <- match.call()

  return(return_list)
}
