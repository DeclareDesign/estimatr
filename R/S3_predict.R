#' Predict method for `lm_robust` object
#'
#' Produces predicted values, obtained by evaluating the regression function in
#' the frame `newdata` for fits from [lm_robust()], [lm_lin()], and
#' [iv_robust()]. For an [iv_robust()] fit the regression function is the
#' second stage on the regressors themselves, as `predict()` on an `ivreg()`
#' fit computes it. If `se.fit`
#' is `TRUE`, standard errors of the predictions are calculated. Setting
#' `interval` adds confidence or prediction (tolerance) intervals at the level
#' set by `alpha`, sometimes called narrow and wide intervals respectively.
#'
#' Called without `newdata`, the method returns the in-sample fitted values,
#' and neither `se.fit` nor `interval` is available.
#'
#' The equation used for the standard error of a prediction given a row of data
#' \eqn{x} is:
#'
#' \eqn{\sqrt(x \Sigma x')},
#'
#' where \eqn{\Sigma} is the estimated variance-covariance matrix from
#' [lm_robust()].
#'
#' The prediction intervals are for a single observation at each case in
#' `newdata` with error variance(s) `pred.var`. The default is to assume that
#' future observations have the same error variance as those used for fitting,
#' which is taken from the fitted [lm_robust()] object. If `weights` is
#' supplied, the inverse of those weights scales the variance. If the fit was
#' weighted, the default is to assume constant prediction variance, with a
#' warning.
#'
#' @param object An object of class `"lm_robust"` or `"iv_robust"`.
#' @param newdata A data frame in which to look for the variables to predict
#'   from. If omitted, the fitted values are returned.
#' @param se.fit Logical. Whether to return standard errors. `FALSE` by default.
#' @param interval Type of interval calculation, which can be abbreviated.
#'   `"none"` by default.
#' @param alpha Numeric. The test size for confidence intervals.
#' @param na.action Function determining what to do with missing values in
#'   `newdata`. The default is to predict `NA`.
#' @param pred.var The variance(s) to assume for future observations when
#'   building prediction intervals.
#' @param weights Variance weights for prediction, either a numeric vector or
#'   the bare (unquoted) name of the weights variable in `newdata`.
#' @param ... (optional) Ignored.
#'
#' @return A numeric vector of predictions, or a data frame with the
#'   predictions and their standard errors and interval bounds when `se.fit` or
#'   `interval` is set.
#'
#' @seealso [lm_robust()], [lm_lin()], [iv_robust()]
#'
#' @examples
#'
#' # Set seed
#' set.seed(42)
#'
#' # Simulate data
#' n <- 10
#' dat <- data.frame(y = rnorm(n), x = rnorm(n))
#'
#' # Fit lm
#' lm_out <- lm_robust(y ~ x, data = dat)
#' # In-sample fitted values
#' predict(lm_out)
#' # Get predicted fits
#' fits <- predict(lm_out, newdata = dat)
#' # With standard errors and confidence intervals
#' fits <- predict(lm_out, newdata = dat, se.fit = TRUE, interval = "confidence")
#'
#' # Use new data as well
#' new_dat <- data.frame(x = runif(n, 5, 8))
#' predict(lm_out, newdata = new_dat)
#'
#' # You can also supply custom variance weights for prediction intervals
#' new_dat$w <- runif(n)
#' predict(lm_out, newdata = new_dat, weights = w, interval = "prediction")
#'
#' # Works for 'lm_lin' models as well
#' dat$z <- sample(1:3, size = nrow(dat), replace = TRUE)
#' lmlin_out1 <- lm_lin(y ~ z, covariates = ~ x, data = dat)
#' predict(lmlin_out1, newdata = dat, interval = "prediction")
#'
#' @importFrom rlang eval_tidy
#' @importFrom stats predict
#' @export
predict.lm_robust <- function(object,
                              newdata,
                              se.fit = FALSE,
                              interval = c("none", "confidence", "prediction"),
                              alpha = 0.05,
                              na.action = na.pass,
                              pred.var = NULL,
                              weights,
                              ...) {

  # predict() with no newdata is the in-sample fit (estimatr #403)
  if (missing(newdata)) {
    if (se.fit || match.arg(interval) != "none") {
      stop("`newdata` is required for `se.fit` or `interval`.")
    }
    return(object[["fitted.values"]])
  }

  X <- get_X(object, newdata, na.action)

  coefs <- as.matrix(coef(object))

  # An lm_lin fit is the only one carrying scaled_center.
  if (!is.null(object$scaled_center)) X <- lin_design(object, X)

  if (isTRUE(object[["fes"]])) X <- drop_absorbed_intercept(X, rownames(coefs))

  beta_na <- is.na(coefs[, 1])

  preds <- X[, !beta_na, drop = FALSE] %*% coefs[!beta_na, ]
  if (isTRUE(object[["fes"]])) preds <- add_fes(preds, object, newdata)
  predictor <- drop(preds)

  df_resid <- object$df.residual
  interval <- match.arg(interval)

  if (se.fit || interval != "none") {
    if (ncol(coefs) > 1) {
      stop("Can't set `se.fit` == TRUE with multivariate outcome")
    }
    if (isTRUE(object[["fes"]])) {
      stop("Can't set `se.fit` or `interval` with `fixed_effects`: the ",
           "absorbed group effects have no variance estimate.")
    }

    ret <- list()

    var_fit <- apply(
      X[, !beta_na, drop = FALSE],
      1,
      function(x) tcrossprod(crossprod(x, object$vcov), x)
    )

    if (interval != "none") {
      tval <- qt(alpha / 2, df_resid, lower.tail = FALSE)

      if (interval == "prediction") {

        if (missing(weights)) {
          if (object$weighted && is.null(pred.var)) {
            warning("Assuming constant prediction variance even though model fit is weighted")
          }

          weights <- 1
        } else {
          weights <- eval(substitute(weights), newdata)
        }

        if (is.null(pred.var)) {
          pred.var <- object$res_var / weights
        }

        hwid <- tval * sqrt(var_fit + pred.var)
      } else if (interval == "confidence") {
        hwid <- tval * sqrt(var_fit)
      }

      predictor <-
        matrix(
          c(
            predictor,
            predictor - hwid,
            predictor + hwid
          ),
          ncol = 3,
          dimnames = list(NULL, c("fit", "lwr", "upr"))
        )
    }

    ret[["fit"]] <- predictor

    if (se.fit) {
      ret[["se.fit"]] <- sqrt(var_fit)
    }

    return(ret)
  } else {
    return(predictor)
  }
}

# get_X() builds an iv_robust design from the regressors' terms, so lm_robust's
# method is ivreg's: without newdata the fitted values, and se.fit and the
# intervals from the fit's own vcov. The separate method this replaced had
# neither, ignoring se.fit and erroring without newdata.
#' @rdname predict.lm_robust
#' @export
predict.iv_robust <- predict.lm_robust

# Rebuild an lm_lin design exactly as lm_lin() built it: centred covariates,
# then every treatment column crossed with every covariate, covariate-major.
# `X` is the plain model matrix of treatment and covariates, as get_X() or
# model.matrix() returns it with its "assign" attribute intact.
lin_design <- function(object, X) {
  coef_names <- rownames(as.matrix(coef(object)))

  demeaned_covars <-
    scale(
      X[
        ,
        names(object$scaled_center),
        drop = FALSE
      ],
      center = object$scaled_center,
      scale = FALSE
    )
  colnames(demeaned_covars) <- lin_covar_names(names(object$scaled_center))

  # `assign == 1` is the treatment, which is one column for a binary treatment
  # and one per level for a factor. The previous code looked the treatment up
  # by term label, so a factor `z` sent it to X[, "z"] when the design matrix
  # holds `zb` and `zc` and there is no such column: subscript out of bounds
  # for every multi-valued treatment.
  treatment <- X[, attr(X, "assign") == 1, drop = FALSE]

  # A numeric treatment taking values other than 0/1 is expanded into
  # indicators against the levels seen at fit time. They have to come from the
  # fit, not from `newdata`, or predicting on a subset silently builds a
  # different design.
  if (!is.null(object$treatment_vals)) {
    treatment <- outer(
      drop(treatment),
      object$treatment_vals,
      function(x, y) as.numeric(x == y)
    )
  }

  n_covars <- ncol(demeaned_covars)
  n_treat_cols <- ncol(treatment)
  interacted_covars <- matrix(0, nrow = nrow(X), ncol = n_covars * n_treat_cols)
  interacted_names <- character(n_covars * n_treat_cols)
  for (i in seq_len(n_covars)) {
    cols <- (i - 1) * n_treat_cols + seq_len(n_treat_cols)
    interacted_covars[, cols] <- treatment * demeaned_covars[, i]
    interacted_names[cols] <-
      paste0(colnames(treatment), ":", colnames(demeaned_covars)[i])
  }
  colnames(interacted_covars) <- interacted_names

  X <- cbind(
    X[, attr(X, "assign") == 0, drop = FALSE],
    treatment,
    demeaned_covars,
    interacted_covars
  )

  # Line the design up with the coefficients by name rather than by position.
  # If the two ever disagree this errors here instead of returning a number
  # built from the wrong columns.
  missing_cols <- setdiff(coef_names, colnames(X))
  if (length(missing_cols)) {
    stop(
      "Cannot rebuild the lm_lin design from `newdata`. Missing: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  X[, coef_names, drop = FALSE]
}

# With fixed_effects the intercept is absorbed and dropped from the
# coefficients, but the terms object still carries it, so the design matrix
# rebuilt from newdata has a column the coefficients do not (estimatr #404).
drop_absorbed_intercept <- function(X, coef_names) {
  keep <- colnames(X) %in% coef_names
  X[, keep, drop = FALSE]
}

# Put the absorbed group effects back. Stored by lm_robust() for one-way FE.
add_fes <- function(preds, object, newdata) {
  fe_effects <- object[["fixed_effects"]]
  if (is.null(fe_effects)) {
    stop(
      "Can't use `predict` on this model: the absorbed group effects are ",
      "stored only for a single outcome with one set of `fixed_effects`. ",
      "With several `fixed_effects` the effects are identified in sum but ",
      "not separately, so a new observation cannot be assigned its share. ",
      "In-sample fits are in `fitted.values`."
    )
  }

  fe_formula <- eval_tidy(as.list(object[["call"]])[["fixed_effects"]])
  fe_frame <- stats::model.frame.default(fe_formula, data = newdata, na.action = NULL)
  fe_names <- paste0(names(fe_frame)[1L], as.factor(fe_frame[[1L]]))

  unknown <- setdiff(unique(fe_names), names(fe_effects))
  if (length(unknown)) {
    stop("Can't have new levels in the `fixed_effects` variable of `newdata`: ",
         paste0(unknown, collapse = ", "), ".")
  }

  # tapply returns a 1d array; as.vector keeps the matrix arithmetic conformable
  preds + as.vector(fe_effects[fe_names])
}

get_X <- function(object, newdata, na.action) {
  if (is.null(object[["terms_regressors"]])) {
    rhs_terms <- delete.response(object[["terms"]])
  } else {
    rhs_terms <- delete.response(object[["terms_regressors"]])
  }
  mf <- model.frame(
    rhs_terms,
    newdata,
    na.action = na.action,
    xlev = object[["xlevels"]]
  )

  if (!is.null(cl <- attr(rhs_terms, "dataClasses"))) .checkMFClasses(cl, mf)

  X <- model.matrix(rhs_terms, mf, contrasts.arg = object$contrasts)

  return(X)
}

#' @export
model.frame.iv_robust <- function(formula, ...) {
  # The stored formula is the two-part `y ~ x | z`, which model.frame.default
  # cannot parse: it reads `|` as a logical operator and returns one garbage
  # column (estimatr #397). clean_model_data() reads it with the Formula
  # package, as the fit did, and drops the same rows, where re-evaluating the
  # call without its clusters had kept a row dropped for a missing cluster.
  # With fixed_effects the absorbed factors join both parts, so that the frame
  # holds the variables of terms(), as lm_robust's does.
  object <- formula
  f <- object[["formula"]]
  if (isTRUE(object[["fes"]])) {
    fe_formula <- eval_tidy(as.list(object[["call"]])[["fixed_effects"]])
    fe <- paste(attr(terms(fe_formula), "term.labels"), collapse = " + ")
    parts <- vapply(1:2, function(i) {
      deparse1(formula(Formula::as.Formula(f), lhs = 0, rhs = i)[[2L]])
    }, "")
    f <- stats::as.formula(
      paste(deparse1(f[[2L]]), "~", paste(parts, fe, sep = " + ", collapse = " | ")),
      env = environment(f)
    )
  }
  rebuild_model_data(object, f, estimator = "iv")[["model_frame"]]
}

#' @importFrom stats variable.names
#' @export
variable.names.lm_robust <- function(object, ...) object[["term"]]

#' @export
model.frame.lm_robust <- function(formula, ...) {
  # The variables of terms(), so that model.frame() and terms() describe the
  # same design: the absorbed factors of a fixed_effects fit, and the centred
  # covariate expressions of an lm_lin() fit.
  object <- formula
  mf <- rebuild_model_data(object, stats::formula(terms(object)))[["model_frame"]]
  if (is.null(object[["scaled_center"]])) return(mf)
  # lm_lin()'s frame also carries the variables as the call named them, so that
  # model.matrix(formula(fit), model.frame(fit)) works, and consumers that read
  # covariates off the frame see x rather than the centred I(x - c).
  call_mf <- rebuild_model_data(object, stats::formula(object[["terms"]]))[["model_frame"]]
  stopifnot(nrow(call_mf) == nrow(mf))
  # The call's variables come first: insight names a column by the variable
  # inside it, so the first column it reads as `x` keeps that name.
  out <- cbind(call_mf, mf[setdiff(names(mf), names(call_mf))])
  attr(out, "terms") <- attr(mf, "terms")
  attr(out, "na.action") <- attr(mf, "na.action")
  out
}

# A fit stores no model frame, so it is rebuilt from the call as lm() does with
# `model = FALSE`. It goes through clean_model_data() rather than
# model.frame.default() so that the rows are the ones the fit used, including
# those dropped for a missing cluster, weight or fixed effect.
rebuild_model_data <- function(object, formula, estimator = "") {
  call <- object[["call"]]
  env <- environment(object[["terms"]])
  call_args <- c(
    weights = "weights",
    subset = "subset",
    cluster = "clusters",
    fixed_effects = "fixed_effects"
  )
  call_args <- call_args[call_args %in% names(call)]
  datargs <- c(
    list(formula = rlang::new_quosure(formula, env)),
    lapply(call_args, function(arg) rlang::new_quosure(call[[arg]], env))
  )
  data <- rlang::new_quosure(
    if (is.null(call[["data"]])) rlang::missing_arg() else call[["data"]],
    env
  )
  # The fit has already given any warning about dropped rows.
  suppressWarnings(clean_model_data(data = data, datargs, estimator = estimator))
}

#' @importFrom stats model.matrix
#' @export
model.matrix.lm_robust <- function(object, ...) {
  X <- model.matrix(
    delete.response(object[["terms"]]),
    rebuild_model_data(object, stats::formula(object[["terms"]]))[["model_frame"]],
    contrasts.arg = object[["contrasts"]]
  )
  if (!is.null(object$scaled_center)) X <- lin_design(object, X)
  if (isTRUE(object[["fes"]])) {
    X <- drop_absorbed_intercept(X, rownames(as.matrix(coef(object))))
  }
  X
}

# ivreg's model.matrix.ivreg, whose "projected" component is these same
# lm.wfit() fitted values, zero weights included. Computing them another way
# leaves an exogenous column off by rounding where it is exactly zero, and
# sandwich's meatHC() divides estfun() by model.matrix().
#' @export
model.matrix.iv_robust <- function(object,
                                   component = c("regressors", "projected", "instruments"),
                                   ...) {
  component <- match.arg(component)
  model_data <- rebuild_model_data(object, object[["formula"]], estimator = "iv")
  X <- model_data[["design_matrix"]]
  Z <- model_data[["instrument_matrix"]]
  rownames(X) <- rownames(Z) <- rownames(model_data[["model_frame"]])
  if (component == "regressors") return(X)
  if (component == "instruments") return(Z)
  w <- object[["weights"]]
  fit <- if (is.null(w)) lm.fit(Z, X) else lm.wfit(Z, X, w)
  fit[["fitted.values"]]
}

# Without these methods sandwich::bread() falls through to its default, which
# is vcov(x) * nobs(x): the fit's robust variance, not the inverse Gram matrix,
# so anything pairing it with sandwich::meat() returned a wrong number silently.
# Each is sandwich's bread for the matching model, (X'WX / n)^-1 on the
# estimated coefficients, with n the observations carrying positive weight.
#' @exportS3Method sandwich::bread
bread.lm_robust <- function(x, ...) {
  bread_from(sqrtw_design(x, "bread()"), x[["weights"]])
}

#' @exportS3Method sandwich::bread
bread.iv_robust <- function(x, ...) {
  bread_from(sqrtw_design(x, "bread()"), x[["weights"]])
}

# A difference in means and a Horvitz-Thompson estimate are not built from a
# design matrix, so there is no inverse Gram matrix to return, and sandwich's
# default handed back vcov * n as if it were one. estfun() and hatvalues()
# already error on these fits for want of a method.
refuse_bread <- function(estimator) {
  function(x, ...) {
    stop("`bread()` is not available for `", estimator, "()` fits: the ",
         "estimator is not built from a design matrix, so there is no bread ",
         "to pair with `sandwich::meat()`. Fit the regression with ",
         "`lm_robust()` for one.", call. = FALSE)
  }
}

#' @exportS3Method sandwich::bread
bread.difference_in_means <- refuse_bread("difference_in_means")

#' @exportS3Method sandwich::bread
bread.horvitz_thompson <- refuse_bread("horvitz_thompson")

# sandwich's estfun.lm and ivreg's estfun.ivreg, w * e * X, with the residuals
# on the data's scale. Without them no sandwich estimator ran on a fit.
#' @exportS3Method sandwich::estfun
estfun.lm_robust <- function(x, ...) {
  estfun_from(sqrtw_design(x, "estfun()"), x)
}

#' @exportS3Method sandwich::estfun
estfun.iv_robust <- function(x, ...) {
  estfun_from(sqrtw_design(x, "estfun()"), x)
}

# ivreg's vcovHC.ivreg: meatHC() recovers the residuals as estfun() divided by
# model.matrix(), which holds only when the latter is the projected regressors.
#' @exportS3Method sandwich::vcovHC
vcovHC.iv_robust <- function(x, ...) {
  class(x) <- c("iv_robust_projected", class(x))
  sandwich::vcovHC.default(x, ...)
}

#' @export
model.matrix.iv_robust_projected <- function(object, ...) {
  model.matrix.iv_robust(object, component = "projected")
}

# The design the coefficients were estimated on, restricted to the estimated
# ones and multiplied by the square roots of the weights. For iv_robust it is
# ivreg's projected regressors, X's projection on the instruments taken in the
# weighted metric.
sqrtw_design <- function(x, what) {
  coefs <- as.matrix(coef(x))
  if (ncol(coefs) > 1) {
    stop("`", what, "` is not available for a fit with several outcomes.",
         call. = FALSE)
  }
  if (isTRUE(x[["fes"]])) {
    stop("`", what, "` is not available with `fixed_effects`: the absorbed ",
         "design is not stored on the fit.", call. = FALSE)
  }
  keep <- !is.na(coefs[, 1])
  w <- x[["weights"]]
  sw <- if (is.null(w)) 1 else sqrt(w)
  component <- if (inherits(x, "iv_robust")) "projected" else NULL
  model.matrix(x, component = component)[, keep, drop = FALSE] * sw
}

# lm's hatvalues(), and ivreg's stage-2 ones on the projected regressors: the
# diagonal of the weighted hat matrix, for the rows with positive weight, as
# lm.influence() drops the others.
#' @importFrom stats hatvalues
#' @export
hatvalues.lm_robust <- function(model, ...) {
  X_w <- sqrtw_design(model, "hatvalues()")
  h <- rowSums((X_w %*% solve(crossprod(X_w))) * X_w)
  names(h) <- rownames(model.frame(model))
  w <- model[["weights"]]
  if (is.null(w)) h else h[w > 0]
}

#' @export
hatvalues.iv_robust <- hatvalues.lm_robust

bread_from <- function(X_w, w) {
  n <- if (is.null(w)) nrow(X_w) else sum(w > 0)
  solve(crossprod(X_w)) * n
}

# w * e * X is e * sqrt(w) * X_w, which needs no division by sqrt(w) at a zero
# weight.
estfun_from <- function(X_w, x) {
  w <- x[["weights"]]
  sw <- if (is.null(w)) 1 else sqrt(w)
  X_w * as.vector(x[["residuals"]] * sw)
}

#' @export
variable.names.iv_robust <- function(object, ...) object[["term"]]

#' @importFrom generics augment
#' @export
generics::augment

#' Augment a Model Object with Fitted Values and Residuals
#'
#' Returns the model frame with `.fitted` and `.resid` columns appended, the
#' form downstream packages expect from [broom::augment()]. Supplying
#' `newdata` returns that instead, with `.fitted` only.
#'
#' @param x An `lm_robust` or `iv_robust` object.
#' @param data The data to augment, defaulting to the model frame.
#' @param newdata Optional new data to predict on instead.
#' @param ... (optional) Ignored.
#'
#' @return A `data.frame`.
#'
#' @examples
#' set.seed(55)
#' dat <- data.frame(x = rnorm(50), z = rep(0:1, 25))
#' dat$y <- dat$x + 0.4 * dat$z + rnorm(50)
#' fit <- lm_robust(y ~ x + z, data = dat)
#'
#' head(augment(fit))
#'
#' # Supplying newdata returns predictions on it, with .fitted only
#' head(augment(fit, newdata = dat[1:5, ]))
#'
#' @export
augment.lm_robust <- function(x, data = NULL, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    newdata[[".fitted"]] <- predict(x, newdata = newdata)
    return(as_tibble(newdata))
  }

  # The variables as the call named them, not terms()'s centred expressions.
  if (is.null(data)) {
    data <- rebuild_model_data(x, stats::formula(x[["terms"]]))[["model_frame"]]
  }

  fitted <- x[["fitted.values"]]
  resid <- x[["residuals"]]
  if (is.null(fitted) || NROW(fitted) != nrow(data)) {
    stop(
      "Cannot augment this model: fitted values are not available for every ",
      "row of the data. Multivariate outcomes are not supported."
    )
  }
  if (NCOL(fitted) > 1L) {
    stop("Cannot augment a model with multiple outcomes.")
  }

  data[[".fitted"]] <- as.vector(fitted)
  if (!is.null(resid)) data[[".resid"]] <- as.vector(resid)
  as_tibble(data)
}

#' @rdname augment.lm_robust
#' @export
augment.iv_robust <- augment.lm_robust
