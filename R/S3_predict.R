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

  if (!is.null(object$scaled_center)) {
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

    treat_name <- attr(object$terms, "term.labels")[1]
    interacted_covars <- X[, treat_name] * demeaned_covars

    X <- cbind(
      X[, attr(X, "assign") <= 1, drop = FALSE],
      demeaned_covars,
      interacted_covars
    )
  }

  coefs <- as.matrix(coef(object))

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

#' @export
predict.iv_robust <- function(object,
                              newdata,
                              na.action = na.pass,
                              ...) {

  X <- get_X(object, newdata, na.action)

  coefs <- as.matrix(coef(object))

  beta_na <- is.na(coefs[, 1])

  preds <- X[, !beta_na, drop = FALSE] %*% coefs[!beta_na, ]
  return(drop(preds))
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
  # The stored terms hold the two-part `y ~ x | z` formula, which
  # model.frame.default cannot parse: it reads `|` as a logical operator and
  # returns one garbage column (estimatr #397). Rebuild a one-part formula
  # naming every variable and defer to the default method.
  object <- formula
  vars <- all.vars(object[["terms"]])
  flat <- stats::reformulate(unique(vars[-1L]), response = vars[1L])
  environment(flat) <- environment(object[["terms"]])
  object[["terms"]] <- stats::terms(flat)
  call <- object[["call"]]
  call[["formula"]] <- flat
  object[["call"]] <- call
  stats::model.frame.default(object, ...)
}

#' @importFrom stats variable.names
#' @export
variable.names.lm_robust <- function(object, ...) object[["term"]]

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
#' @param x an `lm_robust` or `iv_robust` object
#' @param data the data to augment, defaulting to the model frame
#' @param newdata optional new data to predict on instead
#' @param ... ignored
#'
#' @return A `data.frame`.
#'
#' @export
augment.lm_robust <- function(x, data = NULL, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    newdata[[".fitted"]] <- predict(x, newdata = newdata)
    return(newdata)
  }

  if (is.null(data)) data <- stats::model.frame(x)

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
  data
}

#' @rdname augment.lm_robust
#' @export
augment.iv_robust <- augment.lm_robust
