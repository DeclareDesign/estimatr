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

  beta_na <- is.na(coefs[, 1])

  preds <- X[, !beta_na, drop = FALSE] %*% coefs[!beta_na, ]
  predictor <- drop(preds)

  df_resid <- object$df.residual
  interval <- match.arg(interval)

  if (se.fit || interval != "none") {
    if (ncol(coefs) > 1) {
      stop("Can't set `se.fit` == TRUE with multivariate outcome")
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
