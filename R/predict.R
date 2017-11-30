#' Predict method for \code{\link{lm_robust}} object
#'
#' @param object an object of class 'lm_robust'
#' @param newdata the data for which a prediction is preferred
#' @param se.fit a boolean for whether standard errors are required, default = FALSE
#' @param interval type of interval calculation. Can be abbreviated, default = none
#' @param alpha numeric denoting the test size for confidence intervals
#' @param na.action function determining what should be done with missing values in newdata. The default is to predict NA.
#' @param pred.var the variance(s) for future observations to be assumed for prediction intervals.
#' @param weights variance weights for prediction. This can be a numeric vector or a bare (unquoted) name of the weights variable in the supplied newdata.
#' @param ... other arguments, unused
#'
#' @details predict.lm_roust produces predicted values, obtained by evaluating the regression function in the frame newdata (which defaults to model.frame(object). If the logical se.fit is TRUE, standard errors of the predictions are calculated. Setting intervals specifies computation of confidence or prediction (tolerance) intervals at the specified level, sometimes referred to as narrow vs. wide intervals.
#'
#'
#' The prediction intervals are for a single observation at each case in newdata with error variance(s) pred.var. The the default is to assume that future observations have the same error variance as those used for fitting, which is gotten from the fit lm_robust object. If weights is supplied, the inverse of this is used as a scale factor. If the fit was weighted, the default is to assume constant prediction variance, with a warning.
#'
#' @export
predict.lm_robust <- function(
  object,
  newdata,
  se.fit = FALSE,
  interval = c("none", "confidence", "prediction"),
  alpha = 0.05,
  na.action = na.pass,
  pred.var = NULL,
  weights,
  ...) {

  # Get model matrix
  rhs_terms <- delete.response(object$terms)
  mf <- model.frame(rhs_terms, newdata, na.action = na.action)

  # Check class of columns in newdata match those in model fit
  if(!is.null(cl <- attr(rhs_terms, "dataClasses"))) .checkMFClasses(cl, mf)

  X <- model.matrix(rhs_terms, mf, contrasts.arg = object$contrasts)

  # Get predicted values
  predictor <- drop(X %*% object$est)

  df_resid <- object$n - object$rank
  interval <- match.arg(interval)

  if (se.fit || interval != "none") {

    ret <- list()

    var_fit <-
      apply(X, 1, function(x) tcrossprod(crossprod(x, object$vcov), x))

    if (interval != "none") {

      tval <- qt(alpha/2, df_resid)

      if (interval == "prediction") {

        if (missing(weights) &&
            object$weighted &&
            is.null(pred.var)) {
          warning("Assuming constant prediction variance even though model fit is weighted\\n")
        }

        # Get weights
        if (missing(weights)) {
          weights <- 1
        } else {
          weights <- eval(substitute(weights), newdata)
        }

        pred_var <-
          if (is.null(pred.var)) {
            object$res_var / weights
          } else {
            pred.var
          }

        hwid <- tval * sqrt(var_fit + pred_var)
      } else if (interval == "confidence") {
        hwid <- tval * sqrt(var_fit)
      }

      predictor <- cbind(predictor,
                         predictor + hwid %o% c(1, -1))
      colnames(predictor) <- c("fit", "lwr", "upr")
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

