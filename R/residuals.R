#' Extract model residuals for \code{\link{lm_robust}} object
#'
#' @param object an object of class "lm_robust"
#' @param type the type of residuals which should be returned, defaults to \code{"working"} and can be abbreviated
#' @param ... other arguments, unused
#'
#' @details predict.lm_robust produces predicted values, obtained by evaluating the regression function in the frame newdata for fits from \code{lm_robust} and \code{lm_lin}. If the logical se.fit is TRUE, standard errors of the predictions are calculated. Setting intervals specifies computation of confidence or prediction (tolerance) intervals at the specified level, sometimes referred to as narrow vs. wide intervals.
#'
#' The equation used for the standard error of a prediction given a row of data x is:
#'
#' \eqn{\sqrt(x \Sigma x')},
#'
#' where \eqn{\Sigma} is the estimated variance-covariance matrix from \code{lm_robust}.
#'
#' The prediction intervals are for a single observation at each case in newdata with error variance(s) pred.var. The the default is to assume that future observations have the same error variance as those used for fitting, which is gotten from the fit lm_robust object. If weights is supplied, the inverse of this is used as a scale factor. If the fit was weighted, the default is to assume constant prediction variance, with a warning.
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

  # lm_lin scaling
  if (!is.null(object$scaled_center)) {
    demeaned_covars <-
      scale(
        X[
          ,
          names(object$scaled_center),
          drop = F
          ],
        center = object$scaled_center,
        scale = FALSE
      )


    # Interacted with treatment
    treat_name <- attr(object$terms, "term.labels")[1]
    interacted_covars <- X[, treat_name] * demeaned_covars

    X <- cbind(X[, attr(X, "assign") <= 1, drop = F],
               demeaned_covars,
               interacted_covars)
  }

  # Get predicted values
  predictor <- drop(X %*% object$est)

  df_resid <- object$n - object$rank
  interval <- match.arg(interval)

  if (se.fit || interval != "none") {

    ret <- list()

    var_fit <-
      apply(X, 1, function(x) tcrossprod(crossprod(x, object$vcov), x))

    if (interval != "none") {

      tval <- qt(alpha/2, df_resid, lower.tail = FALSE)

      if (interval == "prediction") {

        # Get weights
        if (missing(weights)) {

          if (object$weighted && is.null(pred.var)) {
            warning("Assuming constant prediction variance even though model fit is weighted\\n")
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
          c(predictor,
            predictor - hwid,
            predictor + hwid),
          ncol = 3,
          dimnames = list(NULL, c('fit', 'lwr', 'upr'))
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

