#' Predict method for \code{lm_robust} object
#'
#' @param object an object of class 'lm_robust'
#' @param newdata a data frame in which to look for variables with which to predict
#' @param se.fit logical. Whether standard errors are required, default = FALSE
#' @param interval type of interval calculation. Can be abbreviated, default = none
#' @param alpha numeric denoting the test size for confidence intervals
#' @param na.action function determining what should be done with missing
#' values in newdata. The default is to predict NA.
#' @param pred.var the variance(s) for future observations to be assumed for
#' prediction intervals.
#' @param weights variance weights for prediction. This can be a numeric
#' vector or a bare (unquoted) name of the weights variable in the supplied
#' newdata.
#' @param ... other arguments, unused
#'
#' @details Produces predicted values, obtained by evaluating the regression
#' function in the frame \code{newdata} for fits from \code{lm_robust} and
#' \code{lm_lin}. If the logical se.fit is TRUE, standard errors of the
#' predictions are calculated. Setting intervals specifies computation of
#' confidence or prediction (tolerance) intervals at the specified level,
#' sometimes referred to as narrow vs. wide intervals.
#'
#' The equation used for the standard error of a prediction given a row of
#' data \eqn{x} is:
#'
#' \eqn{\sqrt(x \Sigma x')},
#'
#' where \eqn{\Sigma} is the estimated variance-covariance matrix from
#' \code{lm_robust}.
#'
#' The prediction intervals are for a single observation at each case in
#' \code{newdata} with error variance(s) \code{pred.var}. The the default is to assume
#' that future observations have the same error variance as those used for
#' fitting, which is gotten from the fit \code{\link{lm_robust}} object. If
#' weights is supplied, the inverse of this is used as a scale factor. If the
#' fit was weighted, the default is to assume constant prediction variance,
#' with a warning.
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
#' # Works for 'lm_lin' models as well
#' dat$z <- sample(1:3, size = nrow(dat), replace = TRUE)
#' lmlin_out1 <- lm_lin(y ~ z, covariates = ~ x, data = dat)
#' predict(lmlin_out1, newdata = dat, interval = "prediction")
#'
#' # Predictions from Lin models are equivalent with and without an intercept
#' # and for multi-level treatments entered as numeric or factor variables
#' lmlin_out2 <- lm_lin(y ~ z - 1, covariates = ~ x, data = dat)
#' lmlin_out3 <- lm_lin(y ~ factor(z), covariates = ~ x, data = dat)
#' lmlin_out4 <- lm_lin(y ~ factor(z) - 1, covariates = ~ x, data = dat)
#'
#' predict(lmlin_out2, newdata = dat, interval = "prediction")
#' predict(lmlin_out3, newdata = dat, interval = "prediction")
#' predict(lmlin_out4, newdata = dat, interval = "prediction")
#'
#' # In Lin models, predict will stop with an error message if new
#' # treatment levels are supplied in the new data
#' new_dat$z <- sample(0:3, size = nrow(new_dat), replace = TRUE)
#' # predict(lmlin_out, newdata = new_dat)
#'
#'
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

  # Get coefs
  coefs <- as.matrix(coef(object))

  # Get prediction
  beta_na <- is.na(coefs[, 1])

  # Get predicted values
  preds <- X[, !beta_na, drop = FALSE] %*% coefs[!beta_na, ]
  if (object[["fes"]]) {
    preds <- add_fes(preds, object, newdata)
  }
  predictor <- drop(preds)

  df_resid <- object$df.residual
  interval <- match.arg(interval)

  if (se.fit || interval != "none") {
    if (ncol(coefs) > 1) {
      stop("Can't set `se.fit` == TRUE with multivariate outcome")
    }

    if (object[["fes"]]) {
      stop("Can't set `se.fit` == TRUE with `fixed_effects`")
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

  # Get predicted values
  preds <- X[, !beta_na, drop = FALSE] %*% coefs[!beta_na, ]
  if (object[["fes"]]) {
    preds <- add_fes(preds, object, newdata)
  }
  return(drop(preds))
}

get_X <- function(object, newdata, na.action) {
  # Get model matrix
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

  # Check class of columns in newdata match those in model fit
  if (!is.null(cl <- attr(rhs_terms, "dataClasses"))) .checkMFClasses(cl, mf)

  if (object[["fes"]]) {
    attr(rhs_terms, "intercept") <- 0
  }

  X <- model.matrix(rhs_terms, mf, contrasts.arg = object$contrasts)

  # lm_lin scaling (moved down from predict.lm_robust)
  if (!is.null(object$scaled_center)) {
    # Covariates
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

    # Handle treatment variable reconstruction
    treat_name <- attr(object$terms, "term.labels")[1]
    treatment <- mf[, treat_name]
    vals <- sort(unique(treatment))
    old_vals <- object$treatment_levels

    # Ensure treatment levels in newdata are subset of those for model fit
    if (!all(as.character(vals) %in% as.character(old_vals))) {
      stop(
        "Levels of treatment variable in `newdata` must be a subset of those ",
        "in the model fit."
      )
    }
    treatment <- model.matrix(~ factor(treatment, levels = old_vals) - 1)

    colnames(treatment) <- paste0(treat_name, "_", old_vals)
    # Drop out first group if there is an intercept
    if (attr(rhs_terms, "intercept") == 1) treatment <- treatment[, -1, drop = FALSE]

    # Interactions matching original fitting logic
    n_treat_cols <- ncol(treatment)
    n_covars <- ncol(demeaned_covars)

    interaction_matrix <- matrix(0, nrow = nrow(X), ncol = n_covars * n_treat_cols)

    for (i in 1:n_covars) {
      cols <- (i - 1) * n_treat_cols + (1:n_treat_cols)
      interaction_matrix[, cols] <- treatment * demeaned_covars[, i]
    }

    X <- cbind(
      if (attr(rhs_terms, "intercept") == 1) {
        matrix(1, nrow = nrow(X), ncol = 1, dimnames = list(NULL, "(Intercept)"))
      },
      treatment,
      if (attr(rhs_terms, "intercept") == 1 || ncol(treatment) == 1) demeaned_covars,
      interaction_matrix
    )
  }

  return(X)
}


add_fes <- function(preds, object, newdata) {

  # Add factors!
  args <- as.list(object[["call"]])

  if (length(all.vars(rlang::f_rhs(args[["fixed_effects"]]))) > 1) {
    stop(
      "Can't use `predict.lm_robust` with more than one set of ",
      "`fixed_effects`. Can recover fits in `fitted.values` in the model ",
      "object."
    )
  }

  # cast as factor
  femat <- model.matrix(
    ~ 0 + .,
    data = as.data.frame(
      lapply(
        stats::model.frame.default(
          args[["fixed_effects"]],
          data = newdata,
          na.action = NULL
        ),
        FUN = as.factor
      )
    )
  )

  keep_facs <- intersect(names(object[["fixed_effects"]]), colnames(femat))
  extra_facs <- setdiff(colnames(femat), names(object[["fixed_effects"]]))
  if (length(extra_facs)) {
    stop(
      "Can't have new levels in `newdata` `fixed_effects` variable, such as: ",
      paste0(extra_facs, collapse = ", ")
    )
  }
  preds <- preds +
    femat[, keep_facs] %*%
    object[["fixed_effects"]][keep_facs]

  return(preds)
}
