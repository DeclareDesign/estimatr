# An lh_robust() fit is a container: the lm_robust() fit under `lm_robust`,
# the hypothesis table under `lh`, and the joint test. The fit-level generics
# are the lm_robust() part's, so that coef() and vcov() name the same
# parameters and consumers that pair them work. The hypothesis rows are
# reported beside the coefficients by tidy(), confint(), summary(), and print(),
# and on their own by coef() and vcov() of `fit$lh`. Without these methods
# coef() and df.residual() returned NULL, vcov() and predict() had no method,
# and model.frame(), terms(), and formula() failed on the container's fields.

#' @importFrom stats coef
#' @export
coef.lh_robust <- function(object, ...) coef(object[["lm_robust"]], ...)

#' @importFrom stats vcov
#' @export
vcov.lh_robust <- function(object, ...) vcov(object[["lm_robust"]], ...)

#' @export
df.residual.lh_robust <- function(object, ...) {
  df.residual(object[["lm_robust"]], ...)
}

#' @export
residuals.lh_robust <- function(object, ...) {
  residuals(object[["lm_robust"]], ...)
}

#' @export
fitted.lh_robust <- function(object, ...) fitted(object[["lm_robust"]], ...)

#' @export
predict.lh_robust <- function(object, ...) predict(object[["lm_robust"]], ...)

#' @export
model.frame.lh_robust <- function(formula, ...) {
  model.frame(formula[["lm_robust"]], ...)
}

#' @export
model.matrix.lh_robust <- function(object, ...) {
  model.matrix(object[["lm_robust"]], ...)
}

#' @export
terms.lh_robust <- function(x, ...) terms(x[["lm_robust"]], ...)

#' @importFrom stats formula
#' @export
formula.lh_robust <- function(x, ...) formula(x[["lm_robust"]], ...)

#' @importFrom stats weights
#' @export
weights.lh_robust <- function(object, ...) weights(object[["lm_robust"]], ...)

#' @export
sigma.lh_robust <- function(object, ...) sigma(object[["lm_robust"]], ...)

#' @export
deviance.lh_robust <- function(object, ...) deviance(object[["lm_robust"]], ...)

#' @importFrom stats hatvalues
#' @export
hatvalues.lh_robust <- function(model, ...) hatvalues(model[["lm_robust"]], ...)

#' @export
case.names.lh_robust <- function(object, ...) {
  case.names(object[["lm_robust"]], ...)
}

#' @export
variable.names.lh_robust <- function(object, ...) {
  variable.names(object[["lm_robust"]], ...)
}

#' @export
labels.lh_robust <- function(object, ...) labels(object[["lm_robust"]], ...)

#' @export
augment.lh_robust <- function(x, ...) augment(x[["lm_robust"]], ...)

#' @exportS3Method sandwich::bread
bread.lh_robust <- function(x, ...) sandwich::bread(x[["lm_robust"]], ...)

#' @exportS3Method sandwich::estfun
estfun.lh_robust <- function(x, ...) sandwich::estfun(x[["lm_robust"]], ...)

#' @exportS3Method sandwich::vcovHC
vcovHC.lh_robust <- function(x, ...) sandwich::vcovHC(x[["lm_robust"]], ...)

# summary() on an lh_robust() fit is the container, printed one part at a
# time. Its coefficient table is the fit's with the hypothesis rows appended,
# as tidy() reports them; the column names are the fit's, so a one-sided
# hypothesis's p-value sits in the same column as the fit's two-sided ones,
# as it does in tidy()'s p.value column.
#' @export
coef.summary.lh_robust <- function(object, ...) {
  rbind(coef(summary(object[["lm_robust"]])),
        unclass(summary(object[["lh"]])))
}

# The hypothesis table on its own: one estimate per hypothesis, named by the
# hypothesis as written, and the variance of those linear combinations, which
# car::linearHypothesis() computed from the fit's vcov().
#' @export
coef.lh <- function(object, ...) {
  setNames(object[["coefficients"]], object[["term"]])
}

#' @export
vcov.lh <- function(object, ...) {
  attr(attr(object, "linear_hypothesis"), "vcov")
}
