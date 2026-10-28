# terms() describes the whole design the fit used, so that a consumer pairing
# model.matrix(terms(fit), model.frame(fit)) with the fit's residuals gets the
# fit's regressors. formula() stays the call's formula, which is what update()
# and print need. The terms the fit itself predicts and summarises with are
# kept in `object$terms`, and estimatr's own methods read that field directly.

# With fixed_effects the absorbed factors are part of the design but not of the
# coefficients, so model.matrix(terms(fit)) has more columns than coef(fit) and
# cannot be paired with them by mistake.
#' @importFrom stats terms
#' @export
terms.lm_robust <- function(x, ...) {
  trms <- x[["terms"]]
  if (!isTRUE(x[["fes"]])) return(trms)
  with_absorbed_factors(x, trms, attr(trms, "dataClasses"))
}

# iv_robust's stored terms are the two-part formula, instruments included. Its
# design is the regressors, as ivreg's terms() default to.
#' @export
terms.iv_robust <- function(x, ...) {
  trms <- x[["terms_regressors"]]
  classes <- attr(x[["terms"]], "dataClasses")
  if (!isTRUE(x[["fes"]])) return(with_data_classes(trms, classes))
  with_absorbed_factors(x, trms, classes)
}

with_absorbed_factors <- function(x, trms, classes) {
  fe_formula <- eval_tidy(as.list(x[["call"]])[["fixed_effects"]])
  full <- reformulate(
    c(attr(trms, "term.labels"), attr(terms(fe_formula), "term.labels")),
    response = if (attr(trms, "response")) attr(trms, "variables")[[2L]],
    intercept = attr(trms, "intercept"),
    env = environment(trms)
  )
  fe_vars <- vapply(attr(terms(fe_formula), "variables")[-1L], deparse1, "")
  # The absorbed factors are coerced to factors whatever their class in data.
  with_data_classes(
    terms(full),
    c(classes, setNames(rep("factor", length(fe_vars)), fe_vars))
  )
}

# Consumers such as margins read the variable classes off terms(fit), as they
# can off an lm fit's terms.
with_data_classes <- function(trms, classes) {
  vars <- vapply(attr(trms, "variables")[-1L], deparse1, "")
  # Keyed as terms() deparses them, which rounds the %.17g constants.
  names(classes) <- vapply(names(classes), function(v) deparse1(str2lang(v)), "")
  attr(trms, "dataClasses") <- classes[vars]
  trms
}

# lm_lin's coefficients are on the centred covariates, which are not data
# columns. Its terms write each centred design column out with its centre,
# I(x - c), and interact every treatment column with each, in the order
# lm_lin() built them, so model.matrix() on these terms is the fit's design
# column for column. Only the coefficient names differ.
#' @export
terms.lm_lin <- function(x, ...) {
  trms <- x[["terms"]]
  labels <- attr(trms, "term.labels")
  treat_label <- labels[1L]

  treat_terms <- if (is.null(x[["treatment_vals"]])) {
    treat_label
  } else {
    # Numeric, since a logical term is coded as a factor.
    sprintf(
      "I(as.numeric(%s == %s))",
      treat_label,
      format_constant(x[["treatment_vals"]])
    )
  }

  covar_exprs <- lin_covariate_exprs(x, trms, labels[-1L])
  centre <- x[["scaled_center"]]
  covar_terms <- sprintf(
    "I(%s %s %s)",
    covar_exprs,
    ifelse(centre < 0, "+", "-"),
    format_constant(abs(centre))
  )

  has_intercept <- attr(trms, "intercept") == 1L
  # Without an intercept and with more than one treatment column, lm_lin()
  # leaves out the main covariate effects; whether it did is read off the
  # coefficient names, which carry a dropped column as NA.
  coef_names <- rownames(as.matrix(coef(x)))
  main_covars <- if (all(lin_covar_names(names(centre)) %in% coef_names)) covar_terms

  rhs <- c(
    if (!has_intercept) "0",
    treat_terms,
    main_covars,
    as.vector(outer(treat_terms, covar_terms, paste, sep = ":"))
  )
  lhs <- if (attr(trms, "response")) deparse1(attr(trms, "variables")[[2L]])
  f <- stats::as.formula(
    paste(lhs, "~", paste(rhs, collapse = " + ")),
    env = environment(trms)
  )
  # Every written-out treatment and covariate column is numeric.
  classes <- attr(trms, "dataClasses")
  with_data_classes(
    terms(f, keep.order = TRUE),
    c(
      classes[c(lhs, if (is.null(x[["treatment_vals"]])) treat_label)],
      setNames(rep("numeric", length(treat_terms) + length(covar_terms)),
               c(if (!is.null(x[["treatment_vals"]])) treat_terms, covar_terms))
    )
  )
}

# One expression per centred design column, in the order of `scaled_center`.
# Numeric terms, products of numeric variables, and factor, character, and
# logical main effects under treatment contrasts can be written out; anything
# else is refused rather than written approximately. Every expression is
# checked against the fit's own column names before it is used.
lin_covariate_exprs <- function(object, trms, covar_labels) {
  data_classes <- attr(trms, "dataClasses")
  factors <- attr(trms, "factors")
  xlevels <- object[["xlevels"]]
  contrasts <- object[["contrasts"]]

  exprs <- character(0)
  for (label in covar_labels) {
    vars <- rownames(factors)[factors[, label] > 0]
    classes <- data_classes[vars]
    if (all(classes == "numeric")) {
      exprs[label] <- if (length(vars) == 1L) {
        label
      } else {
        paste0("(", paste(vars, collapse = " * "), ")")
      }
    } else if (length(vars) == 1L && classes == "logical") {
      exprs[paste0(label, "TRUE")] <- paste0("(", label, ")")
    } else if (length(vars) == 1L &&
               classes %in% c("factor", "character") &&
               contrasts_are_treatment(contrasts[[vars]])) {
      # Usually every level but the first; all of them when R codes the factor
      # in full, as it does for the first factor of a fit with no intercept.
      levs <- xlevels[[vars]]
      levs <- levs[paste0(label, levs) %in% names(object[["scaled_center"]])]
      exprs[paste0(label, levs)] <-
        sprintf("(%s == %s)", label, encodeString(levs, quote = "\""))
    } else {
      stop(lin_terms_refusal(label), call. = FALSE)
    }
  }

  if (!identical(names(exprs), names(object[["scaled_center"]]))) {
    stop(lin_terms_refusal(NULL), call. = FALSE)
  }
  unname(exprs)
}

contrasts_are_treatment <- function(contr) {
  is.null(contr) || identical(contr, "contr.treatment")
}

lin_terms_refusal <- function(label) {
  paste0(
    "`terms()` cannot write this `lm_lin()` design out as a formula: ",
    if (is.null(label)) {
      "the written-out columns do not match the fit's centred covariates."
    } else {
      paste0(
        "the covariate `", label, "` is not a numeric variable, a product of ",
        "numeric variables, or a factor under treatment contrasts, so its ",
        "centred columns have no exact expression."
      )
    },
    " `model.matrix()` on the fit returns the design itself."
  )
}

# Seventeen significant digits round-trip a double exactly.
format_constant <- function(x) sprintf("%.17g", x)

#' @export
labels.lm_robust <- function(object, ...) attr(terms(object), "term.labels")

#' @export
labels.iv_robust <- labels.lm_robust

# performance 0.17.0 registers a Gaussian logLik() for iv_robust fits, so
# estimatr does not: whichever package loaded last would decide what AIC(fit)
# returns.
# These rest on a likelihood or on leave-one-out refits that estimatr does not
# compute, and lm's versions assume the homoskedastic model a robust fit does
# not. Without a method each failed with "no applicable method"; now each says
# why.
refuse_method <- function(what, why, generic = function(model, ...) NULL,
                          fits = "estimatr fits") {
  force(what)
  force(why)
  force(fits)
  f <- function(...) {
    stop("`", what, "()` is not available for ", fits, ": ", why,
         call. = FALSE)
  }
  # R CMD check requires a method to carry its generic's arguments
  formals(f) <- formals(generic)
  f
}

no_likelihood <- paste(
  "the estimators make no distributional assumption about the errors, so",
  "there is no likelihood, and no AIC or BIC built on one."
)
no_influence <- paste(
  "case-deletion diagnostics are not computed. `hatvalues()` gives the",
  "leverages; refit with `lm()` for the rest."
)

#' @importFrom stats logLik
#' @export
logLik.lm_robust <- refuse_method("logLik", no_likelihood, stats::logLik)
#' @importFrom stats extractAIC
#' @export
extractAIC.lm_robust <- refuse_method("extractAIC", no_likelihood,
                                      stats::extractAIC)

#' @importFrom stats influence
#' @export
influence.lm_robust <- refuse_method("influence", no_influence)
#' @importFrom stats cooks.distance
#' @export
cooks.distance.lm_robust <- refuse_method("cooks.distance", no_influence)
#' @importFrom stats rstandard
#' @export
rstandard.lm_robust <- refuse_method("rstandard", no_influence)
#' @importFrom stats rstudent
#' @export
rstudent.lm_robust <- refuse_method("rstudent", no_influence)
#' @importFrom stats dfbeta
#' @export
dfbeta.lm_robust <- refuse_method("dfbeta", no_influence)
#' @importFrom stats dfbetas
#' @export
dfbetas.lm_robust <- refuse_method("dfbetas", no_influence)
#' @export
influence.iv_robust <- influence.lm_robust
#' @export
cooks.distance.iv_robust <- cooks.distance.lm_robust
#' @export
rstandard.iv_robust <- rstandard.lm_robust
#' @export
rstudent.iv_robust <- rstudent.lm_robust
#' @export
dfbeta.iv_robust <- dfbeta.lm_robust
#' @export
dfbetas.iv_robust <- dfbetas.lm_robust

# A difference in means stores neither residuals nor a design, so the defaults
# returned NULL, a vector of length zero, or (labels) the names of the fit's
# fields, which a consumer took for a value.
dim_fits <- "`difference_in_means()` fits"
dim_no_residuals <- paste(
  "the estimator has no residuals or fitted values. The estimate and its",
  "standard error are in `tidy(fit)`."
)
dim_no_design <- "the estimator is not built from a design matrix."

#' @importFrom stats residuals
#' @export
residuals.difference_in_means <- refuse_method("residuals", dim_no_residuals,
                                               stats::residuals, dim_fits)
#' @importFrom stats fitted
#' @export
fitted.difference_in_means <- refuse_method("fitted", dim_no_residuals,
                                            stats::fitted, dim_fits)
#' @export
deviance.difference_in_means <- refuse_method("deviance", dim_no_residuals,
                                              stats::deviance, dim_fits)
#' @export
sigma.difference_in_means <- refuse_method("sigma", dim_no_residuals,
                                           stats::sigma, dim_fits)
#' @export
variable.names.difference_in_means <- refuse_method("variable.names",
                                                    dim_no_design,
                                                    stats::variable.names,
                                                    dim_fits)
#' @export
labels.difference_in_means <- refuse_method("labels", dim_no_design,
                                            base::labels, dim_fits)
