#' @importFrom Formula Formula
#' @importFrom stats getCall
#' @export
update.iv_robust <- function(object, formula., ..., evaluate = TRUE) {
  if (is.null(call <- getCall(object)))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.))
    call$formula <- formula(update(Formula(formula(object)), formula.))
  call <- update_extras(call, extras)
  if (evaluate)
    eval(call, parent.frame())
  else call
}

# formula() on an lm_lin() fit is `y ~ z + x`, the treatment and the
# covariates together, so update.default() expanded `. ~ . + w` against it and
# handed the whole right-hand side to `formula`, which lm_lin() refuses with a
# message about the user's own formula. The updated formula is split as
# lm_lin() reads it: the treatment stays in `formula` and every other term is
# the new `covariates`.
#' @export
update.lm_lin <- function(object, formula., ..., evaluate = TRUE) {
  call <- getCall(object)
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) {
    new <- update(formula(object), formula.)
    new_terms <- terms(new)
    labels <- attr(new_terms, "term.labels")
    treatment <- attr(object[["terms"]], "term.labels")[1L]
    if (!treatment %in% labels) {
      stop("`update()` cannot remove the treatment `", treatment, "` from an ",
           "`lm_lin()` fit: `formula` must keep it on the right-hand side.",
           call. = FALSE)
    }
    covariates <- setdiff(labels, treatment)
    if (length(covariates) == 0) {
      stop("`update()` cannot remove every covariate from an `lm_lin()` fit: ",
           "with nothing to centre, the model is `lm_robust(",
           deparse1(new), ")`.", call. = FALSE)
    }
    env <- environment(new)
    call$formula <- stats::as.formula(
      paste(deparse1(new[[2L]]), "~",
            paste(c(if (attr(new_terms, "intercept") == 0L) "0", treatment),
                  collapse = " + ")),
      env = env
    )
    call$covariates <- stats::as.formula(
      paste("~", paste(covariates, collapse = " + ")),
      env = env
    )
  }
  call <- update_extras(call, extras)
  if (evaluate)
    eval(call, parent.frame())
  else call
}

# An lh_robust() fit's call is kept on its lm_robust part, which is where the
# hypothesis table's coefficients come from.
#' @export
getCall.lh_robust <- function(x, ...) x[["lm_robust"]][["call"]]

# The formula reaches lm_robust() through lh_robust()'s dots, so the stored
# call may carry it unnamed; it is then the first unnamed argument.
#' @export
update.lh_robust <- function(object, formula., ..., evaluate = TRUE) {
  call <- getCall(object)
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) {
    arg_names <- names(call)
    slot <- if ("formula" %in% arg_names) {
      "formula"
    } else {
      which(c(FALSE, arg_names[-1L] == ""))[1L]
    }
    call[[slot]] <- update(formula(object), formula.)
  }
  call <- update_extras(call, extras)
  if (evaluate)
    eval(call, parent.frame())
  else call
}

# update.default()'s treatment of the arguments after `formula.`: one the call
# already has is replaced, a new one is appended.
update_extras <- function(call, extras) {
  if (!length(extras)) return(call)
  existing <- !is.na(match(names(extras), names(call)))
  for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
  if (any(!existing)) {
    call <- as.call(c(as.list(call), extras[!existing]))
  }
  call
}
