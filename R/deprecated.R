#' Functions removed in 2.0
#'
#' `commarobust()` and `starprep()` were helpers for producing robust standard
#' errors outside the package's own estimators and formatting them for
#' \pkg{stargazer}. Both are removed.
#'
#' They are kept here as names that error rather than deleted outright, so that
#' a script written against estimatr 1.x says what happened and what to do
#' instead of failing with `could not find function`.
#'
#' `commarobust()` recomputed robust standard errors on a fitted `lm`. Fit the
#' model with [lm_robust()] instead, which is what it was reimplementing.
#'
#' `starprep()` prepared a list of fits for \pkg{stargazer}, which has not been
#' maintained for years. Table-building now goes through
#' \pkg{modelsummary}, which reads [tidy()] and [glance()] and therefore
#' works on every estimator in this package without any adapter.
#'
#' @param ... Ignored.
#' @return Never returns; both functions signal an error.
#' @name estimatr-defunct
#' @examples
#' # Both of these error. The replacements:
#' set.seed(1)
#' dat <- data.frame(y = rnorm(20), z = rep(0:1, 10))
#'
#' # was: commarobust(lm(y ~ z, data = dat))
#' lm_robust(y ~ z, data = dat)
#'
#' # was: starprep(fit1, fit2) |> stargazer::stargazer()
#' # now: modelsummary::modelsummary(list(fit1, fit2))
NULL

#' @rdname estimatr-defunct
#' @export
commarobust <- function(...) {
  stop(
    "`commarobust()` was removed in estimatr 2.0.\n",
    "It recomputed robust standard errors on a fitted `lm`. ",
    "Fit the model with `lm_robust()` instead:\n",
    "  lm_robust(y ~ x, data = dat, se_type = \"HC2\", clusters = cl)",
    call. = FALSE
  )
}

#' @rdname estimatr-defunct
#' @export
starprep <- function(...) {
  stop(
    "`starprep()` was removed in estimatr 2.0.\n",
    "It prepared fits for stargazer, which is no longer maintained. ",
    "Use modelsummary, which reads `tidy()` and `glance()` and so works on ",
    "every estimator in this package:\n",
    "  modelsummary::modelsummary(list(fit1, fit2))",
    call. = FALSE
  )
}
