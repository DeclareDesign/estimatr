#' Linear Hypothesis Test for OLS with Robust Standard Errors
#'
#' @param ... Other arguments passed to [lm_robust()]
#' @param data A `data.frame`
#' @param linear_hypothesis A character string or matrix specifying the
#'   hypothesis, passed to `car::linearHypothesis`
#'
#' @return An object of class `"lh_robust"` with components `lm_robust` and `lh`.
#'
#' @importFrom rlang quos eval_tidy
#' @export
lh_robust <- function(..., data, linear_hypothesis) {

  requireNamespace("car")

  lmr <- lm_robust(..., data = data)

  alpha <- eval_tidy(quos(...)$alpha)
  if (is.null(alpha)) {
    alpha <- 0.05
  }

  car_lht <- car::linearHypothesis(
    lmr, hypothesis.matrix = linear_hypothesis, level = 1 - alpha)

  estimate <- drop(attr(car_lht, "value"))
  std.error <- sqrt(diag(attr(car_lht, "vcov")))

  df <- lmr$df.residual

  statistic <- estimate / std.error
  p.value <- 2 * pt(abs(statistic), df, lower.tail = FALSE)
  ci <- estimate + std.error %o% qt(c(alpha / 2, 1 - alpha / 2), df)

  return_lh_robust <- data.frame(
    coefficients = estimate,
    std.error = std.error,
    statistic = statistic,
    p.value = p.value,
    alpha = alpha,
    conf.low = ci[, 1],
    conf.high = ci[, 2],
    df = df,
    term = linear_hypothesis,
    outcome = lmr$outcome
  )

  attr(return_lh_robust, "linear_hypothesis") <- car_lht
  class(return_lh_robust) <- c("lh", "data.frame")

  return_lmr <- lmr
  return_lmr[["call"]] <- match.call()

  return(structure(
    list(lm_robust = return_lmr, lh = return_lh_robust),
    class = "lh_robust"
  ))

}
