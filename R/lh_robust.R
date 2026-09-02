#' Linear Hypothesis Test for OLS with Robust Standard Errors
#'
#' Tests a linear combination of coefficients, or several of them jointly,
#' from a model fitted by [lm_robust()]. The robust variance and the
#' degrees of freedom of the fit are carried through, so a clustered fit is
#' tested on its cluster-adjusted degrees of freedom rather than on the
#' residual ones.
#'
#' @param ... Other arguments passed to [lm_robust()]
#' @param data A `data.frame`
#' @param linear_hypothesis A character string or matrix specifying the
#'   hypothesis, passed to `car::linearHypothesis`
#'
#' @return An object of class `"lh_robust"` with three components:
#'   `lm_robust`, the underlying fit; `lh`,
#'   and `joint_hypothesis`.
#'
#' @importFrom rlang quos eval_tidy
#' @examples
#' set.seed(35)
#' dat <- data.frame(x = rnorm(100), z = rbinom(100, 1, 0.5),
#'                   cl = rep(1:10, each = 10))
#' dat$y <- dat$x + 0.5 * dat$z + rnorm(100)
#'
#' # One linear combination of coefficients
#' fit <- lh_robust(y ~ x + z, data = dat, linear_hypothesis = "z + 2*x = 0")
#' fit
#' tidy(fit)
#'
#' # Degrees of freedom follow the fit, so a clustered model tests against the
#' # cluster-adjusted df rather than the residual df
#' lh_robust(y ~ x + z, data = dat, clusters = cl,
#'           linear_hypothesis = "z + 2*x = 0")
#'
#' # Several restrictions at once give one joint Wald test as well
#' joint <- lh_robust(y ~ x + z, data = dat, linear_hypothesis = c("x = 0", "z = 0"))
#' joint$joint_hypothesis
#'
#' @export
lh_robust <- function(..., data, linear_hypothesis) {

  requireNamespace("car")

  lmr <- lm_robust(..., data = data)

  # With several outcomes the coefficients are named "<outcome>:<term>", so
  # car::linearHypothesis cannot match a hypothesis written in terms of the
  # variables and reports it as malformed (estimatr #297).
  if (length(lmr[["outcome"]]) > 1) {
    stop(
      "`lh_robust` does not support multiple outcomes: coefficients of a ",
      "multivariate fit are named '<outcome>:<term>', which a hypothesis ",
      "such as 'cyl = 2' cannot refer to. Fit one outcome at a time."
    )
  }

  alpha <- eval_tidy(quos(...)$alpha)
  if (is.null(alpha)) {
    alpha <- 0.05
  }

  car_lht <- car::linearHypothesis(
    lmr, hypothesis.matrix = linear_hypothesis, level = 1 - alpha)

  estimate  <- drop(attr(car_lht, "value"))
  vcov_lh   <- attr(car_lht, "vcov")
  std.error <- sqrt(diag(vcov_lh))

  # Resolve the df for each hypothesis. Hypothesis names look like "x=0" or
  # "x - z=0". Extract the LHS, check if it's a single coefficient name, and
  # use that coefficient's (Satterthwaite-adjusted) df. For complex combinations
  # fall back to the conservative min across all per-coefficient dfs.
  coef_dfs <- lmr$df
  df_vec <- vapply(names(estimate), function(nm) {
    lhs <- trimws(sub("\\s*=.*", "", nm))
    if (lhs %in% names(coef_dfs)) {
      unname(coef_dfs[lhs])
    } else {
      min(coef_dfs, na.rm = TRUE)
    }
  }, numeric(1))

  statistic  <- estimate / std.error
  p.value    <- 2 * pt(abs(statistic), df_vec, lower.tail = FALSE)
  half_width <- std.error * qt(1 - alpha / 2, df_vec)
  ci_low     <- estimate - half_width
  ci_high    <- estimate + half_width

  return_lh_robust <- data.frame(
    coefficients = estimate,
    std.error    = std.error,
    statistic    = statistic,
    p.value      = p.value,
    alpha        = alpha,
    conf.low     = unname(ci_low),
    conf.high    = unname(ci_high),
    df           = df_vec,
    term         = linear_hypothesis,
    outcome      = lmr$outcome
  )
  attr(return_lh_robust, "linear_hypothesis") <- car_lht
  class(return_lh_robust) <- c("lh", "data.frame")

  # Joint Wald F-test: W = t(Lβ) (L Vcov L')^{-1} (Lβ) / m ~ F(m, df_joint)
  m          <- length(estimate)
  wald       <- drop(t(estimate) %*% solve(vcov_lh) %*% estimate)
  joint_F    <- wald / m
  df_joint   <- min(df_vec)
  joint_pval <- pf(joint_F, m, df_joint, lower.tail = FALSE)
  joint_hypothesis <- c(
    value   = joint_F,
    numdf   = m,
    dendf   = df_joint,
    p.value = joint_pval
  )

  return_lmr <- lmr
  return_lmr[["call"]] <- match.call()

  return(structure(
    list(lm_robust = return_lmr, lh = return_lh_robust, joint_hypothesis = joint_hypothesis),
    class = "lh_robust"
  ))

}
