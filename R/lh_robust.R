#' Linear Hypothesis Test for OLS with Robust Standard Errors
#'
#' Tests a linear combination of coefficients, or several of them jointly,
#' from a model fitted by [lm_robust()]. The robust variance and the
#' degrees of freedom of the fit are carried through, so a clustered fit is
#' tested on its cluster-adjusted degrees of freedom rather than on the
#' residual ones.
#'
#' @param ... (optional) Other arguments passed to [lm_robust()]
#' @param data (optional) A `data.frame`
#' @param linear_hypothesis (required) A character string or matrix specifying the
#'   hypothesis, passed to `car::linearHypothesis`
#'
#' @return An object of class `"lh_robust"` with three components:
#'   `lm_robust`, the underlying fit; `lh`, one row per hypothesis holding
#'   `coefficients`, `std.error`, `statistic`, `p.value`, `alpha`, `conf.low`,
#'   `conf.high`, `df`, `term`, and `outcome`; and `joint_hypothesis`, the Wald
#'   F test of all of them at once, as `value`, `numdf`, `dendf`, and
#'   `p.value`. Under `se_type = "CR2"` each hypothesis's `df` is its own
#'   Satterthwaite approximation, as `clubSandwich::linear_contrast()` computes
#'   it, and `dendf` is the smallest of them.
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

  lmr <- lm_robust_hypotheses(..., data = data, linear_hypothesis = linear_hypothesis)

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

  ci <- eval_tidy(quos(...)$ci)
  if (is.null(ci)) {
    ci <- TRUE
  }

  # car::linearHypothesis reaches for the fit's vcov and refuses without one,
  # naming `return_vcov`, which is not the argument the caller set. A linear
  # combination of coefficients has no standard error when the fit has no
  # variance, so say that instead.
  if (identical(eval_tidy(quos(...)$se_type), "none")) {
    stop(
      "`lh_robust` needs a variance: `se_type = \"none\"` fits without one, ",
      "so a linear combination of the coefficients has no standard error. ",
      "Choose a variance estimator, or call `lm_robust(se_type = \"none\")` ",
      "for the coefficients alone."
    )
  }

  car_lht <- car::linearHypothesis(
    lmr, hypothesis.matrix = linear_hypothesis, level = 1 - alpha)

  estimate  <- drop(attr(car_lht, "value"))
  vcov_lh   <- attr(car_lht, "vcov")
  std.error <- sqrt(diag(vcov_lh))

  # Under CR2 each hypothesis has its own Satterthwaite degrees of freedom,
  # which the fit computed from the same components as the coefficients'. A
  # combination used to take the smallest per-coefficient df, described here as
  # conservative; it is not a bound, and on the data in test_vs_clubsandwich.R
  # it was 16.10 where the combination's is 15.75. Every other se_type gives
  # every coefficient the same df, so any combination takes that.
  df_vec <- if (!is.null(lmr[["hypothesis_df"]])) {
    setNames(lmr[["hypothesis_df"]], names(estimate))
  } else if (any(!is.na(lmr$df))) {
    setNames(rep(min(lmr$df, na.rm = TRUE), length(estimate)), names(estimate))
  } else {
    # `ci = FALSE` leaves the fit with no degrees of freedom, and min() of
    # nothing is Inf with a warning about missing arguments that says nothing
    # to the caller.
    setNames(rep(NA_real_, length(estimate)), names(estimate))
  }

  statistic <- estimate / std.error
  if (isTRUE(ci)) {
    p.value    <- 2 * pt(abs(statistic), df_vec, lower.tail = FALSE)
    half_width <- std.error * qt(1 - alpha / 2, df_vec)
    ci_low     <- estimate - half_width
    ci_high    <- estimate + half_width
  } else {
    p.value <- rep(NA_real_, length(estimate))
    ci_low  <- rep(NA_real_, length(estimate))
    ci_high <- rep(NA_real_, length(estimate))
  }

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

  # Joint Wald F-test: W = t(Lβ) (L Vcov L')^{-1} (Lβ) / m ~ F(m, df_joint).
  # The denominator takes the smallest of the hypotheses' df. That is a
  # convention: clubSandwich's small-sample joint test (HTZ) is a different
  # approximation, and nothing here claims to match it.
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
