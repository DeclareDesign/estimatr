
## will this work for lm, glm, ...?
#' Estimate and Summarize a Model
#'
#' @param ... arguments passed to model
#'
#' @param model a modeling function like \code{\link{lm}} or \code{\link{glm}}.  Must have coef() and vcov() methods.
#' @param data A data.frame.
#' @param alpha The significance level, 0.05 by default. (not yet implemented)
#' @param coefficient_name a character or character vector that indicates which coefficients should be reported. Defaults to "Z".
#'
#' @export
use_model <- function(...,
                      model = lm,
                      data,
                      alpha,
                      coefficient_name = "Z") {

  fit <- do.call(model, args = c(list(...), list(data = data)))

  coef <- coef(fit)
  se <- sqrt(diag(vcov(fit)))

  variable_names <- names(coef)

  ## how to generally get the p-value and CI?

  ##p <- 2 * pt(abs(coef / se), df = df, lower.tail = FALSE)
  ##ci_lower <- coef - qt(1 - alpha / 2, df = df) * se
  ##ci_upper <- coef + qt(1 - alpha / 2, df = df) * se

  return_frame <-
    data.frame(
      variable_names = variable_names,
      est = coef,
      se = se, 
      ##p = p,
      ##ci_lower = ci_lower,
      ##ci_upper = ci_upper
      stringsAsFactors = FALSE
    )

  return(subset(return_frame, variable_names %in% coefficient_name))

}
