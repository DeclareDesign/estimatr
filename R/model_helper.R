
## will this work for lm, glm, ...?
#' @export
use_model <- function(...,
                      model = lm,
                      data,
                      coefficient_name = "Z") {

  fit <- do.call(model, args = c(list(...), list(data = data)))

  coef <- coef(fit)
  se <- sqrt(diag(vcov(fit)))

  variable_names <- names(coef)

  ## how to generally get the p-value and CI?

  ##p <- 2 * pt(abs(coef), df = df, lower.tail = FALSE)
  ##ci_lower <- coef - qt(1 - alpha / 2, df = df) * se
  ##ci_upper <- coef + qt(1 - alpha / 2, df = df) * se

  return_frame <-
    data.frame(
      variable_names = variable_names,
      est = coef,
      se = se ##,
      ##p = p,
      ##ci_lower = ci_lower,
      ##ci_upper = ci_upper
    )

  return(subset(return_frame, variable_names %in% coefficient_name))

}
