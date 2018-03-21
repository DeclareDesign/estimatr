# Internal method takes the results and adds p-values and confidence intervals
add_cis_pvals <- function(return_frame, alpha, ci, ttest = TRUE) {
  if (ci) {
    if (alpha <= 0 || alpha >= 1) {
      stop("`alpha` must be numeric between 0 and 1")
    }

    if (ttest) {
      if (any(return_frame$df <= 0, na.rm = TRUE)) {
        warning(
          "Some degrees of freedom have been estimated as negative or zero\n",
          "p-values and confidence intervals may not be calculated"
        )

        return_frame$df <- ifelse(return_frame$df <= 0, NA, return_frame$df)
      }

      return_frame$p.value <- with(
        return_frame,
        2 * pt(abs(coefficients / std.error), df = df, lower.tail = FALSE)
      )

      crit_se <- with(return_frame, qt(1 - alpha / 2, df = df) * std.error)
    } else {
      return_frame$p.value <- with(
        return_frame,
        2 * pnorm(abs(coefficients / std.error), lower.tail = FALSE)
      )

      crit_se <- with(return_frame, qnorm(1 - alpha / 2) * std.error)

      return_frame$df <- NA
    }

    return_frame$ci.lower <- with(return_frame, coefficients - crit_se)
    return_frame$ci.upper <- with(return_frame, coefficients + crit_se)

    return(as.list(return_frame))
  } else {
    return_frame$p.value <- NA
    return_frame$ci.lower <- NA
    return_frame$ci.upper <- NA

    return(as.list(return_frame))
  }
}
