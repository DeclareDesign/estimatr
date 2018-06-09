# Internal method takes the results and adds p-values and confidence intervals
add_cis_pvals <- function(return_frame, alpha, ci, ttest = TRUE) {
  if (ci) {
    if (alpha <= 0 || alpha >= 1) {
      stop("`alpha` must be numeric between 0 and 1")
    }

    return_frame$statistic <- with(return_frame, coefficients / std.error)

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
        2 * pt(abs(statistic), df = df, lower.tail = FALSE)
      )

      crit_se <- with(return_frame, qt(1 - alpha / 2, df = df) * std.error)
    } else {
      return_frame$p.value <- with(
        return_frame,
        2 * pnorm(abs(statistic), lower.tail = FALSE)
      )

      crit_se <- with(return_frame, qnorm(1 - alpha / 2) * std.error)

      return_frame$df <- NA
    }

    return_frame$conf.low <- with(return_frame, coefficients - crit_se)
    return_frame$conf.high <- with(return_frame, coefficients + crit_se)

    return(as.list(return_frame))
  } else {
    return_frame$p.value <- NA
    return_frame$statistic <- NA
    return_frame$conf.low <- NA
    return_frame$conf.high <- NA

    return(as.list(return_frame))
  }
}
