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

        return_frame$dof <- ifelse(return_frame$df <= 0, NA, return_frame$df)
      } else {
        return_frame$dof <- return_frame$df
      }

      return_frame$p <- with(
        return_frame,
        2 * pt(abs(coefficients / se), df = dof, lower.tail = FALSE)
      )

      crit_se <- with(return_frame, qt(1 - alpha / 2, df = dof) * se)
    } else {
      return_frame$p <- with(
        return_frame,
        2 * pnorm(abs(coefficients / se), lower.tail = FALSE)
      )

      crit_se <- with(return_frame, qnorm(1 - alpha / 2) * se)

      return_frame$df <- NA
    }

    return_frame$ci_lower <- with(return_frame, coefficients - crit_se)
    return_frame$ci_upper <- with(return_frame, coefficients + crit_se)

    if (is.data.frame(return_frame)) {
      return_frame <- return_frame[, !(names(return_frame) == "dof"), drop = FALSE]
    }

    return(as.list(return_frame))
  } else {
    return_frame$p <- NA
    return_frame$ci_lower <- NA
    return_frame$ci_upper <- NA

    return(as.list(return_frame))
  }
}
