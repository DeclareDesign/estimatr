# Internal method takes the results and adds pvalues and confidence intervals
add_cis_pvals <- function(return_frame, alpha) {

  if (alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be numeric between 0 and 1")
  }

  return_frame$p <- with(return_frame,
                         2 * pt(abs(est / se), df = df, lower.tail = FALSE))
  return_frame$ci_lower <- with(return_frame,
                                est - qt(1 - alpha / 2, df = df) * se)
  return_frame$ci_upper <- with(return_frame,
                                est + qt(1 - alpha / 2, df = df) * se)

  return(as.list(return_frame))
}
