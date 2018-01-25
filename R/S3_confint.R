#' @export
confint.lm_robust <-
  function(
           object,
           parm = NULL,
           level = NULL,
           ...) {
    cis <- get_ci_mat(object, level)

    if (!is.null(parm)) {
      cis <- cis[parm, , drop = FALSE]
    }

    return(cis)
  }


#' @export
confint.difference_in_means <-
  function(
           object,
           parm = NULL,
           level = NULL,
           ...) {
    cis <- get_ci_mat(object, level)

    return(cis)
  }

#' @export
confint.horvitz_thompson <-
  function(
           object,
           parm = NULL,
           level = NULL,
           ...) {
    cis <- get_ci_mat(object, level, ttest = FALSE)

    return(cis)
  }

## internal method that builds confidence intervals and labels the matrix to be returned
get_ci_mat <- function(object, level, ttest = TRUE) {
  if (!is.null(level)) {
    if (!is.null(object[["alpha"]])) {
      object[["alpha"]] <- NULL
    }
    object <- add_cis_pvals(object, alpha = 1 - level, ci = TRUE, ttest = ttest)
    cis <- cbind(object$ci_lower, object$ci_upper)
  } else {
    cis <- cbind(object$ci_lower, object$ci_upper)
    level <- 1 - object$alpha
  }

  dimnames(cis) <-
    list(
      object$coefficient_name,
      paste((1 - level) / 2 * c(100, -100) + c(0, 100), "%")
    )

  return(cis)
}
