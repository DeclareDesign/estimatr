#' Summarizing \code{\link{lm_robust}} objects
#'
#' @param object an object of class \code{\link{lm_robust}}
#' @param ... arguments passed to \code{\link{tidy.lm_robust}}, unused
#'
#' @export
summary.lm_robust <-
  function(
           object,
           ...) {

    return_list <-
      object[c(
        "call",
        "k",
        "rank",
        "df.residual",
        "r.squared",
        "adj.r.squared",
        "fstatistic",
        "weighted",
        "se_type"
      )]

    return_list[["coefficients"]] = summarize_tidy(object)
    return_list[["N"]] = nobs(object)

    class(return_list) <- "summary.lm_robust"

    print(return_list)
  }


#' Summarizing \code{\link{difference_in_means}} objects
#'
#' @param object an object of class \code{\link{difference_in_means}}
#' @param ... arguments passed to \code{\link{tidy.difference_in_means}}, unused
#'
#' @export
summary.difference_in_means <-
  function(
           object,
           ...) {
    return(list(
      coefficients = summarize_tidy(object),
      design = object$design
    ))
  }


#' Summarizing \code{\link{horvitz_thompson}} objects
#'
#' @param object an object of class \code{\link{horvitz_thompson}}
#' @param ... arguments passed to \code{\link{tidy.horvitz_thompson}}, unused
#'
#' @export
summary.horvitz_thompson <-
  function(
           object,
           ...) {
    return(list(coefficients = summarize_tidy(object, "z")))
  }

summarize_tidy <- function(object, test = "t", ...) {
  remove_cols <- c("coefficient_name", "outcome")

  # This is ugly SO THAT summary(fit)$coefficients returns something like lm does.
  tidy_out <- tidy(object, ...)
  colnames(tidy_out)[2:7] <-
    c("Estimate",
      "Std. Error",
      paste0("Pr(>|", test, "|)"),
      "CI Lower",
      "CI Upper",
      "DF")
  tidy_mat <- as.matrix(tidy_out[, !(names(tidy_out) %in% remove_cols)])
  rownames(tidy_mat) <- tidy_out$coefficient_name

  return(tidy_mat)
}
