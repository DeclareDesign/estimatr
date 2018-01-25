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
        "res_var",
        "weighted",
        "se_type"
      )]

    return_list[["coefficients"]] <- summarize_tidy(object)
    return_list[["N"]] <- nobs(object)

    class(return_list) <- "summary.lm_robust"

    return_list
  }


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
    c(
      "Estimate",
      "Std. Error",
      paste0("Pr(>|", test, "|)"),
      "CI Lower",
      "CI Upper",
      "DF"
    )
  tidy_mat <- as.matrix(tidy_out[, !(names(tidy_out) %in% remove_cols)])
  rownames(tidy_mat) <- tidy_out$coefficient_name

  return(tidy_mat)
}
