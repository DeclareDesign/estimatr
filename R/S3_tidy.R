# some code taken from the "broom" package
# https://github.com/tidyverse/broom

#' Tidy the result of an estimator into a data.frame
#' @rdname tidy
#'
#' @export
tidy <- function(object, ...) {
  if (requireNamespace("broom", quietly = TRUE))
    broom::tidy(object, ...)
  else
    UseMethod("tidy")
}


#' @rdname tidy
#'
#' @param object An object returned by one of the estimators
#' @param ... extra arguments (not used)
#'
#' @export
tidy.NULL <- function(object, ...) {
  data.frame()
}

#' @rdname tidy
#'
#' @export
tidy.default <- function(object, ...) {
  warning(paste(
    "No method for tidying an S3 object of class",
    class(object),
    ", using as.data.frame"
  ))
  as.data.frame(object)
}


#' @rdname tidy
#'
#' @return A data.frame with with coefficient names, estimates, standard
#' errors, confidence intervals, p-values, degrees of freedom, and the
#' name of the outcome variable
#'
#' @export tidy.lm_robust
#' @export
tidy.lm_robust <- function(object, ...) {
  warn_singularities(object)
  tidy_data_frame(object)
}

#' @rdname tidy
#'
#' @export tidy.iv_robust
#' @export
tidy.iv_robust <- function(object, ...) {
  warn_singularities(object)
  tidy_data_frame(object)
}


#' @rdname tidy
#'
#' @export tidy.difference_in_means
#' @export
tidy.difference_in_means <- function(object, ...) {
  tidy_data_frame(object)
}


#' @rdname tidy
#'
#' @export tidy.horvitz_thompson
#' @export
tidy.horvitz_thompson <- function(object, ...) {
  tidy_data_frame(object)
}


tidy_data_frame <- function(object, digits = NULL) {
  vec_cols <- c(
      "coefficients",
      "std.error",
      "statistic",
      "p.value",
      "conf.low",
      "conf.high",
      "df"
  )

  tidy_mat <- do.call("cbind", lapply(object[vec_cols], as.vector))
  vec_cols[vec_cols == "coefficients"] <- "estimate"
  colnames(tidy_mat) <- vec_cols
  return_frame <- data.frame(
    term = object[["term"]],
    tidy_mat,
    stringsAsFactors = FALSE
  )

  return_frame$outcome <- rep(object[["outcome"]], each = length(object[["term"]]))

  rownames(return_frame) <- NULL
  return(return_frame)
}

warn_singularities <- function(object) {
  if (object$rank < object$k) {
    singularities <- object$k - object$rank
    what <- ifelse(singularities > 1, " coefficients ", " coefficient ")
    message(
      singularities, what,
      " not defined because the design matrix is rank deficient\n"
    )
  }
}
