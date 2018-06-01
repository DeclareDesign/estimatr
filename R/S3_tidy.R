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
  return_frame <- tidy_data_frame(object)

  warn_singularities(object)

  return(return_frame)
}

#' @rdname tidy
#'
#' @export tidy.iv_robust
#' @export
tidy.iv_robust <- function(object, ...) {
  return_frame <- tidy_data_frame(object)

  warn_singularities(object)

  return(return_frame)
}


#' @rdname tidy
#'
#' @export tidy.difference_in_means
#' @export
tidy.difference_in_means <- function(object, ...) {
  return_frame <- tidy_data_frame(object)
  return(return_frame)
}


#' @rdname tidy
#'
#' @export tidy.horvitz_thompson
#' @export
tidy.horvitz_thompson <- function(object, ...) {
  return_frame <- tidy_data_frame(object)
  return(return_frame)
}


tidy_data_frame <- function(object, digits = NULL) {
  vec_cols <-
    c(
      "coefficients",
      "std.error",
      "p.value",
      "ci.lower",
      "ci.upper",
      "df"
    )

  tidy_mat <- do.call("cbind", lapply(vec_cols, function(x) {
    as.vector(object[[x]])
  }))
  vec_cols[which(vec_cols == "coefficients")] <- "estimate"
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
    plural <- ifelse(singularities > 1, "s", "")
    message(
      singularities, " coefficient", plural,
      " not defined because the design matrix is rank deficient\n"
    )
  }
}
