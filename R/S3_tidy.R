#' @importFrom generics tidy
#' @export
generics::tidy

tidy_data_frame <- function(x, ...) {
  vec_cols <- c(
    "coefficients",
    "std.error",
    "statistic",
    "p.value",
    "conf.low",
    "conf.high",
    "df"
  )

  tidy_mat <- do.call("cbind", lapply(x[vec_cols], as.vector))
  vec_cols[vec_cols == "coefficients"] <- "estimate"
  colnames(tidy_mat) <- vec_cols
  return_frame <- data.frame(
    term = x[["term"]],
    tidy_mat,
    outcome = rep(x[["outcome"]], each = length(x[["term"]])),
    stringsAsFactors = FALSE
  )

  rownames(return_frame) <- NULL

  return(return_frame)
}

warn_singularities <- function(x) {
  if (x$rank < x$k) {
    singularities <- x$k - x$rank
    what <- ifelse(singularities > 1, " coefficients ", " coefficient ")
    message(
      singularities, what,
      " not defined because the design matrix is rank deficient\n"
    )
  }
}

#' Tidy an estimatr object
#' @name estimatr_tidiers
#' @templateVar class lm_robust
#' @return A data.frame with columns for coefficient names, estimates, standard
#' errors, confidence intervals, p-values, degrees of freedom, and the
#' name of the outcome variable
#'
#' @param x An object returned by one of the estimators
#' @param ... extra arguments (not used)
#'
#' @export
#' @family estimatr tidiers
#' @seealso [generics::tidy()], [estimatr::lm_robust()], [estimatr::iv_robust()],  [estimatr::difference_in_means()], [estimatr::horvitz_thompson()]
#' @md
tidy.lm_robust <- function(x, ...) {
  warn_singularities(x)
  tidy_data_frame(x)
}

#' @rdname estimatr_tidiers
#' @templateVar class iv_robust
#'
#' @export
#' @family estimatr tidiers
tidy.iv_robust <- function(x, ...) {
  warn_singularities(x)
  tidy_data_frame(x)
}

#' @rdname estimatr_tidiers
#' @templateVar class difference_in_means
#'
#' @export
#' @family estimatr tidiers
tidy.difference_in_means <- tidy_data_frame

#' @rdname estimatr_tidiers
#' @templateVar class horvitz_thompson
#'
#' @export
#' @family estimatr tidiers
tidy.horvitz_thompson <- tidy_data_frame
