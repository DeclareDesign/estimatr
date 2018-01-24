#' Printing \code{\link{lm_robust}} objects
#'
#' @param x an object of class \code{\link{lm_robust}}
#' @param ... arguments passed to \code{\link{tidy.lm_robust}}, unused
#'
#' @export
print.lm_robust <-
  function(
           x,
           ...) {
    print(summarize_tidy(x))
  }

#' Printing \code{\link{summary.lm_robust}} objects
#'
#' @param x an object of class \code{\link{lm_robust}}
#' @param ... arguments passed to \code{\link{tidy.lm_robust}}, unused
#'
#' @export
print.summary.lm_robust <-
  function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    ...) {

    cat("\nCall:\n", paste(deparse(x$call, nlines = 5), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    if (x$weighted) {
      cat("Weighted, ")
    }
    cat("Standard error type = ", x$se_type, "\n")

    if (x$rank < x$k) {
      singularities <- x$k - x$rank
      cat("\nCoefficients: (", nsingular, " not defined the design matrix is rank deficient)\n",
          sep = "")
    } else {
      cat("\nCoefficients:\n")
    }

    print(x$coefficients, digits = digits)

    if (!is.null(x$fstatistic)) {

      cat(
        "\nMultiple R-squared: ", formatC(x$r.squared, digits = digits),
        ",\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits),
        "\nF-statistic:", formatC(x$fstatistic[1L],digits = digits),
        "on", x$fstatistic[2L], "and", x$fstatistic[3L],
        "DF,  p-value:",
        format.pval(pf(
          x$fstatistic[1L],
          x$fstatistic[2L],
          x$fstatistic[3L],
          lower.tail = FALSE
        ), digits = digits)
      )
      cat("\n")
    }

    invisible(x)
  }


#' Printing \code{\link{difference_in_means}} objects
#'
#' @param x an object of class \code{\link{difference_in_means}}
#' @param ... arguments passed to \code{\link{tidy.difference_in_means}}, unused
#'
#' @export
print.difference_in_means <-
  function(
           x,
           ...) {
    cat("Design: ", x$design, "\n")
    print(summarize_tidy(x))
  }


#' Printing \code{\link{horvitz_thompson}} objects
#'
#' @param x an object of class \code{\link{horvitz_thompson}}
#' @param ... arguments passed to \code{\link{tidy.horvitz_thompson}}, unused
#'
#' @export
print.horvitz_thompson <-
  function(x,
           ...) {
    print(summarize_tidy(x))
  }
