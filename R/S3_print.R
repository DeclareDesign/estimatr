#' @export
print.lm_robust <- function(x, ...) {
  print(summarize_tidy(x))
}

#' @export
print.iv_robust <- function(x, ...) {
  print(summarize_tidy(x))
}

print_summary_lm_like <- function(x,
                                  digits,
                                  ...) {
  cat(
    "\nCall:\n",
    paste(deparse(x$call, nlines = 5), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )
  if (x$weighted) {
    cat("Weighted, ")
  }
  cat("Standard error type: ", x$se_type, "\n")

  if (x$rank < x$k) {
    singularities <- x$k - x$rank
    cat(
      "\nCoefficients: (",
      singularities,
      " not defined because the design matrix is rank deficient)\n",
      sep = ""
    )
  } else {
    cat("\nCoefficients:\n")
  }

  print(coef(x), digits = digits)

  fstat <- if (is.numeric(x$fstatistic)) {
    paste(
      "\nF-statistic:", formatC(x$fstatistic[1L], digits = digits),
      "on", x$fstatistic[2L], "and", x$fstatistic[3L],
      "DF,  p-value:",
      format.pval(pf(
        x$fstatistic[1L],
        x$fstatistic[2L],
        x$fstatistic[3L],
        lower.tail = FALSE
      ), digits = digits)
    )
  } else NULL

  cat(
    "\nMultiple R-squared: ", formatC(x$r.squared, digits = digits),
    ",\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits),
    fstat
  )

  if (!is.null(x$proj_fstatistic)) {
    cat(
      "\nMultiple R-squared (proj. model): ",
      formatC(x$proj_r.squared, digits = digits),
      ",\tAdjusted R-squared (proj. model): ",
      formatC(x$proj_adj.r.squared, digits = digits),
      "\nF-statistic (proj. model):",
      formatC(x$proj_fstatistic[1L], digits = digits),
      "on", x$proj_fstatistic[2L], "and", x$proj_fstatistic[3L],
      "DF,  p-value:",
      format.pval(pf(
        x$proj_fstatistic[1L],
        x$proj_fstatistic[2L],
        x$proj_fstatistic[3L],
        lower.tail = FALSE
      ), digits = digits)
    )
  }
  cat("\n")

  invisible(x)
}

#' @export
print.summary.lm_robust <- function(x,
                                    digits = max(3L, getOption("digits") - 3L),
                                    ...) {
  print_summary_lm_like(x, digits, ...)
}

#' @export
print.summary.iv_robust <- function(x,
                                    digits = max(3L, getOption("digits") - 3L),
                                    ...) {
  print_summary_lm_like(x, digits, ...)
}

#' @export
print.difference_in_means <- function(x, ...) {
  cat("Design: ", x$design, "\n")
  print(summarize_tidy(x))
}


#' @export
print.horvitz_thompson <- function(x, ...) {
  print(summarize_tidy(x))
}
