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
                                  signif.stars = getOption("show.signif.stars"),
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

  fstat <- if (is.numeric(x[["fstatistic"]])) {
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

  if (is.numeric(x[["proj_fstatistic"]])) {
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

  if (is.numeric(x[["diagnostic_endogeneity_test"]])) {
    cat("\nDiagnostics:\n")
    printCoefmat(
      build_ivreg_diagnostics_mat(x),
      cs.ind = 1L:2L,
      tst.ind = 3L,
      has.Pvalue = TRUE,
      P.values = TRUE,
      digits = digits,
      signif.stars = signif.stars,
      na.print = "NA",
      ...
    )
  }
  invisible(x)
}

#' @export
print.summary.lm_robust <- function(x,
                                    digits = max(3L, getOption("digits") - 3L),
                                    signif.stars = getOption("show.signif.stars"),
                                    ...) {
  print_summary_lm_like(x, digits, ...)
}

#' @export
print.summary.iv_robust <- function(x,
                                    digits = max(3L, getOption("digits") - 3L),
                                    signif.stars = getOption("show.signif.stars"),
                                    ...) {
  print_summary_lm_like(x, digits, signif.stars, ...)
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

#' @export
print.lh <- function(x, ...) {
  print(summary(x))
}

#' @export
print.lh_robust <- function(x, ...) {
  lnames <- names(x)
  for (i in seq_along(x)) {
    cat("$", lnames[i], "\n", sep = "")
    print(x[[i]])
    cat("\n")
  }
  invisible(x)
}

#' @export
print.summary.lh_robust <- function(x,
                                    digits = max(3L, getOption("digits") - 3L),
                                    ...){
  lnames <- names(x)
  for (i in seq_along(x)) {
    cat("$", lnames[i], "\n", sep = "")
    print(summary(x[[i]]), digits = digits)
    cat("\n")
  }
}

#' @export
print.summary.lh <- function(x,
                             digits = max(3L, getOption("digits") - 3L),
                             ...){
  class(x) <- NULL
  print(x, digits = digits)
}
