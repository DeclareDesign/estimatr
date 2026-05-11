#' @importFrom generics tidy
#' @export
generics::tidy

tidy_data_frame <- function(x,
                            conf.int = TRUE,
                            conf.level = NULL,
                            ...) {
  vec_cols <- c(
    "coefficients",
    "std.error",
    "statistic",
    "p.value",
    "conf.low",
    "conf.high",
    "df"
  )

  if (!conf.int) {
    vec_cols <- c(
      "coefficients",
      "std.error",
      "statistic",
      "p.value",
      "df"
    )
  }

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

  if (!is.null(conf.level) && conf.int) {
    ci <- stats::confint(x, level = conf.level, ...)
    if (all(row.names(ci) == return_frame$term)) {
      return_frame$conf.low <- ci[, 1]
      return_frame$conf.high <- ci[, 2]
    }
  }
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

#' Tidy an estimatrZero object
#' @name estimatrZero_tidiers
#' @templateVar class lm_robust
#'
#' @param x An object returned by one of the estimators
#' @param conf.int Logical, whether to include confidence intervals.
#' @param conf.level The confidence level for intervals.
#' @param ... extra arguments (not used)
#'
#' @export
tidy.lm_robust <- function(x,
                           conf.int = TRUE,
                           conf.level = NULL,
                           ...) {
  warn_singularities(x)
  tidy_data_frame(x, conf.int = conf.int, conf.level = conf.level, ...)
}

#' @rdname estimatrZero_tidiers
#' @export
tidy.iv_robust <- function(x, conf.int = TRUE, conf.level = NULL, ...) {
  warn_singularities(x)
  tidy_data_frame(x, conf.int = conf.int, conf.level = conf.level, ...)
}

#' @rdname estimatrZero_tidiers
#' @export
tidy.difference_in_means <- tidy_data_frame

#' @rdname estimatrZero_tidiers
#' @export
tidy.horvitz_thompson <- tidy_data_frame

#' @rdname estimatrZero_tidiers
#' @export
tidy.lh_robust <- function(x,
                           conf.int = TRUE,
                           conf.level = NULL,
                           ...) {
  rbind(tidy(x$lm_robust, conf.int = conf.int, conf.level = conf.level, ...),
        tidy(x$lh, conf.int = conf.int, conf.level = conf.level, ...))
}

#' @rdname estimatrZero_tidiers
#' @export
tidy.lh <- function(x,
                    conf.int = TRUE,
                    conf.level = NULL,
                    ...) {
  tidy_data_frame(simplify_lh_outcome(x), conf.int = conf.int, conf.level = conf.level, ...)
}

simplify_lh_outcome <- function(x) {
  x_list <- as.list(x)
  x_list[["outcome"]] <- unique(x_list[["outcome"]])
  class(x_list) <- "lh"
  x_list
}
