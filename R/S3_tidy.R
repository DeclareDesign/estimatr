#' @importFrom generics tidy
#' @importFrom tibble as_tibble
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

  # confint() builds its rows from the same as.vector() of the same fields, so
  # they line up with these by position. They used to be matched by name
  # instead, and a multivariate fit's intervals are named "<outcome>:<term>",
  # so the names never matched and conf.level was ignored without a word,
  # returning the fit's own 95% intervals under a request for 90% (1.0.6 too).
  if (!is.null(conf.level) && conf.int) {
    ci <- stats::confint(x, level = conf.level, ...)
    return_frame$conf.low <- ci[, 1]
    return_frame$conf.high <- ci[, 2]
  }
  as_tibble(return_frame)
}

warn_singularities <- function(x) {
  if (x$rank < x$k) {
    singularities <- x$k - x$rank
    what <- if (singularities > 1) "coefficients" else "coefficient"
    message(
      singularities, " ", what,
      " not defined because the design matrix is rank deficient"
    )
  }
}

#' Tidy an estimatr object
#' @name estimatr_tidiers
#' @templateVar class lm_robust
#'
#' @param x An object returned by one of the estimators
#' @param conf.int Logical, whether to include confidence intervals.
#' @param conf.level The confidence level for intervals.
#' @param ... (optional) Ignored.
#'
#' @return A tibble with one row per term (and per outcome, for a multivariate
#'   fit): `term`, `estimate`, `std.error`, `statistic`, `p.value`,
#'   `conf.low`, `conf.high`, `df`, and `outcome`. A tibble rather than a
#'   plain data frame, as broom's tidiers return; 1.x returned a data frame.
#'
#' @examples
#' set.seed(50)
#' dat <- data.frame(x = rnorm(50), z = rep(0:1, 25))
#' dat$y <- dat$x + 0.4 * dat$z + rnorm(50)
#' fit <- lm_robust(y ~ x + z, data = dat)
#'
#' # One row per term, with the interval the fit was built with
#' tidy(fit)
#' tidy(fit, conf.int = FALSE)
#' tidy(fit, conf.level = 0.9)
#'
#' # The same shape for every estimator in the package
#' tidy(difference_in_means(y ~ z, data = dat))
#'
#' @export
tidy.lm_robust <- function(x,
                           conf.int = TRUE,
                           conf.level = NULL,
                           ...) {
  warn_singularities(x)
  tidy_data_frame(x, conf.int = conf.int, conf.level = conf.level, ...)
}

#' @rdname estimatr_tidiers
#' @export
tidy.iv_robust <- function(x, conf.int = TRUE, conf.level = NULL, ...) {
  warn_singularities(x)
  tidy_data_frame(x, conf.int = conf.int, conf.level = conf.level, ...)
}

#' @rdname estimatr_tidiers
#' @export
tidy.difference_in_means <- tidy_data_frame

#' @rdname estimatr_tidiers
#' @export
tidy.horvitz_thompson <- tidy_data_frame

#' @rdname estimatr_tidiers
#' @export
tidy.lh_robust <- function(x,
                           conf.int = TRUE,
                           conf.level = NULL,
                           ...) {
  as_tibble(rbind(
    tidy(x$lm_robust, conf.int = conf.int, conf.level = conf.level, ...),
    tidy(x$lh, conf.int = conf.int, conf.level = conf.level, ...)
  ))
}

#' @rdname estimatr_tidiers
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
