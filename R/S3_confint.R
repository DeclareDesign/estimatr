confint_lm_like <- function(object,
                            parm = NULL,
                            level = NULL,
                            ...) {
  cis <- get_ci_mat(object, level)

  if (!is.null(parm)) {
    cis <- cis[parm, , drop = FALSE]
  }

  return(cis)
}

#' @export
confint.lm_robust <- confint_lm_like

#' @export
confint.iv_robust <- confint_lm_like

#' @importFrom stats confint
#' @export
confint.lh_robust <- function(object, parm = NULL, level = NULL, ...) {
  rbind(confint(object$lm_robust, parm = parm, level = level, ...), tidy(object$lh, ...))
}

#' @export
confint.lh <- confint_lm_like

#' @export
confint.difference_in_means <- function(object,
                                        parm = NULL,
                                        level = NULL,
                                        ...) {
  cis <- get_ci_mat(object, level)
  return(cis)
}

#' @export
confint.horvitz_thompson <- function(object,
                                     parm = NULL,
                                     level = NULL,
                                     ...) {
  # `ttest = FALSE` because a Horvitz-Thompson fit has no degrees of freedom to
  # spend: `df` is NA and the interval is normal. Passing the default would
  # rebuild the interval off a t quantile with df = NA and return NA bounds
  # whenever `level` is supplied.
  cis <- get_ci_mat(object, level, ttest = FALSE)
  return(cis)
}

get_ci_mat <- function(object, level, ttest = TRUE) {
  if (!is.null(level)) {
    if (!is.null(object[["alpha"]])) {
      object[["alpha"]] <- NULL
    }
    object <- add_cis_pvals(object, alpha = 1 - level, ci = TRUE, ttest = ttest)
  } else {
    level <- 1 - object$alpha
  }

  cis <- cbind(
    as.vector(object$conf.low),
    as.vector(object$conf.high)
  )

  if (is.matrix(object$conf.low)) {
    ny <- ncol(object$conf.low)
    p <- nrow(object$conf.low)
    rownames(cis) <- paste0(
      rep(object$outcome, each = p),
      ":",
      rep(object$term, times = ny)
    )
  } else {
    rownames(cis) <- object$term
  }

  colnames(cis) <- paste((1 - level) / 2 * c(100, -100) + c(0, 100), "%")

  return(cis)
}
