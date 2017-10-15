#' Summarizing \code{\link{lm_robust}} objects
#'
#' @param object an object of class \code{\link{lm_robust}}
#' @param ... arguments passed to \code{\link{tidy.lm_robust}}, unused
#'
#' @export
summary.lm_robust <-
  function(
    object,
    ...
  ) {

    # This is ugly SO THAT summary(fit)$coefficients returns something like lm does.
    tidy_out <- tidy(object, ...)
    colnames(tidy_out)[2:4] <- c("Estimate", "Std. Error", "Pr(>|t|)")
    tidy_mat <- as.matrix(tidy_out[,-1])
    rownames(tidy_mat) <- tidy_out$coefficient_name
    return(list(coefficients = tidy_mat))

  }


#' Summarizing \code{\link{difference_in_means}} objects
#'
#' @param object an object of class \code{\link{difference_in_means}}
#' @param ... arguments passed to \code{\link{tidy.difference_in_means}}, unused
#'
#' @export
summary.difference_in_means <-
  function(
    object,
    ...
  ) {

    tidy(object, ...)

  }


#' Summarizing \code{\link{horvitz_thompson}} objects
#'
#' @param object an object of class \code{\link{horvitz_thompson}}
#' @param ... arguments passed to \code{\link{tidy.horvitz_thompson}}, unused
#'
#' @export
summary.horvitz_thompson <-
  function(
    object,
    ...
  ) {

    tidy(object, ...)

  }
