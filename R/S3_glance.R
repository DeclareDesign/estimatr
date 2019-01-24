#' @importFrom generics glance
#' @export
generics::glance

glance_lm_robust <- function(x, ...) {
  ret <- cbind(
    data.frame(
      r.squared = x[["r.squared"]],
      adj.r.squared = x[["adj.r.squared"]]
    ),
    if (exists("fstatistic", x)) {
      data.frame(
        statistic = x[["fstatistic"]][1],
        p.value = pf(x[["fstatistic"]][1], x[["fstatistic"]][2], x[["fstatistic"]][3], lower.tail = FALSE)
      )
    } else {
      data.frame(statistic = NA_real_, p.value = NA_real_)
    },
    data.frame(
      df.residual = x[["df"]][1],
      N = x[["N"]]
    )
  )

  rownames(ret) <- NULL

  as.data.frame(ret)
}

#' Glance at an estimatr object
#' @name estimatr_glancers
#' @templateVar class lm_robust
#' @return A data.frame with columns:
#'   \item{r.squared}{The \eqn{R^2},
#'   \deqn{R^2 = 1 - Sum(e[i]^2) / Sum((y[i] - y^*)^2),} where \eqn{y^*}
#'   is the mean of \eqn{y[i]} if there is an intercept and zero otherwise,
#'   and \eqn{e[i]} is the ith residual.}
#'   \item{adj.r.squared}{The \eqn{R^2} but penalized for having more parameters, \code{rank}}
#'   \item{statistic}{a vector with the value of the F-statistic with the numerator and denominator degrees of freedom}
#'   \item{p.value}{p-value from the F test}
#'   \item{df.residual}{residual degrees of freedom}
#'   \item{N}{the number of observations used}
#'
#' @param x An object returned by one of the estimators
#' @param ... extra arguments (not used)
#'
#' @export
#' @family estimatr glancers
#' @seealso [generics::glance()], [estimatr::lm_robust()], [estimatr::lm_lin()]
#' @md
glance.lm_robust <- glance_lm_robust
