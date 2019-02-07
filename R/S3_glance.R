#' @importFrom generics glance
#' @export
generics::glance

#' Glance at an estimatr object
#' @name estimatr_glancers
#' @templateVar class lm_robust
#' @return A data.frame with columns:
#'   \item{r.squared}{the \eqn{R^2},
#'   \deqn{R^2 = 1 - Sum(e[i]^2) / Sum((y[i] - y^*)^2),} where \eqn{y^*}
#'   is the mean of \eqn{y[i]} if there is an intercept and zero otherwise,
#'   and \eqn{e[i]} is the ith residual.}
#'   \item{adj.r.squared}{the \eqn{R^2} but penalized for having more parameters, \code{rank}}
#'   \item{se_type}{the standard error type specified by the user}
#'   \item{statistic}{the value of the F-statistic}
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
glance.lm_robust <- function(x, ...) {
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
      N = x[["N"]],
      se_type = x[["se_type"]]
    )
  )

  rownames(ret) <- NULL

  as.data.frame(ret)
}

#' Glance at an estimatr object
#' @name estimatr_glancers
#' @templateVar class iv_robust
#' @return A data.frame with columns:
#'   \item{r.squared}{The \eqn{R^2} of the second stage regression}
#'   \item{adj.r.squared}{The \eqn{R^2} but penalized for having more parameters, \code{rank}}
#'   \item{df.residual}{residual degrees of freedom}
#'   \item{N}{the number of observations used}
#'   \item{se_type}{the standard error type specified by the user}
#'   \item{statistic}{the value of the F-statistic}
#'   \item{p.value}{p-value from the F test}
#'   \item{statistic.weakinst}{the value of the first stage F-statistic, useful for the weak instruments test}
#'   \item{p.value.weakinst}{p-value from the first-stage F test, a test of weak instruments}
#'
#' @param x An object returned by one of the estimators
#' @param ... extra arguments (not used)
#'
#' @export
#' @family estimatr glancers
#' @seealso [generics::glance()], [estimatr::iv_robust()]
#' @md
glance.iv_robust <- function(x, ...) {
  ret <- cbind(
    data.frame(
      r.squared = x[["r.squared"]],
      adj.r.squared = x[["adj.r.squared"]],
      df.residual = x[["df.residual"]],
      N = x[["N"]],
      se_type = x[["se_type"]]
    ),
    if (exists("fstatistic", x)) {
      data.frame(
        statistic = x[["fstatistic"]][1],
        p.value = pf(x[["fstatistic"]][1], x[["fstatistic"]][2], x[["fstatistic"]][3], lower.tail = FALSE)
      )
    } else {
      data.frame(statistic = NA_real_, p.value = NA_real_)
    },
    if (exists("firststage_fstatistic", x)) {
      data.frame(
        statistic.weakinst = x[["firststage_fstatistic"]][1],
        p.value.weakinst = pf(x[["firststage_fstatistic"]][1], x[["firststage_fstatistic"]][2], x[["firststage_fstatistic"]][3], lower.tail = FALSE)
      )
    } else {
      data.frame(statistic.weakinst = NA_real_, p.value.weakinst = NA_real_)
    }
  )

  as.data.frame(ret)
}

#' Glance at an estimatr object
#' @name estimatr_glancers
#' @templateVar class difference_in_means
#' @return A data.frame with columns:
#'   \item{r.squared}{The \eqn{R^2} of the second stage regression}
#'   \item{adj.r.squared}{The \eqn{R^2} but penalized for having more parameters, \code{rank}}
#'   \item{df.residual}{residual degrees of freedom}
#'   \item{N}{the number of observations used}
#'   \item{statistic}{the value of the F-statistic}
#'   \item{p.value}{p-value from the F test}
#'   \item{statistic.weakinst}{the value of the first stage F-statistic, useful for the weak instruments test}
#'   \item{p.value.weakinst}{p-value from the first-stage F test, a test of weak instruments}
#'
#' @param x An object returned by one of the estimators
#' @param ... extra arguments (not used)
#'
#' @export
#' @family estimatr glancers
#' @seealso [generics::glance()], [estimatr::difference_in_means()]
#' @md
glance.difference_in_means <- function(x, ...) {
  ret <- cbind(
    data.frame(
      design = x[["design"]],
      df = x[["df"]],
      N = x[["N"]]
    ),
    if (exists("N_blocks", x)) {
      data.frame(N_blocks = x[["N_blocks"]])
    } else {
      data.frame(N_blocks = NA_real_)
    },
    if (exists("N_clusters", x)) {
      data.frame(N_clusters = x[["N_clusters"]])
    } else {
      data.frame(N_clusters = NA_real_)
    }
  )

  as.data.frame(ret)
}
