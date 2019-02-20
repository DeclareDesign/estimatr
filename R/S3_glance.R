#' @importFrom generics glance
#' @export
generics::glance

#' Glance at an estimatr object
#' @name estimatr_glancers
#' @templateVar class lm_robust
#' @return For \code{glance.lm_robust}, a data.frame with columns:
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
#' @seealso [generics::glance()], [estimatr::lm_robust()], [estimatr::lm_lin()], [estimatr::iv_robust()], [estimatr::difference_in_means()], [estimatr::horvitz_thompson()]
#' @md
glance.lm_robust <- function(x, ...) {

  if (length(x[["outcome"]]) > 1) {
    stop("Cannot use `glance` on linear models with multiple responses.")
  }

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

#' @name estimatr_glancers
#' @templateVar class iv_robust
#' @return For \code{glance.iv_robust}, a data.frame with columns:
#'   \item{r.squared}{The \eqn{R^2} of the second stage regression}
#'   \item{adj.r.squared}{The \eqn{R^2} but penalized for having more parameters, \code{rank}}
#'   \item{df.residual}{residual degrees of freedom}
#'   \item{N}{the number of observations used}
#'   \item{se_type}{the standard error type specified by the user}
#'   \item{statistic}{the value of the F-statistic}
#'   \item{p.value}{p-value from the F test}
#'   \item{statistic.weakinst}{the value of the first stage F-statistic, useful for the weak instruments test; only reported if there is only one endogenous variable}
#'   \item{p.value.weakinst}{p-value from the first-stage F test, a test of weak instruments; only reported if there is only one endogenous variable}
#'   \item{statistic.endogeneity}{the value of the F-statistic for the test of endogeneity; often called the Wu-Hausman statistic, with robust standard errors, we employ the regression based test}
#'   \item{p.value.endogeneity}{p-value from the F-test for endogeneity}
#'   \item{statistic.overid}{the value of the chi-squared statistic for the test of instrument correlation with the error term; only reported with overidentification}
#'   \item{p.value.overid}{p-value from the chi-squared test; only reported with overidentification}
#'
#' @inheritParams glance.lm_robust
#'
#' @export
#' @family estimatr glancers
#' @md
glance.iv_robust <- function(x, ...) {

  if (length(x[["outcome"]]) > 1) {
    stop("Cannot use `glance` on linear models with multiple responses.")
  }

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
    if (exists("diagnostic_firststage_fstatistic", x) && length(x[["diagnostic_firststage_fstatistic"]] == 4)) {
      data.frame(
        statistic.weakinst = x[["diagnostic_firststage_fstatistic"]]["value"],
        p.value.weakinst = x[["diagnostic_firststage_fstatistic"]]["p.value"]
      )
    } else {
      data.frame(statistic.weakinst = NA_real_, p.value.weakinst = NA_real_)
    },
    if (exists("diagnostic_endogeneity_fstatistic", x)) {
      data.frame(
        statistic.endogeneity = x[["diagnostic_endogeneity_fstatistic"]]["value"],
        p.value.endogeneity = x[["diagnostic_endogeneity_fstatistic"]]["p.value"]
      )
    } else {
      data.frame(statistic.endogeneity = NA_real_, p.value.endogeneity = NA_real_)
    },
    if (exists("diagnostic_overid_fstatistic", x)) {
      data.frame(
        statistic.overid = x[["diagnostic_overid_fstatistic"]]["value"],
        p.value.overid = x[["diagnostic_overid_fstatistic"]]["p.value"]
      )
    } else {
      data.frame(statistic.overid = NA_real_, p.value.overid = NA_real_)
    }
  )

  as.data.frame(ret)
}

#' @name estimatr_glancers
#' @templateVar class difference_in_means
#' @return For \code{glance.difference_in_means}, a data.frame with columns:
#'   \item{design}{the design used, and therefore the estimator used}
#'   \item{df}{the degrees of freedom}
#'   \item{N}{the number of observations used}
#'   \item{N_blocks}{the number of blocks, if used}
#'   \item{N_clusters}{the number of clusters, if used}
#'   \item{condition2}{the second, "treatment", condition}
#'   \item{condition1}{the first, "control", condition}
#'
#' @inheritParams glance.lm_robust
#'
#' @export
#' @family estimatr glancers
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
    },
    data.frame(
      condition2 = x[["condition2"]],
      condition1 = x[["condition1"]]
    )
  )

  as.data.frame(ret)
}

#' @name estimatr_glancers
#' @templateVar class horvitz_thompson
#' @return For \code{glance.horvitz_thompson}, a data.frame with columns:
#'   \item{N}{the number of observations used}
#'   \item{se_type}{the type of standard error estimator used}
#'   \item{condition2}{the second, "treatment", condition}
#'   \item{condition1}{the first, "control", condition}
#'
#' @inheritParams glance.lm_robust
#'
#' @export
#' @family estimatr glancers
#' @md
glance.horvitz_thompson <- function(x, ...) {
  ret <- data.frame(
    N = x[["N"]],
    se_type = x[["se_type"]],
    condition2 = x[["condition2"]],
    condition1 = x[["condition1"]]
  )

  as.data.frame(ret)
}
