#' @importFrom generics glance
#' @export
generics::glance

retrieve_value <- function(x, what) if (exists(what, x)) x[[what]] else NA_real_

retrieve_fstatistic <- function(x) {
  if (exists("fstatistic", x)) {
    data.frame(
      statistic = x[["fstatistic"]][1],
      p.value = pf(
        x[["fstatistic"]][1],
        x[["fstatistic"]][2],
        x[["fstatistic"]][3],
        lower.tail = FALSE
      )
    )
  } else {
    data.frame(statistic = NA_real_, p.value = NA_real_)
  }
}

#' @export
glance.lm_robust <- function(x, ...) {

  if (length(x[["outcome"]]) > 1) {
    stop("Cannot use `glance` on linear models with multiple responses.")
  }

  ret <- cbind(
    data.frame(
      r.squared = x[["r.squared"]],
      adj.r.squared = x[["adj.r.squared"]]
    ),
    retrieve_fstatistic(x),
    data.frame(
      # x[["df"]] is the per-coefficient degrees of freedom, which under CR2
      # is Satterthwaite and is not the residual df at all: on a 10-cluster
      # fit it read 8.56 where the residual df is 98. glance.iv_robust() has
      # always used df.residual; both do now.
      df.residual = x[["df.residual"]],
      nobs = as.integer(x[["nobs"]]),
      se_type = x[["se_type"]],
      stringsAsFactors = FALSE
    )
  )

  rownames(ret) <- NULL
  ret
}

#' @export
glance.lh_robust <- function(x, ...) {
  glance(x[["lm_robust"]])
}

# The first-stage F test has one entry per endogenous regressor, named
# "<var>:value" once there is more than one, so indexing by "value" returns an
# NA with an NA name and the data.frame call fails on the row name (estimatr
# #389). glance() must be one row, so with several endogenous regressors we
# report the weakest first stage, which is the quantity a weak-instrument
# diagnostic is asking about. Read the per-regressor tests in
# `diagnostic_first_stage_fstatistic` directly.
weakinst_row <- function(fstat) {
  if (is.null(fstat)) {
    return(data.frame(statistic.weakinst = NA_real_, p.value.weakinst = NA_real_))
  }
  stat_i <- grep("(^|:)value$", names(fstat))
  p_i <- grep("(^|:)p\\.value$", names(fstat))
  weakest <- which.min(fstat[stat_i])
  data.frame(
    statistic.weakinst = unname(fstat[stat_i][weakest]),
    p.value.weakinst = unname(fstat[p_i][weakest])
  )
}

#' @export
glance.iv_robust <- function(x, ...) {

  if (length(x[["outcome"]]) > 1) {
    stop("Cannot use `glance` on linear models with multiple responses.")
  }

  ret <- cbind(
    data.frame(
      r.squared = x[["r.squared"]],
      adj.r.squared = x[["adj.r.squared"]],
      df.residual = x[["df.residual"]],
      nobs = as.integer(x[["nobs"]]),
      se_type = x[["se_type"]],
      stringsAsFactors = FALSE
    ),
    retrieve_fstatistic(x),
    weakinst_row(x[["diagnostic_first_stage_fstatistic"]]),
    if (exists("diagnostic_endogeneity_test", x)) {
      data.frame(
        statistic.endogeneity = x[["diagnostic_endogeneity_test"]]["value"],
        p.value.endogeneity = x[["diagnostic_endogeneity_test"]]["p.value"]
      )
    } else {
      data.frame(statistic.endogeneity = NA_real_, p.value.endogeneity = NA_real_)
    },
    if (exists("diagnostic_overid_test", x)) {
      data.frame(
        statistic.overid = x[["diagnostic_overid_test"]]["value"],
        p.value.overid = x[["diagnostic_overid_test"]]["p.value"]
      )
    } else {
      data.frame(statistic.overid = NA_real_, p.value.overid = NA_real_)
    }
  )

  ret
}

#' @export
glance.difference_in_means <- function(x, ...) {
  data.frame(
    design = x[["design"]],
    df = x[["df"]],
    nobs = as.integer(x[["nobs"]]),
    nblocks = retrieve_value(x, "nblocks"),
    nclusters = retrieve_value(x, "nclusters"),
    condition2 = x[["condition2"]],
    condition1 = x[["condition1"]],
    stringsAsFactors = FALSE
  )
}

#' @export
glance.horvitz_thompson <- function(x, ...) {
  # Same four columns estimatr returns, so a table built over a mix of old and
  # new fits binds rather than erroring. Without this method `modelsummary()`
  # does not fail: it drops the goodness-of-fit rows and prints a coefficient
  # table that looks complete, which is the worse outcome of the two.
  data.frame(
    nobs = as.integer(x[["nobs"]]),
    se_type = x[["se_type"]],
    condition2 = x[["condition2"]],
    condition1 = x[["condition1"]],
    stringsAsFactors = FALSE
  )
}
