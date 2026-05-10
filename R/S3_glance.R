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
      df.residual = x[["df"]][1],
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
    if (exists("diagnostic_first_stage_fstatistic", x) && length(x[["diagnostic_first_stage_fstatistic"]] == 4)) {
      data.frame(
        statistic.weakinst = x[["diagnostic_first_stage_fstatistic"]]["value"],
        p.value.weakinst = x[["diagnostic_first_stage_fstatistic"]]["p.value"]
      )
    } else {
      data.frame(statistic.weakinst = NA_real_, p.value.weakinst = NA_real_)
    },
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
