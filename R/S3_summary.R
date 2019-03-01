#' @export
summary.lm_robust <- function(object, ...) {
  if (is.matrix(coef(object))) {
    ny <- ncol(coef(object))

    ret <- setNames(
      vector("list", ny),
      paste("Response", object$outcome)
    )

    mat_objs <- c(
      "coefficients",
      "std.error",
      "statistic",
      "df",
      "conf.low",
      "conf.high",
      "p.value"
    )

    vec_objs <- c(
      "outcome",
      "r.squared",
      "adj.r.squared",
      "res_var"
    )

    all_models <- object

    for (i in seq(ny)) {
      for (nm in names(object)) {
        if (nm %in% mat_objs) {
          object[[nm]] <- all_models[[nm]][, i, drop = TRUE]
        } else if (nm %in% vec_objs) {
          object[[nm]] <- all_models[[nm]][i]
        } else if (nm == "fstatistic") {
          object[[nm]] <- all_models[[nm]][c(i, ny + 1:2)]
        }
      }
      object$call$formula[[2L]] <- object$terms[[2L]] <-
        as.name(all_models$outcome[i])
      ret[[i]] <- summary(object, ...)
    }

    class(ret) <- "listof"
  } else {
    ret <- summary_lm_model(object)
  }

  ret
}

#' @export
summary.iv_robust <- function(object, ...) {
  summary_lm_model(object)
}


summary_lm_model <- function(object) {

  out_values <- c(
    "call",
    "k",
    "rank",
    "df.residual",
    "res_var",
    "weighted",
    "se_type",
    "fes",
    "r.squared",
    "adj.r.squared",
    "fstatistic"
  )
  # Different returns if fixed effects in the output
  if (object[["fes"]]) {
    out_values <- c(
      out_values,
      "proj_r.squared",
      "proj_adj.r.squared",
      "proj_fstatistic"
    )
  }

  # Different returns if fixed effects in the output
  if (is.numeric(object[["diagnostic_endogeneity_test"]])) {
    out_values <- c(
      out_values,
      "diagnostic_first_stage_fstatistic",
      "diagnostic_endogeneity_test",
      "diagnostic_overid_test"
    )
  }

  return_list <- object[out_values]

  # Split into two lists if multivariate linear model

  return_list[["coefficients"]] <- summarize_tidy(object)
  return_list[["N"]] <- nobs(object)

  class(return_list) <- "summary.lm_robust"
  return(return_list)
}


#' @export
summary.difference_in_means <- function(object, ...) {
  return(list(
    coefficients = summarize_tidy(object),
    design = object$design
  ))
}


#' @export
summary.horvitz_thompson <- function(object, ...) {
  return(list(coefficients = summarize_tidy(object, "z")))
}

summarize_tidy <- function(object, test = "t", ...) {
  remove_cols <- c("term", "outcome")

  # Ugly so that summary(fit)$coefficients matches lm()
  tidy_out <- tidy(object, ...)
  colnames(tidy_out)[2:8] <-
    c(
      "Estimate",
      "Std. Error",
      paste0(test, " value"),
      paste0("Pr(>|", test, "|)"),
      "CI Lower",
      "CI Upper",
      "DF"
    )
  tidy_mat <- as.matrix(tidy_out[, !(names(tidy_out) %in% remove_cols)])

  ny <- length(object$outcome)
  p <- length(object$term)
  if (length(object$outcome) > 1) {
    rownames(tidy_mat) <- paste0(
      rep(object$outcome, each = p),
      ":",
      rep(object$term, times = ny)
    )
  } else {
    rownames(tidy_mat) <- object$term
  }

  return(tidy_mat)
}



#' @export
summary.lh_robust <- function(object,...){
  class(object) <- "summary.lh_robust"
  object
}

#' @export
summary.lh <- function(object,...){
  attr(object, "linearHypothesis")
}

