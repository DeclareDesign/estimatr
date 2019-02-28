# This file has helper functions for returning the lists from various estimators
lm_return <- function(return_list, model_data, formula) {
  if (!is.null(model_data)) {
    return_list[["contrasts"]] <- attr(model_data$design_matrix, "contrasts")
    return_list[["terms"]] <- model_data$terms
    return_list[["xlevels"]] <- model_data$xlevels
    return_list[["felevels"]] <- model_data$felevels
    return_list[["weights"]] <- model_data$weights
    if (is.matrix(model_data$outcome) &&
        is.character(colnames(model_data$outcome))) {
      return_list[["outcome"]] <- colnames(model_data$outcome)
    } else {
      return_list[["outcome"]] <- deparse(formula[[2]], nlines = 5)
    }
  }

  # Name and flatten objects
  if (is.matrix(return_list[["std.error"]]) &&
      ncol(return_list[["std.error"]]) > 1) {
    dimnames(return_list[["std.error"]]) <- dimnames(return_list[["coefficients"]])
  } else {
    return_list[["coefficients"]] <- drop(return_list[["coefficients"]])
    nms <- c("std.error", "statistic", "p.value", "df", "conf.low", "conf.high")
    for (nm in nms) {
      if (length(return_list[[nm]]) > 1 || !is.na(return_list[[nm]])) {
        return_list[[nm]] <- setNames(
          drop(return_list[[nm]]),
          names(return_list[["coefficients"]])
        )
      }
    }
  }
  if (return_list[["weighted"]]) {
    names(return_list[["weights"]]) <- if (is.matrix(return_list[["fitted.values"]]))
      rownames(return_list[["fitted.values"]])
      else names(return_list[["fitted.values"]])
  }
  return_list[["fitted.values"]] <- drop(return_list[["fitted.values"]])
  return_list[["ei.iv"]] <- drop(return_list[["ei.iv"]])
  return_list[["residuals"]] <- drop(return_list[["residuals"]])
  return(return_list)
}

dim_like_return <- function(return_list, alpha, formula, conditions) {
  return_list[["alpha"]] <- alpha

  # get "max" condition to account for case with only 1 condition
  treat_condition <- conditions[[2]]

  # now we add the condition 2 value to coefficient name like it were a factor

  # Only add label if conditions aren't 0/1
  add_label <- !(conditions[[2]] == 1 && conditions[[1]] == 0)
  # If horvitz_thompson and there is only one treatment, add_label will be NA
  # In this case, we add the non-null value if it's condition 2
  if (is.na(add_label)) {
    add_label <- !is.null(conditions[[2]])
  }

  fterms <- terms(formula)
  coef_name <- labels(fterms)

  if (add_label) {
    return_list[["term"]] <- paste0(
      coef_name,
      conditions[[2]]
    )
  } else {
    return_list[["term"]] <- coef_name
  }

  return_list[["outcome"]] <- deparse(formula[[2]], nlines = 5)

  names(return_list[["coefficients"]]) <-
    names(return_list[["std.error"]]) <-
    names(return_list[["p.value"]]) <-
    names(return_list[["df"]]) <- return_list[["term"]]

  return_list[["condition2"]] <- conditions[[2]]
  return_list[["condition1"]] <- conditions[[1]]

  return(return_list)
}
