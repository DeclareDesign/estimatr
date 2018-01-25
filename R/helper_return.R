# This file has helper functions for returning the lists from various estimators
lm_return <- function(return_list, model_data, formula) {
  return_list[["contrasts"]] <- attr(model_data$design_matrix, "contrasts")
  return_list[["terms"]] <- model_data$terms
  return_list[["weights"]] <- model_data$weights
  return_list[["outcome"]] <- deparse(formula[[2]], nlines = 5)
  names(return_list[["coefficients"]]) <- return_list[["coefficient_name"]]

  return(return_list)
}

dim_like_return <- function(return_list, alpha, formula, conditions) {
  return_list[["alpha"]] <- alpha

  # get "max" condition to account for case with only 1 condition
  treat_condition <- conditions[[2]]

  # now we add the condition 2 value to coefficient name like it were a factor

  # Only add label if conditions aren't 0/1
  add_label <- !(conditions[[2]] == 1 && conditions[[1]] == 0)
  # If one of the two conditions is NULL, add_label will be NA
  # In this case, we always add label only if condntion2 is the nonnull
  # value
  if (is.na(add_label)) {
    add_label <- !is.null(conditions[[2]])
  }

  fterms <- terms(formula)
  coef_name <- labels(fterms)

  if (add_label) {
    return_list[["coefficient_name"]] <- paste0(
      coef_name,
      conditions[[2]]
    )
  } else {
    return_list[["coefficient_name"]] <- coef_name
  }

  return_list[["outcome"]] <- deparse(formula[[2]], nlines = 5)

  names(return_list[["coefficients"]]) <- return_list[["coefficient_name"]]

  return(return_list)
}
