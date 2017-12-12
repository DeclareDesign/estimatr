# This file has helper functions for returning the lists from various estimators
lm_return <- function(return_list, model_data, formula) {

  return_list[["contrasts"]] <- attr(model_data$design_matrix, "contrasts")
  return_list[["terms"]] <- model_data$terms
  return_list[["weights"]] <- model_data$weights
  return_list[["outcome"]] <- all.vars(formula[[2]])

  return(return_list)

}

dim_like_return <- function(return_list, alpha, formula) {

  return_list[["alpha"]] <- alpha
  return_list[["coefficient_name"]] <- all.vars(formula[[3]])
  return_list[["outcome"]] <- all.vars(formula[[2]])

  return(return_list)

}
