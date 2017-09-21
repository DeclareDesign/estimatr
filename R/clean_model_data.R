# Internal method to check for missingness on auxiliary variables and warn
find_warn_missing <- function(x, type) {
  x_missing <- is.na(x)

  if (any(x_missing)) {
    warning(
      sprintf(
        "Some observations have missingness in the %s variable. These observations have been dropped.",
        type
      )
    )
  }

  return(x_missing)
}


# Internal method to process data
clean_model_data <- function(formula,
                             data,
                             weights,
                             condition_call,
                             cluster_variable_name) {

  if (!is.null(condition_call)) {
    r <- eval(condition_call, data)
    data <- data[r,]
  }

  mf <- model.frame.default(formula, data = data)

  mf_rows_to_drop <- list(cluster = integer(0),
                          weights = integer(0))
  ## Parse cluster variable
  if (!is.null(cluster_variable_name)) {
    # get cluster variable from subset of data
    cluster <- as.factor(eval(cluster_variable_name,
                              data[row.names(mf), ]))

    mf_rows_to_drop$cluster <- which(is.na(cluster))

    if (length(mf_rows_to_drop$cluster) != 0) {
      warning(
        "Some observations have missingness in the cluster variable but have not in the outcome or covariates. These observations have been dropped."
      )
    }

  } else {
    cluster <- NULL
  }

  if (!is.null(weights)) {

    if (!is.null(cluster)) {
      stop("weights not yet supported with clustered standard errors")
    }

    weights <- eval(weights, data[row.names(mf), ])

    mf_rows_to_drop$weights <- which(is.na(weights))

    if (length(mf_rows_to_drop$weights) != 0) {
      warning(
        "Some observations have missingness in the weights variable but not in the outcome or covariates. These observations have been dropped."
      )
    }

  } else {
    weights <- NULL
  }

  all_rows_to_drop <- unlist(mf_rows_to_drop, F, F)

  if (length(all_rows_to_drop) != 0) {
    mf <- mf[-all_rows_to_drop, ]
    cluster <- cluster[-all_rows_to_drop]
    weights <- weights[-all_rows_to_drop]
  }

  outcome <- model.response(mf)
  design_matrix <- model.matrix.default(formula, data = mf)

  return(list(outcome=outcome,
              design_matrix=design_matrix,
              cluster=cluster,
              weights=weights))

}
