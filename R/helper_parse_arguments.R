# This function parses condition names for HT and DiM estimators
parse_conditions <- function(treatment, condition1, condition2, estimator) {
  if (is.factor(treatment)) {
    condition_names <- levels(droplevels(treatment))
  } else {
    condition_names <- sort(unique(treatment))
  }

  if (any(!(c(condition1, condition2) %in% condition_names))) {
    stop("`condition1` and `condition2` must be values found in the treatment")
  }

  n_conditions <- length(condition_names)

  conditions <- list(NULL, NULL)

  if (n_conditions > 2) {
    if (is.null(condition1) || is.null(condition2)) {
      stop(
        "Treatment has > 2 values; must specify both `condition1` and ",
        "`condition2` or use a treatment with only 2 values"
      )
    } else {
      conditions[1:2] <- c(condition1, condition2)
    }
  } else if (n_conditions == 2) {
    if (is.null(condition1) && is.null(condition2)) {
      conditions[1:2] <- condition_names
    } else if (!is.null(condition2) && is.null(condition1)) {
      conditions[1:2] <- c(setdiff(condition_names, condition2), condition2)
    } else if (!is.null(condition1) && is.null(condition2)) {
      conditions[1:2] <- c(condition1, setdiff(condition_names, condition1))
    } else {
      conditions[1:2] <- c(condition1, condition2)
    }
  } else if (n_conditions == 1) {
    # Allowable for HT estimator
    if (estimator != "horvitz_thompson") {
      stop(
        "Must have more than one value in treatment unless using Horvitz-",
        "Thompson estimator"
      )
    }

    if (is.null(condition1) && is.null(condition2)) {
      conditions[2] <- condition_names
    } else if (!is.null(condition2)) {
      conditions[2] <- condition2
    } else if (!is.null(condition1)) {
      conditions[1] <- condition1
    }
  }

  return(conditions)
}

# This function ensures that blocks and  clusters have been specified correctly
check_clusters_blocks <- function(data) {
  if (!is.null(data$cluster)) {
    one_block_per_clust <-
      tapply(data$block, data$cluster, function(x) all(x == x[1]))

    # Check that clusters nest within blocks
    if (any(!one_block_per_clust)) {
      stop("All `clusters` must be contained within `blocks`")
    }

    # get number of clusters per block
    clust_per_block <- tapply(
      data$cluster,
      data$block,
      function(x) length(unique(x))
    )
  } else {
    clust_per_block <- tabulate(as.factor(data$block))
  }

  return(clust_per_block)
}
