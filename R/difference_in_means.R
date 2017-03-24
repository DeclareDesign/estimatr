#' Built-in Estimators: Difference-in-means
#'
#' @param formula An object of class "formula", such as Y ~ Z
#' @param data A data.frame, often created by \code{\link{draw_population}}.
#' @param weights An optional vector of weights (not yet implemented).
#' @param subset An optional vector specifying a subset of observations to be used.
#' @param alpha The significance level, 0.05 by default.
#'
#' @export
difference_in_means <-
  function(formula,
           condition1 = NULL,
           condition2 = NULL,
           data,
           weights = NULL,
           subset = NULL,
           alpha = .05) {
    if (length(all.vars(formula[[3]])) > 1)
      stop(
        "The formula should only include one variable on the right-hand side: the treatment variable."
      )

    if (!is.null(subset)){
      condition_call <- substitute(condition)
      r <- eval(condition_call, data)
      data <- data[r, ]
    }

    Y <- data[, all.vars(formula[[2]])]
    t <- data[, all.vars(formula[[3]])]

    if (is.factor(t)) {
      condition_names <- levels(t)
    } else{
      condition_names <- sort(unique(t))
    }

    if (is.null(condition1) & is.null(condition2)) {
      condition1 <- condition_names[1]
      condition2 <- condition_names[2]
    }

    N <- length(Y)

    Y2 <- Y[t == condition2]
    Y1 <- Y[t == condition1]

    weights <- eval(substitute(weights), data)

    if (is.null(weights)) {
      diff <- mean(Y2) - mean(Y1)

      se <- sqrt(var(Y2) / length(Y2) + var(Y1) / length(Y1))

    } else {

      w2 <- weights[t == condition2]
      w1 <- weights[t == condition1]

      mean2 <- weighted.mean(Y2, w2)
      mean1 <- weighted.mean(Y1, w1)
      diff <-  mean2 - mean1

      se <- sqrt(weighted_var_internal(w2, Y2, mean2) + weighted_var_internal(w1, Y1, mean1))
    }

    df <- N - 2
    p <- 2 * pt(abs(diff / se), df = df, lower.tail = FALSE)
    ci_lower <- diff - qt(1 - alpha / 2, df = df) * se
    ci_upper <- diff + qt(1 - alpha / 2, df = df) * se

    return_df <- data.frame(
      est = diff,
      se = se,
      p = p,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      df = df
    )

    return(return_df)

  }



weighted_var_internal <- function(w, x, xWbar){
  wbar <- mean(w)
  n <- length(w)
  return(n / ((n - 1) * sum(w) ^ 2) * (sum((w * x - wbar * xWbar) ^ 2) -
                                         2 * xWbar * sum((w - wbar) * (w * x - wbar * xWbar)) + xWbar ^ 2 * sum((w -
                                                                                                                   wbar) ^ 2)))
}

#' Built-in Estimators: Blocked Difference-in-means
#'
#' @param formula An object of class "formula", such as Y ~ Z
#' @param data A data.frame, often created by \code{\link{draw_population}}.
#' @param weights An optional vector of weights (not yet implemented).
#' @param subset An optional vector specifying a subset of observations to be used.
#' @param block_variable_name The name of the blocking variable.
#' @param cluster_variable_name The (optional) name of a clustered variable. The function will first collapse the data by cluster before calculating the difference-in-means.
#' @param cluster_collapse_function The (optional) function to use to collapse the data by cluster. Default is \code{mean}.
#' @param alpha The significance level, 0.05 by default.
#'
#' @export
difference_in_means_blocked <-
  function(formula,
           block_variable_name = NULL,
           condition1 = NULL,
           condition2 = NULL,
           data,
           weights = NULL,
           subset = NULL,
           alpha = .05) {

    if (!is.null(subset)){
      condition_call <- substitute(condition)
      r <- eval(condition_call, data)
      data <- data[r, ]
    }

    blocks <- eval(substitute(block_variable_name), data)

    block_estimates <-
      data %>%
      split(blocks) %>%
      map(~difference_in_means_internal(formula, data = ., weights = weights, alpha = alpha)) %>%
      bind_rows()

    N_overall <- with(block_estimates, sum(N))
    diff <- with(block_estimates, sum(est * N/N_overall))
    se <- with(block_estimates, sqrt(sum(se^2 * (N/N_overall)^2)))

    ## we don't know if this is correct!
    df <- N_overall - 2
    p <- 2 * pt(abs(diff / se), df = df, lower.tail = FALSE)
    ci_lower <- diff - qt(1 - alpha / 2, df = df) * se
    ci_upper <- diff + qt(1 - alpha / 2, df = df) * se

    return_df <- data.frame(
      est = diff,
      se = se,
      p = p,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      df = df
    )

    return(return_df)

  }


difference_in_means_internal <-
  function(formula,
           condition1 = NULL,
           condition2 = NULL,
           data,
           weights = NULL,
           alpha = .05) {

    Y <- data[, all.vars(formula[[2]])]
    t <- data[, all.vars(formula[[3]])]

    if (is.factor(t)) {
      condition_names <- levels(t)
    } else{
      condition_names <- sort(unique(t))
    }

    if (is.null(condition1) & is.null(condition2)) {
      condition1 <- condition_names[1]
      condition2 <- condition_names[2]
    }

    N <- length(Y)

    Y2 <- Y[t == condition2]
    Y1 <- Y[t == condition1]

    weights <- eval(substitute(weights), data)

    if (is.null(weights)) {
      diff <- mean(Y2) - mean(Y1)

      se <- sqrt(var(Y2) / length(Y2) + var(Y1) / length(Y1))

    } else {

      w2 <- weights[t == condition2]
      w1 <- weights[t == condition1]

      mean2 <- weighted.mean(Y2, w2)
      mean1 <- weighted.mean(Y1, w1)
      diff <-  mean2 - mean1

      se <- sqrt(weighted_var_internal(w2, Y2, mean2) + weighted_var_internal(w1, Y1, mean1))
    }

    return_df <- data.frame(
      est = diff,
      se = se,
      N = N
    )

    return(return_df)

  }
