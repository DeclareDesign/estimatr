

#' Ordinary Least Squares with Robust Standard Errors
#'
#' @param formula an object of class formula, as in \code{\link{lm}}.
#'
#' @param data A data.frame.
#' @param weights the bare (unquoted) names of the weights variable in the supplied data.
#' @param subset An optional bare (unquoted) expression specifying a subset of observations to be used.
#' @param cluster_variable_name An optional bare (unquoted) name of the factor variable that corresponds to the clusters in the data. Will return Bell-McCaffrey standard errors, overriding \code{se_type}.
#' @param se_type The sort of standard error sought. Without clustering: "HCO", "HC1", "HC2" (default), "HC3", or "classical". With clustering: "BM" (default).
#' @param alpha The significance level, 0.05 by default.
#' @param coefficient_name a character or character vector that indicates which coefficients should be reported. Defaults to "Z".
#'
#' @export
#'
lm_robust_se <- function(formula,
                         data,
                         weights = NULL,
                         subset = NULL,
                         cluster_variable_name = NULL,
                         se_type = NULL,
                         alpha = .05,
                         coefficient_name = "Z") {

  condition_call <- substitute(subset)

  if (is.null(se_type)) {
    if (is.null(cluster_variable_name)) {

    }
  }

  if (!is.null(condition_call)) {
    r <- eval(condition_call, data)
    data <- data[r,]
  }

  design_matrix <- model.matrix.default(formula, data = data)
  variable_names <- colnames(design_matrix)
  outcome <- data[, all.vars(formula[[2]])]


  if (!is.null(substitute(cluster_variable_name))) {

    # get cluster variable from subset of data
    cluster <- as.factor(eval(substitute(cluster_variable_name), data))

    # set/check se_type
    if (is.null(se_type)) {
      se_type <- "BM"
    } else if (se_type != "BM") {
      stop("Only se_type = 'BM' is allowed for clustered standard errors.")
    }

  } else {
    cluster <- NULL

    # set/check se_type
    if (is.null(se_type)) {
      se_type <- "HC2"
    } else if (se_type == "BM") {
      stop("se_type = 'BM' is only allowed for clustered standard errors.")
    }

  }

  if (!is.null(substitute(weights))) {
    weights <- eval(substitute(weights), data)
    outcome <- sqrt(weights) * outcome
    design_matrix <- sqrt(weights) * design_matrix
  }


  fit <-
    lm_robust_helper(
      y = outcome,
      X = design_matrix,
      cluster = cluster,
      type = se_type)

  est <- fit$beta_hat

  if(se_type == "none"){

    se = NA
    p = NA
    ci_lower = NA
    ci_upper = NA

  } else {

    if(se_type == "BM"){
      df <- fit$df
    } else {
      N <- nrow(design_matrix)
      k <- ncol(design_matrix)
      df <- N - k
    }

    se <- sqrt(diag(fit$Vcov_hat))
    p <- 2 * pt(abs(est / se), df = df, lower.tail = FALSE)
    ci_lower <- est - qt(1 - alpha / 2, df = df) * se
    ci_upper <- est + qt(1 - alpha / 2, df = df) * se

  }

  return_frame <-
    data.frame(
      variable_names = variable_names,
      est = est,
      se = se,
      p = p,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      stringsAsFactors = FALSE
    )

  which_ests <- return_frame$variable_names %in% coefficient_name

  # if ever we can figure out all the use cases in the test....
  # which_ests <- return_frame$variable_names %in% deparse(substitute(coefficient_name))

  return(return_frame[which_ests, ])

}
