#' Internal method that creates linear fits
#'
#' @param y numeric outcome vector
#' @param X numeric design matrix
#' @param weights numeric weights vector
#' @param cluster numeric cluster vector
#' @param ci boolean that when T returns confidence intervals and p-values
#' @param se_type character denoting which kind of SEs to return
#' @param alpha numeric denoting the test size for confidence intervals
#' @param coefficient_name character vector of coefficients to return
#' @param return_vcov a boolean for whether to return the vcov matrix for later usage
#'
#' @export
#'
lm_robust_fit <- function(y,
                          X,
                          weights,
                          cluster,
                          ci,
                          se_type,
                          alpha,
                          coefficient_name,
                          return_vcov) {

  ## allowable se_types with clustering
  cl_se_types <- c("BM", "stata")
  rob_se_types <- c("HC0", "HC1", "HC2", "HC3", "classical")
  ## Parse cluster variable
  if (!is.null(cluster)) {

    # set/check se_type
    if (is.null(se_type)) {
      se_type <- "BM"
    } else if (!(se_type %in% c(cl_se_types, "none"))) {
      stop("Incorrect se_type. Only 'BM' or 'stata' allowed for se_type with clustered standard errors. Also can choose 'none'.")
    }

  } else {

    # set/check se_type
    if (is.null(se_type)) {
      se_type <- "HC2"
    } else if (se_type %in% cl_se_types) {
      stop("Incorrect se_type. 'BM' and 'stata' only allowed for clustered standard errors.")
    } else if (!(se_type %in% c(rob_se_types, "none"))) {
      stop('Incorrect se_type. "HC0", "HC1", "HC2", "HC3", "classical" are the se_type options without clustering. Also can choose "none".')
    }

  }

  variable_names <- colnames(X)

  # Get coefficients to get df adjustments for and return
  if (is.null(coefficient_name)) {

    which_covs <- rep(TRUE, ncol(X))

  } else {

    # subset return to coefficients the user asked for
    which_covs <- variable_names %in% coefficient_name

    # if ever we can figure out all the use cases in the test....
    # which_ests <- return_frame$variable_names %in% deparse(substitute(coefficient_name))
  }

  if (!is.null(weights)) {
    X <- sqrt(weights) * X
    y <- sqrt(weights) * y
  }

  fit <-
    lm_robust_helper(
      y = y,
      X = X,
      cluster = cluster,
      ci = ci,
      type = se_type,
      which_covs = which_covs
    )

  est <- as.vector(fit$beta_hat)
  se <- NA
  p <- NA
  ci_lower <- NA
  ci_upper <- NA
  dof <- NA

  if(se_type != "none"){

    se <- sqrt(diag(fit$Vcov_hat))

    if(ci) {

      if(se_type %in% cl_se_types){

        ## Replace -99 with NA, easy way to flag that we didn't compute
        ## the DoF because the user didn't ask for it
        dof <- ifelse(fit$dof == -99,
                      NA,
                      fit$dof)

      } else {

        N <- nrow(X)
        k <- ncol(X)
        dof <- N - k

      }

      dof <- as.vector(dof)

      p <- 2 * pt(abs(est / se), df = dof, lower.tail = FALSE)
      ci_lower <- est - qt(1 - alpha / 2, df = dof) * se
      ci_upper <- est + qt(1 - alpha / 2, df = dof) * se

    }

  }

  return_list <-
    list(
      coefficient_name = variable_names,
      est = est,
      se = se,
      p = p,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      df = dof
    )

  if (return_vcov) {
    return_list$vcov <- fit$Vcov_hat
    dimnames(return_list$vcov) <- list(return_list$coefficient_name,
                                       return_list$coefficient_name)
  }

  attr(return_list, "class") <- "lm_robust"

  return(return_list)

}
