# Internal method that creates linear fits
lm_fit <- function(y,
                   design_matrix,
                   weights,
                   cluster,
                   ci,
                   se_type,
                   alpha,
                   coefficient_name) {

  ## allowable se_types with clustering
  cl_se_types <- c("BM", "stata")

  ## Parse cluster variable
  if (!is.null(cluster)) {

    # set/check se_type
    if (is.null(se_type)) {
      se_type <- "BM"
    } else if (!(se_type %in% cl_se_types)) {
      stop("Only 'BM' or 'stata' allowed for se_type with clustered standard errors.")
    }

  } else {

    # set/check se_type
    if (is.null(se_type)) {
      se_type <- "HC2"
    } else if (se_type %in% cl_se_types) {
      stop("'BM' and 'stata' only allowed for clustered standard errors.")
    }

  }

  variable_names <- colnames(design_matrix)

  # Get coefficients to get df adjustments for and return
  if (is.null(coefficient_name)) {

    which_covs <- rep(TRUE, ncol(design_matrix))

  } else {

    # subset return to coefficients the user asked for
    which_covs <- variable_names %in% coefficient_name

    # if ever we can figure out all the use cases in the test....
    # which_ests <- return_frame$variable_names %in% deparse(substitute(coefficient_name))
  }

  if (!is.null(weights)) {
    design_matrix <- sqrt(weights) * design_matrix
    y <- sqrt(weights) * y
  }


  fit <-
    lm_robust_helper(
      y = y,
      X = design_matrix,
      cluster = cluster,
      ci = ci,
      type = se_type,
      which_covs = which_covs
    )

  est <- fit$beta_hat
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

        N <- nrow(design_matrix)
        k <- ncol(design_matrix)
        dof <- N - k

      }

      p <- 2 * pt(abs(est / se), df = dof, lower.tail = FALSE)
      ci_lower <- est - qt(1 - alpha / 2, df = dof) * se
      ci_upper <- est + qt(1 - alpha / 2, df = dof) * se

    }

  }

  return_frame <-
    data.frame(
      coefficient_name = variable_names,
      est = est,
      se = se,
      p = p,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      df = dof,
      stringsAsFactors = FALSE
    )

  rownames(return_frame) <- NULL

  return(return_frame[which_covs, ])

}
