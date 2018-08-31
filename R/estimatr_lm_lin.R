#' Linear regression with the Lin (2013) covariate adjustment
#'
#' @description This function is a wrapper for \code{\link{lm_robust}} that
#' is useful for estimating treatment effects with pre-treatment covariate
#' data. This implements the method described by Lin (2013).
#'
#' @param formula an object of class formula, as in \code{\link{lm}}, such as
#' \code{Y ~ Z} with only one variable on the right-hand side, the treatment
#' @param covariates a right-sided formula with pre-treatment covariates on
#' the right hand side, such as \code{ ~ x1 + x2 + x3}.
#' @param data A \code{data.frame}
#' @param weights the bare (unquoted) names of the weights variable in the
#' supplied data.
#' @param subset An optional bare (unquoted) expression specifying a subset
#' of observations to be used.
#' @param clusters An optional bare (unquoted) name of the variable that
#' corresponds to the clusters in the data.
#' @param se_type The sort of standard error sought. If `clusters` is
#' not specified the options are "HC0", "HC1" (or "stata", the equivalent),
#' "HC2" (default), "HC3", or  "classical". If `clusters` is specified the
#' options are "CR0", "CR2" (default), or "stata" are
#' permissible.
#' @param ci logical. Whether to compute and return p-values and confidence
#' intervals, TRUE by default.
#' @param alpha The significance level, 0.05 by default.
#' @param return_vcov logical. Whether to return the variance-covariance
#' matrix for later usage, TRUE by default.
#' @param try_cholesky logical. Whether to try using a Cholesky
#' decomposition to solve least squares instead of a QR decomposition,
#' FALSE by default. Using a Cholesky decomposition may result in speed gains, but should only
#' be used if users are sure their model is full-rank (i.e., there is no
#' perfect multi-collinearity)
#'
#' @details
#'
#' This function is simply a wrapper for \code{\link{lm_robust}} and implements
#' the Lin estimator (see the reference below). This method
#' pre-processes the data by taking the covariates specified in the
#' \code{`covariates`} argument, centering them by subtracting from each covariate
#' its mean, and interacting them with the treatment. If the treatment has
#' multiple values, a series of dummies for each value is created and each of
#' those is interacted with the demeaned covariates. More details can be found
#' in the
#' \href{https://declaredesign.org/r/estimatr/articles/getting-started.html}{Getting Started vignette}
#' and the
#' \href{https://declaredesign.org/r/estimatr/articles/mathematical-notes.html}{mathematical notes}.
#'
#' @return An object of class \code{"lm_robust"}.
#'
#' The post-estimation commands functions \code{summary} and \code{\link{tidy}}
#' return results in a \code{data.frame}. To get useful data out of the return,
#' you can use these data frames, you can use the resulting list directly, or
#' you can use the generic accessor functions \code{coef}, \code{vcov},
#' \code{confint}, and \code{predict}. Marginal effects and uncertainty about
#' them can be gotten by passing this object to
#' \code{\link[margins]{margins}} from the \pkg{margins}.
#'
#' Users who want to print the results in TeX of HTML can use the
#' \code{\link[texreg]{extract}} function and the \pkg{texreg} package.
#'
#' An object of class \code{"lm_robust"} is a list containing at least the
#' following components:
#'   \item{coefficients}{the estimated coefficients}
#'   \item{std.error}{the estimated standard errors}
#'   \item{statistic}{the t-statistic}
#'   \item{df}{the estimated degrees of freedom}
#'   \item{p.value}{the p-values from a two-sided t-test using \code{coefficients}, \code{std.error}, and \code{df}}
#'   \item{conf.low}{the lower bound of the \code{1 - alpha} percent confidence interval}
#'   \item{conf.high}{the upper bound of the \code{1 - alpha} percent confidence interval}
#'   \item{term}{a character vector of coefficient names}
#'   \item{alpha}{the significance level specified by the user}
#'   \item{se_type}{the standard error type specified by the user}
#'   \item{res_var}{the residual variance}
#'   \item{N}{the number of observations used}
#'   \item{k}{the number of columns in the design matrix (includes linearly dependent columns!)}
#'   \item{rank}{the rank of the fitted model}
#'   \item{vcov}{the fitted variance covariance matrix}
#'   \item{r.squared}{The \eqn{R^2},
#'   \deqn{R^2 = 1 - Sum(e[i]^2) / Sum((y[i] - y^*)^2),} where \eqn{y^*}
#'   is the mean of \eqn{y[i]} if there is an intercept and zero otherwise,
#'   and \eqn{e[i]} is the ith residual.}
#'   \item{adj.r.squared}{The \eqn{R^2} but penalized for having more parameters, \code{rank}}
#'   \item{weighted}{whether or not weights were applied}
#'   \item{call}{the original function call}
#'   \item{fitted.values}{the matrix of predicted means}
#' We also return \code{terms} and \code{contrasts}, used by \code{predict},
#' and \code{scaled_center}{the means of each of the covariates used for centering them}
#'
#' @seealso \code{\link{lm_robust}}
#'
#' #' @references
#' Freedman, David A. 2008. "On Regression Adjustments in Experiments with Several Treatments." The Annals of Applied Statistics. JSTOR, 176-96. \url{https://doi.org/10.1214/07-AOAS143}.
#'
#' Lin, Winston. 2013. "Agnostic Notes on Regression Adjustments to Experimental Data: Reexamining Freedman's Critique." The Annals of Applied Statistics 7 (1). Institute of Mathematical Statistics: 295-318. \url{https://doi.org/10.1214/12-AOAS583}.
#'
#' @examples
#' library(fabricatr)
#' library(randomizr)
#' dat <- fabricate(
#'   N = 40,
#'   x = rnorm(N, mean = 2.3),
#'   x2 = rpois(N, lambda = 2),
#'   x3 = runif(N),
#'   y0 = rnorm(N) + x,
#'   y1 = rnorm(N) + x + 0.35
#' )
#'
#' dat$z <- complete_ra(N = nrow(dat))
#' dat$y <- ifelse(dat$z == 1, dat$y1, dat$y0)
#'
#' # Same specification as lm_robust() with one additional argument
#' lmlin_out <- lm_lin(y ~ z, covariates = ~ x, data = dat)
#' tidy(lmlin_out)
#'
#' # Works with multiple pre-treatment covariates
#' lm_lin(y ~ z, covariates = ~ x + x2, data = dat)
#'
#' # Also centers data AFTER evaluating any functions in formula
#' lmlin_out2 <- lm_lin(y ~ z, covariates = ~ x + log(x3), data = dat)
#' lmlin_out2$scaled_center["log(x3)"]
#' mean(log(dat$x3))
#'
#' # Works easily with clusters
#' dat$clusterID <- rep(1:20, each = 2)
#' dat$z_clust <- cluster_ra(clusters = dat$clusterID)
#'
#' lm_lin(y ~ z_clust, covariates = ~ x, data = dat, clusters = clusterID)
#'
#' # Works with multi-valued treatments
#' dat$z_multi <- sample(1:3, size = nrow(dat), replace = TRUE)
#' lm_lin(y ~ z_multi, covariates = ~ x, data = dat)
#'
#' # Stratified estimator with blocks
#' dat$blockID <- rep(1:5, each = 8)
#' dat$z_block <- block_ra(blocks = dat$blockID)
#'
#' lm_lin(y ~ z_block, ~ factor(blockID), data = dat)
#'
#' \dontrun{
#'   # Can also use 'margins' package if you have it installed to get
#'   # marginal effects
#'   library(margins)
#'   lmlout <- lm_lin(y ~ z_block, ~ x, data = dat)
#'   summary(margins(lmlout))
#'
#'   # Can output results using 'texreg'
#'   library(texreg)
#'   texregobj <- extract(lmlout)
#' }
#'
#' @export
lm_lin <- function(formula,
                   covariates,
                   data,
                   weights,
                   subset,
                   clusters,
                   se_type = NULL,
                   ci = TRUE,
                   alpha = .05,
                   return_vcov = TRUE,
                   try_cholesky = FALSE) {

  # ----------
  # Check arguments
  # ----------

  if (length(all.vars(f_rhs(formula))) > 1) {
    stop(
      "The `formula` argument, `", format(formula), "`, must only have the ",
      "treatment variable on the right-hand side of the formula. Covariates ",
      "should be specified in the `covariates` argument like so:\n`covariates = ",
      paste0("~ ", paste(all.vars(f_rhs(formula))[-1], sep = " + ")), "`.",
      "\n\n See ?lm_lin."
    )
  }

  if (class(covariates) != "formula") {
    stop(
      "The `covariates` argument must be specified as a formula:\n",
      "You passed an object of class ", class(covariates)
    )
  }

  cov_terms <- terms(covariates)

  # Check covariates is right hand sided fn
  if (attr(cov_terms, "response") != 0) {
    stop(
      "Must not specify a response variable in `covariates`` formula.\n",
      "`covariates` must be a right-sided formula, such as '~ x1 + x2 + x3'"
    )
  }

  if (length(attr(cov_terms, "order")) == 0) {
    stop(
      "`covariates` must have a variable on the right-hand side, not 0 or 1"
    )
  }

  # ----------
  # Get design matrix including `covariates` for centering
  # ----------

  full_formula <-
    update(
      formula,
      reformulate(
        c(".", labels(cov_terms), response = ".")
      )
    )

  datargs <- enquos(
    formula = full_formula,
    weights = weights,
    subset = subset,
    cluster = clusters
  )
  data <- enquo(data)
  model_data <- clean_model_data(data = data, datargs)

  outcome <- as.matrix(model_data$outcome)
  n <- nrow(outcome)
  design_matrix <- model_data$design_matrix
  weights <- model_data$weights
  cluster <- model_data$cluster

  # Get treatment columns
  has_intercept <- attr(terms(formula), "intercept")
  treat_col <- which(attr(design_matrix, "assign") == 1)
  treatment <- design_matrix[, treat_col, drop = FALSE]
  design_mat_treatment <- colnames(design_matrix)[treat_col]

  # Check case where treatment is not factor and is not binary
  if (any(!(treatment %in% c(0, 1)))) {
    # create dummies for non-factor treatment variable

    # Drop out first group if there is an intercept
    vals <- sort(unique(treatment))
    if (has_intercept) vals <- vals[-1]

    n_treats <- length(vals)
    # Could warn if there are too many values
    # (i.e., if there are as many treatments as observations)

    names(vals) <- paste0(colnames(design_matrix)[treat_col], vals)


    # Create matrix of dummies
    treatment <-
      outer(
        drop(treatment),
        vals,
        function(x, y) as.numeric(x == y)
      )
  }


  # ----------
  # Center and interact variables
  # ----------

  # Initialize as non-demeaned
  demeaned_covars <-
    design_matrix[
      ,
      setdiff(colnames(design_matrix), c(design_mat_treatment, "(Intercept)")),
      drop = FALSE
    ]

  # Choose what to center on!
  if (is.numeric(weights)) {
    center <- apply(demeaned_covars, 2, weighted.mean, weights)
  } else {
    center <- colMeans(demeaned_covars)
  }

  demeaned_covars <- sweep(demeaned_covars, 2, center)

  original_covar_names <- colnames(demeaned_covars)

  # Change name of centered covariates to end in "_c"
  # If covar name has `:` or a `(` not in the first position,
  # wrap the whole var name in parentheses first
  colnames(demeaned_covars) <- paste0(
    ifelse(grepl("\\:|(^.+\\()", colnames(demeaned_covars)),
           paste0("(", colnames(demeaned_covars), ")"),
           colnames(demeaned_covars)),
    "_c"
  )

  n_treat_cols <- ncol(treatment)
  n_covars <- ncol(demeaned_covars)

  # Interacted
  # n_int_covar_cols <- n_covars * (n_treat_cols + has_intercept)
  n_int_covar_cols <- n_covars * (n_treat_cols)
  interacted_covars <- matrix(0, nrow = n, ncol = n_int_covar_cols)
  interacted_covars_names <- character(n_int_covar_cols)
  for (i in 1:n_covars) {
    covar_name <- colnames(demeaned_covars)[i]

    cols <- (i - 1) * n_treat_cols + (1:n_treat_cols)
    interacted_covars[, cols] <- treatment * demeaned_covars[, i]
    interacted_covars_names[cols] <- paste0(colnames(treatment), ":", covar_name)
  }
  colnames(interacted_covars) <- interacted_covars_names

  if (has_intercept) {
    # Have to manually create intercept if treatment wasn't a factor
    X <- cbind(
      matrix(1, nrow = n, ncol = 1, dimnames = list(NULL, "(Intercept)")),
      treatment,
      demeaned_covars,
      interacted_covars
    )
  } else {
    # If no intercept, but treatment is only one column,
    # need to add base terms for covariates
    if (n_treat_cols == 1) {
      X <- cbind(
        treatment,
        demeaned_covars,
        interacted_covars
      )
    } else {
      X <- cbind(
        treatment,
        interacted_covars
      )
    }
  }

  # ----------
  # Estimation
  # ----------

  return_list <-
    lm_robust_fit(
      y = outcome,
      X = X,
      weights = weights,
      cluster = cluster,
      ci = ci,
      se_type = se_type,
      alpha = alpha,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky,
      has_int = has_intercept,
      iv_stage = list(0)
    )

  # ----------
  # Build return object
  # ----------

  return_list <- lm_return(
    return_list,
    model_data = model_data,
    formula = formula
  )

  return_list[["scaled_center"]] <- center
  setNames(return_list[["scaled_center"]], original_covar_names)

  return_list[["call"]] <- match.call()


  return(return_list)
}
