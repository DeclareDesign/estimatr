#' Ordinary Least Squares with Robust Standard Errors
#'
#' @description This formula fits a linear model, provides a variety of
#' options for robust standard errors, and conducts coefficient tests
#'
#' @param formula an object of class formula, as in \code{\link{lm}}
#' @param data A \code{data.frame}
#' @param weights the bare (unquoted) names of the weights variable in the
#' supplied data.
#' @param subset An optional bare (unquoted) expression specifying a subset
#' of observations to be used.
#' @param clusters An optional bare (unquoted) name of the variable that
#' corresponds to the clusters in the data.
#' @param se_type The sort of standard error sought. Without clustering:
#' "HC0", "HC1" (or "stata", the equivalent), "HC2" (default), "HC3", or
#' "classical". With clustering: "CR0", "CR2" (default), or "stata". are
#' permissible.
#' @param ci A boolean for whether to compute and return pvalues and confidence
#' intervals, TRUE by default.
#' @param alpha The significance level, 0.05 by default.
#' @param coefficient_name a character or character vector that indicates which
#' coefficients should be reported. If left unspecified, returns all
#' coefficients. Especially for models with clustering where only one
#' coefficient is of interest, specifying a coefficient of interest may
#' result in improvements in speed
#' @param return_vcov a boolean for whether to return the variance-covariance
#' matrix for later usage, TRUE by default.
#' @param try_cholesky a boolean for whether to try using a Cholesky
#' decomposition to solve LS instead of a QR decomposition, FALSE by default.
#' Using a Cholesky decomposition may result in speed gains, but should only
#' be used if users are sure their model is full-rank (i.e. there is no
#' perfect multi-collinearity)
#'
#' @details
#'
#' This function does linear regression and provides a variety of standard
#' errors. It takes formulae and data much in the same was as \code{\link{lm}}
#' does, and all auxiliary variables, such as clusters and weights, can be
#' passed either as quoted names of columns, as bare column names, or
#' as a self-contained vector. Examples of usage can be seen below and in the
#' \href{http://estimatr.declaredesign.org/articles/getting-started.html}{Getting Started vignette}.
#'
#' To see what the exact equations and computations each different
#' type of standard error entails, users can see the
#' \href{http://estimatr.declaredesign.org/articles/technical-notes.html}{technical notes in this vignette}.
#' The default estimators for the non-clustered and clustered case have been
#' chosen for desirable similarities to randomization estimators, performance
#' in small samples, or some combination of both.
#'
#' The function estimates the coefficients and standard errors in C++, using
#' the \code{RcppEigen} package. By default, we estimate the coefficients
#' using Column-Pivoting QR decomposition from the Eigen C++ library, although
#' users could get faster solutions by setting \code{try_cholesky = TRUE} to
#' use a Cholesky decomposition instead. This will likely result in quicker
#' solutions, but the algorithm does not reliably detect when there are linear
#' dependencies in the model and may fail silently if they exist.
#'
#' @return \code{lm_robust} returns an object of class \code{"lm_robust"}.
#'
#' The functions \code{summary} and \code{\link{tidy}} can be used to get
#' the results as a \code{data.frame}. To get useful data out of the return,
#' you can use these data frames, you can use the resulting list directly, or
#' you can use the generic accessor functions \code{coef}, \code{vcov},
#' \code{confint}, and \code{predict}.
#'
#' An object of class \code{"lm_robust"} is a list containing at least the
#' following components:
#' \describe{
#'   \item{est}{the estimated coefficients}
#'   \item{se}{the estimated standard errors}
#'   \item{df}{the estimated degrees of freedom}
#'   \item{p}{the p-values from the t-test using \code{est}, \code{se}, and \code{df}}
#'   \item{ci_lower}{the lower bound of the \code{1 - alpha} percent confidence interval}
#'   \item{ci_upper}{the upper bound of the \code{1 - alpha} percent confidence interval}
#'   \item{coefficient_name}{a character vector of coefficient names}
#'   \item{alpha}{the significance level specified by the user}
#'   \item{res_var}{the residual variance, used for uncertainty when using \code{predict}}
#'   \item{N}{the number of observations used}
#'   \item{k}{the number of columns in the design matrix (includes linearly dependent columns!)}
#'   \item{rank}{the rank of the fitted model}
#'   \item{vcov}{the fitted variance covariance matrix}
#'   \item{weighted}{whether or not weights were applied}
#' }
#' We also return \code{terms} and \code{contrasts}, used by \code{predict}.
#'
#' @examples
#' library(fabricatr)
#' dat <- fabricate(
#'   N = 40,
#'   y = rpois(N, lambda = 4),
#'   x = rnorm(N),
#'   z = rbinom(N, 1, prob = 0.4)
#' )
#'
#' # Default is to use HC2 robust standard errors
#' lmro <- lm_robust(y ~ x + z, data = dat)
#'
#' # Can `tidy()` or `summary()` the data to easily get output as a data.frame
#' tidy(lmro)
#' summary(lmro)
#' # Can also get coefficients many ways
#' all.equal(
#'   lmro$est,
#'   coef(lmro),
#'   tidy(lmro)$est,
#'   check.attributes = FALSE
#' )
#'
#' # Can recover classical standard errors
#' lmclassic <- lm_robust(y ~ x + z, data = dat, se_type = "classical")
#' tidy(lmclassic)
#'
#' # Can easily match Stata's robust standard errors
#' lmstata <- lm_robust(y ~ x + z, data = dat, se_type = "stata")
#' tidy(lmstata)
#'
#' # Easy to specify clusters for cluster-robust inference
#' dat$clusterID <- sample(1:10, size = 40, replace = TRUE)
#'
#' lmclust <- lm_robust(y ~ x + z, data = dat, clusters = clusterID)
#' tidy(lmclust)
#'
#' # Can also match Stata's clustered standard errors
#' lm_robust(
#'   y ~ x + z,
#'   data = dat,
#'   clusters = clusterID,
#'   se_type = "stata"
#' )
#'
#' # Works just as LM does with functions in the formula
#' dat$blockID <- rep(c("A", "B", "C", "D"), each = 10)
#'
#' lm_robust(y ~ x + z + factor(blockID), data = dat)
#'
#' # Weights are also easily specified
#' dat$w <- runif(40)
#'
#' lm_robust(
#'   y ~ x + z,
#'   data = dat,
#'   weights = w,
#'   clusters = clusterID
#' )
#'
#' # Subsetting works just as in `lm()`
#' lm_robust(y ~ x, data = dat, subset = z == 1)
#'
#' # One can also choose to set the significance level for different CIs
#' lm_robust(y ~ x + z, data = dat, alpha = 0.1)
#'
#' # Or only return coefficients of interest
#' lmro_z <- lm_robust(y ~ z + x, data = dat, coefficient_name = c("z", "x"))
#' # Resulting data.frame doesn't have the intercept!
#' tidy(lmro_z)
#'
#' @export
lm_robust <- function(formula,
                      data,
                      weights,
                      subset,
                      clusters,
                      se_type = NULL,
                      ci = TRUE,
                      alpha = .05,
                      coefficient_name = NULL,
                      return_vcov = TRUE,
                      try_cholesky = FALSE) {
  where <- parent.frame()
  model_data <- eval(substitute(
    clean_model_data(
      formula = formula,
      data = data,
      subset = subset,
      cluster = clusters,
      weights = weights,
      where = where
    )
  ))

  return_list <-
    lm_robust_fit(
      y = model_data$outcome,
      X = model_data$design_matrix,
      weights = model_data$weights,
      cluster = model_data$cluster,
      ci = ci,
      se_type = se_type,
      alpha = alpha,
      coefficient_name = coefficient_name,
      return_vcov = return_vcov,
      try_cholesky = try_cholesky
    )

  return_list <- lm_return(
    return_list,
    model_data = model_data,
    formula = formula
  )

  return(return_list)
}
