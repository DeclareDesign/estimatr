## BMlmSE.R implements Bell-McCaffrey standard errors
## This code is taken from Michal Kolesar
## https://github.com/kolesarm/Robust-Small-Sample-Standard-Errors
## Only changed the name of one function, df

#'  Compute the inverse square root of a symmetric matrix
#' @param A matrix
MatSqrtInverse <- function(A) {
  ei <- eigen(A, symmetric = TRUE)

  if (min(ei$values) <= 0) {
    warning("Gram matrix doesn't appear to be positive definite")
  }

  d <- pmax(ei$values, 0)
  d2 <- 1 / sqrt(d)
  d2[d == 0] <- 0
  ## diag(d2) is d2 x d2 identity if d2 is scalar, instead we want 1x1 matrix
  ei$vectors %*% (if (length(d2) == 1) d2 else diag(d2)) %*% t(ei$vectors)
}

#' Compute Bell-McCaffrey Standard Errors
#' @param model Fitted model returned by the \code{lm} function
#' @param clustervar Factor variable that defines clusters. If \code{NULL} (or
#'     not supplied), the command computes heteroscedasticity-robust standard
#'     errors, rather than cluster-robust standard errors.
#' @param ell A vector of the same length as the dimension of covariates,
#'     specifying which linear combination \eqn{\ell'\beta} of coefficients
#'     \eqn{\beta} to compute. If \code{NULL}, compute standard errors for each
#'     regressor coefficient
#' @param IK Logical flag only relevant if cluster-robust standard errors are
#'     being computed. Specifies whether to compute the degrees-of-freedom
#'     adjustment using the Imbens-Kolesár method (if \code{TRUE}), or the
#'     Bell-McCaffrey method (if \code{FALSE})
#' @return Returns a list with the following components \describe{
#'
#' \item{vcov}{Variance-covariance matrix estimator. For the case without
#' clustering, it corresponds to the HC2 estimator (see MacKinnon and White,
#' 1985 and the reference manual for the \code{sandwich} package). For the case
#' with clustering, it corresponds to a generalization of the HC2 estimator,
#' called LZ2 in Imbens and Kolesár.}
#'
#' \item{dof}{Degrees-of-freedom adjustment}
#'
#' \item{se}{Standard error}
#'
#' \item{adj.se}{Adjusted standard errors. For \beta_j, they are defined as
#' \code{adj.se[j]=sqrt(vcov[j,j]se*qt(0.975,df=dof)} so that the Bell-McCaffrey
#' confidence intervals are given as \code{coefficients(fm)[j] +- 1.96* adj.se=}
#'
#' \item{se.Stata}{Square root of the cluster-robust variance estimator used in
#' STATA}
#'
#' }
#' @examples
#' ## No clustering:
#' set.seed(42)
#' x <- sin(1:10)
#' y <- rnorm(10)
#' fm <- lm(y~x)
#' BMlmSE(fm)
#' ## Clustering, defining the first six observations to be in cluster 1, the
#' #next two in cluster 2, and the last three in cluster three.
#' clustervar <- as.factor(c(rep(1, 6), rep(2, 2), rep(3, 2)))
#' BMlmSE(fm, clustervar)
BMlmSE <- function(model, clustervar=NULL, ell=NULL, IK=TRUE) {
  X <- model.matrix(model)
  sum.model <- summary.lm(model)
  n <- sum(sum.model$df[1:2])
  K <- model$rank
  XXinv <- sum.model$cov.unscaled # XX^{-1}
  u <- residuals(model)

  ## Compute DoF given G'*Omega*G without calling eigen as suggested by
  ## Winston Lin
  DoF <- function(GG)
    sum(diag(GG)) ^ 2 / sum(GG * GG)
  ## Previously:
  ## lam <- eigen(GG, only.values=TRUE)$values
  ## sum(lam)^2/sum(lam^2)

  ## no clustering
  if (is.null(clustervar)) {
    Vhat <- sandwich::vcovHC(model, type = "HC2")
    Vhat.Stata <- Vhat * NA

    M <- diag(n) - X %*% XXinv %*% t(X) # annihilator matrix
    ## G'*Omega*G
    GOG <- function(ell) {
      Xtilde <- drop(X %*% XXinv %*% ell / sqrt(diag(M)))
      crossprod(M * Xtilde)
    }
  } else {
    if (!is.factor(clustervar)) stop("'clustervar' must be a factor")

    ## Stata
    S <- length(levels(clustervar)) # number clusters
    uj <- apply(u * X, 2, function(x) tapply(x, clustervar, sum))
    Vhat.Stata <- S / (S - 1) * (n - 1) / (n - K) *
      sandwich::sandwich(model, meat = crossprod(uj) / n)

    ## HC2
    tXs <- function(s) {
      Xs <- X[clustervar == s, , drop = FALSE]
      MatSqrtInverse(diag(NROW(Xs)) - Xs %*% XXinv %*% t(Xs)) %*% Xs
    }
    tX <- lapply(levels(clustervar), tXs) # list of matrices

    tu <- split(u, clustervar)
    tutX <- sapply(seq_along(tu), function(i) crossprod(tu[[i]], tX[[i]]))
    Vhat <- sandwich::sandwich(model, meat = tcrossprod(tutX) / n)

    ## DOF adjustment
    tHs <- function(s) {
      Xs <- X[clustervar == s, , drop = FALSE]
      index <- which(clustervar == s)
      ss <- outer(rep(0, n), index) # n x ns matrix of 0
      ss[cbind(index, 1:length(index))] <- 1
      ss - X %*% XXinv %*% t(Xs)
    }
    tH <- lapply(levels(clustervar), tHs) # list of matrices

    Moulton <- function() {
      ## Moulton estimates
      ns <- tapply(u, clustervar, length)
      ssr <- sum(u ^ 2)
      rho <- max((sum(sapply(seq_along(tu), function(i)
        sum(tu[[i]] %o% tu[[i]]))) - ssr) / (sum(ns ^ 2) - n), 0)
      c(sig.eps = max(ssr / n - rho, 0), rho = rho)
    }

    GOG <- function(ell) {
      G <- sapply(
        seq_along(tX),
        function(i) tH[[i]] %*% tX[[i]] %*% XXinv %*% ell
      )
      GG <- crossprod(G)

      ## IK method
      if (IK == TRUE) {
        Gsums <- apply(
          G, 2,
          function(x) tapply(x, clustervar, sum)
        ) # Z'*G
        GG <- Moulton()[1] * GG + Moulton()[2] * crossprod(Gsums)
      }
      GG
    }
  }

  if (!is.null(ell)) {
    se <- drop(sqrt(crossprod(ell, Vhat) %*% ell))
    dof <- DoF(GOG(ell))
    se.Stata <- drop(sqrt(crossprod(ell, Vhat.Stata) %*% ell))
  } else {
    se <- sqrt(diag(Vhat))
    dof <- sapply(seq(K), function(k) DoF(GOG(diag(K)[, k])))
    se.Stata <- sqrt(diag(Vhat.Stata))
  }
  names(dof) <- names(se)

  list(
    vcov = Vhat, dof = dof, adj.se = se * qt(0.975, df = dof) / qnorm(0.975),
    se = se, se.Stata = se.Stata
  )
}
