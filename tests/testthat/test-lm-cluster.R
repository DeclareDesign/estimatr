context("Estimator - lm_robust, clustered")

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


rmcall <- function(obj) {
  obj[["call"]] <- NULL
  obj
}


test_that("lm cluster se", {
  N <- 100
  dat <- data.frame(
    Y = rnorm(N),
    Z = rbinom(N, 1, .5),
    X = rnorm(N),
    J = sample(1:10, N, replace = T),
    W = runif(N)
  )


  ## Test functionality
  lm_robust(Y ~ Z, clusters = J, data = dat)

  lm_robust(Y ~ Z + X, clusters = J, data = dat)

  lm_robust(
    Y ~ Z + X,
    clusters = J,
    data = dat
  )


  lm_robust(
    Y ~ Z + X,
    clusters = J,
    se_type = "stata",
    data = dat,
    ci = T
  )

  expect_equivalent(
    as.matrix(
      tidy(
        lm_robust(
          Y ~ X + Z,
          clusters = J,
          ci = FALSE,
          data = dat
        )
      )[, c("p.value", "ci.lower", "ci.upper")]
    ),
    matrix(NA, nrow = 3, ncol = 3)
  )

  ## Test equality
  lm_interact <-
    lm_robust(
      Y ~ Z * X,
      clusters = J,
      data = dat
    )

  lm_interact_stata <-
    lm_robust(
      Y ~ Z * X,
      clusters = J,
      se_type = "stata",
      data = dat
    )

  lm_interact_simple <- lm(Y ~ Z * X, data = dat)

  bm_interact <-
    BMlmSE(
      lm_interact_simple,
      clustervar = as.factor(dat$J),
      IK = FALSE
    )

  bm_interact

  bm_interact_interval <-
    lm_interact_simple$coefficients["Z:X"] +
    qt(0.975, df = bm_interact$dof["Z:X"]) * bm_interact$se["Z:X"] * c(-1, 1)

  bm_interact_stata_interval <-
    lm_interact_simple$coefficients["Z:X"] +
    qt(0.975, df = length(unique(dat$J)) - 1) * bm_interact$se.Stata["Z:X"] * c(-1, 1)

  expect_equivalent(
    tidy(lm_interact)[4, c("std.error", "ci.lower", "ci.upper")],
    c(bm_interact$se["Z:X"], bm_interact_interval)
  )

  expect_equivalent(
    tidy(lm_interact_stata)[4, c("std.error", "ci.lower", "ci.upper")],
    c(bm_interact$se.Stata["Z:X"], bm_interact_stata_interval)
  )


  lm_full <-
    lm_robust(
      Y ~ Z + X,
      clusters = J,
      data = dat
    )

  lm_full_simple <- lm(Y ~ Z + X, data = dat)

  bm_full <-
    BMlmSE(
      lm_full_simple,
      clustervar = as.factor(dat$J),
      IK = FALSE
    )

  bm_full_moe <- qt(0.975, df = bm_full$dof) * bm_full$se
  bm_full_lower <- lm_full_simple$coefficients - bm_full_moe
  bm_full_upper <- lm_full_simple$coefficients + bm_full_moe

  expect_equivalent(
    as.matrix(tidy(lm_full)[, c("std.error", "ci.lower", "ci.upper")]),
    cbind(bm_full$se, bm_full_lower, bm_full_upper)
  )

  ## Works with rank deficient case
  dat$X2 <- dat$X
  lmr_rd <- lm_robust(Y ~ X + Z + X2, data = dat, clusters = J, se_type = "stata")
  lmr_full <- lm_robust(Y ~ X + Z, data = dat, clusters = J, se_type = "stata")
  expect_identical(
    tidy(lmr_rd)[1:3, ],
    tidy(lmr_full)
  )

  lmr_rd_cr2 <- lm_robust(Y ~ X + Z + X2, data = dat, clusters = J, se_type = "CR2")
  lmr_full_cr2 <- lm_robust(Y ~ X + Z, data = dat, clusters = J, se_type = "CR2")
  expect_identical(
    tidy(lmr_rd_cr2)[1:3, ],
    tidy(lmr_full_cr2)
  )

  ## Test error handling
  expect_error(
    lm_robust(
      Y ~ Z,
      clusters = J,
      se_type = "HC2",
      data = dat
    ),
    "CR2"
  )

  expect_error(
    lm_robust(
      Y ~ Z,
      se_type = "CR2",
      data = dat
    ),
    "CR2"
  )

  # To easily do with and without weights
  test_lm_cluster_variance <- function(w) {
    # Test other estimators
    lm_cr0 <- lm_robust(Y ~ Z + X, data = dat, weights = w, clusters = J, se_type = "CR0")
    lm_stata <- lm_robust(Y ~ Z + X, data = dat, weights = w, clusters = J, se_type = "stata")
    lm_cr2 <- lm_robust(Y ~ Z + X, data = dat, weights = w, clusters = J, se_type = "CR2")

    # Stata is the same as CR0 but with finite sample
    expect_equivalent(
      lm_cr0$se ^ 2,
      lm_stata$se ^ 2 * (N - length(lm_stata$coefficients)) * (length(unique(dat$J)) - 1) / ((N - 1) * length(unique(dat$J)))
    )

    expect_false(all(lm_cr0$se == lm_stata$se))
    expect_false(all(lm_cr0$se == lm_cr2$se))
    expect_false(all(lm_stata$se == lm_cr2$se))
    expect_false(all(lm_stata$df == lm_cr2$df))

    expect_equivalent(
      lm_cr0$df,
      lm_stata$df
    )
  }

  # No weights first
  test_lm_cluster_variance(NULL)
  test_lm_cluster_variance(dat$W)

})

test_that("Clustered weighted SEs are correct", {
  skip_if_not_installed("clubSandwich")

  lm_cr2 <- lm_robust(mpg ~ hp, data = mtcars, weights = wt, clusters = cyl, se_type = "CR2")
  lm_stata <- lm_robust(mpg ~ hp, data = mtcars, weights = wt, clusters = cyl, se_type = "stata")
  lm_cr0 <- lm_robust(mpg ~ hp, data = mtcars, weights = wt, clusters = cyl, se_type = "CR0")

  lm_o <- lm(mpg ~ hp, data = mtcars, weights = wt)

  expect_equivalent(
    vcov(lm_cr2),
    as.matrix(clubSandwich::vcovCR(lm_o, cluster = mtcars$cyl, type = "CR2"))
  )
  expect_equivalent(
    vcov(lm_cr0),
    as.matrix(clubSandwich::vcovCR(lm_o, cluster = mtcars$cyl, type = "CR0"))
  )
  expect_equivalent(
    vcov(lm_stata),
    as.matrix(clubSandwich::vcovCR(lm_o, cluster = mtcars$cyl, type = "CR1S"))
  )
})

test_that("lm cluster se with missingness", {
  dat <- data.frame(
    Y = rnorm(100),
    Z = rbinom(100, 1, .5),
    X = rnorm(100),
    J = sample(1:10, 100, replace = T),
    W = runif(100)
  )

  dat$X[23] <- NA
  dat$J[63] <- NA

  expect_warning(
    estimatr_cluster_out <- lm_robust(
      Y ~ Z + X,
      clusters = J,
      data = dat
    ),
    "missingness in the cluster"
  )

  estimatr_cluster_sub <- lm_robust(
    Y ~ Z + X,
    clusters = J,
    data = dat[-c(23, 63), ]
  )

  estimatr_cluster_out[["call"]] <- NULL
  estimatr_cluster_sub[["call"]] <- NULL
  expect_equal(
    estimatr_cluster_out,
    estimatr_cluster_sub
  )
})

test_that("lm works with quoted or unquoted vars and withor without factor clusters", {
  dat <- data.frame(
    Y = rnorm(100),
    Z = rbinom(100, 1, .5),
    X = rnorm(100),
    J = sample(1:10, 100, replace = T),
    W = runif(100)
  )

  lmr <- lm_robust(Y~Z, data = dat, weights = W)
  lmrq <- lm_robust(Y~Z, data = dat, weights = W)
  lmr[["call"]] <- NULL
  lmrq[["call"]] <- NULL
  expect_equal(
    lmr,
    lmrq
  )

  # works with char
  dat$J <- as.character(dat$J)

  lmrc <- lm_robust(Y~Z, data = dat, clusters = J)
  lmrcq <- lm_robust(Y~Z, data = dat, clusters = J)
  lmrc[["call"]] <- NULL
  lmrcq[["call"]] <- NULL
  expect_equal(
    lmrc,
    lmrcq
  )

  # works with num
  dat$J_num <- as.numeric(dat$J)

  lmrc_qnum <- lm_robust(Y~Z, data = dat, clusters = J_num)
  lmrc_qnum[["call"]] <- NULL
  expect_equal(
    lmrc,
    lmrc_qnum
  )


  # works with factor
  dat$J_fac <- as.factor(dat$J)
  expect_equivalent(
    rmcall(lm_robust(Y~Z, data = dat, clusters = J_fac)),
    rmcall(lm_robust(Y~Z, data = dat, clusters = J))
  )

  # works with being cast in the call
  lm_robust(Y~Z, data = dat, clusters = as.factor(J))
})

test_that("Clustered SEs work with clusters of size 1", {
  dat <- data.frame(
    Y = rnorm(100),
    X = rnorm(100),
    J = 1:100
  )

  lm_cr2 <- lm_robust(Y ~ X, data = dat, clusters = J)
  lm_stata <- lm_robust(Y ~ X, data = dat, clusters = J, se_type = "stata")
  lmo <- lm(Y ~ X, data = dat)

  bmo <-
    BMlmSE(
      lmo,
      clustervar = as.factor(dat$J),
      IK = FALSE
    )

  expect_equivalent(
    as.matrix(tidy(lm_cr2)[, c("estimate", "std.error", "df")]),
    cbind(lmo$coefficients, bmo$se, bmo$dof)
  )

  expect_equivalent(
    as.matrix(tidy(lm_stata)[, c("estimate", "std.error")]),
    cbind(lmo$coefficients, bmo$se.Stata)
  )
})

test_that("multiple outcomes", {

  skip_if_not_installed("clubSandwich")
  lmo <- lm(cbind(mpg, hp) ~ wt, data = mtcars)
  lmro <- lm_robust(cbind(mpg, hp) ~ wt, data = mtcars, clusters = cyl)

  expect_equivalent(
    as.matrix(clubSandwich::vcovCR(lmo, cluster = mtcars$cyl, type = "CR2")),
    vcov(lmro)
  )

  expect_equivalent(
    as.matrix(clubSandwich::vcovCR(lmo, cluster = mtcars$cyl, type = "CR0")),
    vcov(lm_robust(cbind(mpg, hp) ~ wt, data = mtcars, clusters = cyl, se_type = "CR0"))
  )

  J <- length(unique(mtcars$cyl))
  n <- nrow(mtcars)
  r <- 2

  # Have to manually do correction because clubSandwich uses n*ny and r*ny in place of n and r in
  # stata correction
  expect_equivalent(
    as.matrix(clubSandwich::vcovCR(lmo, cluster = mtcars$cyl, type = "CR0")) *
      ((J * (n - 1)) / ((J - 1) * (n - r))) ,
    vcov(lm_robust(cbind(mpg, hp) ~ wt, data = mtcars, clusters = cyl, se_type = "stata"))
  )

})
