estimatr: Fast Estimators for Design-Based Inference
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![CRAN
status](https://www.r-pkg.org/badges/version/estimatr)](https://cran.r-project.org/package=estimatr)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/estimatr?color=green)](https://r-pkg.org/pkg/estimatr)
[![Build
status](https://github.com/DeclareDesign/estimatr/workflows/R-CMD-check/badge.svg)](https://github.com/DeclareDesign/estimatr/actions)
[![Code
coverage](https://codecov.io/gh/DeclareDesign/estimatr/branch/master/graph/badge.svg?token=x9MpkuKobc)](https://codecov.io/gh/DeclareDesign/estimatr)
[![Replications](https://softwarecite.com/badge/estimatr)](https://softwarecite.com/package/estimatr)

**estimatr** is an `R` package providing a range of commonly-used linear
estimators, designed for speed and for ease-of-use. Users can easily
recover robust, cluster-robust, and other design appropriate estimates.
We include two functions that implement means estimators,
`difference_in_means()` and `horvitz_thompson()`, and three linear
regression estimators, `lm_robust()`, `lm_lin()`, and `iv_robust()`. In
each case, users can choose an estimator to reflect cluster-randomized,
block-randomized, and block-and-cluster-randomized designs. The [Getting
Started
Guide](https://declaredesign.org/r/estimatr/articles/getting-started.html)
describes each estimator provided by **estimatr** and how it can be used
in your analysis.

You can also see the multiple ways you can [get regression tables out of
estimatr](https://declaredesign.org/r/estimatr/articles/regression-tables.html)
using commonly used `R` packages such as `texreg` and `stargazer`. Fast
estimators also enable fast simulation of research designs to learn
about their properties (see [DeclareDesign](https://declaredesign.org)).

## Installing estimatr

To install the latest stable release of **estimatr**, please ensure that
you are running version 3.5 or later of R and run the following code:

``` r
install.packages("estimatr")
```

## Easy to use

Once the package is installed, getting appropriate estimates and
standard errors is now both fast and easy.

``` r
library(estimatr)

# sample data from cluster-randomized experiment
library(fabricatr)
library(randomizr)
dat <- fabricate(
  N = 100,
  y = rnorm(N),
  clusterID = sample(letters[1:10], size = N, replace = TRUE),
  z = cluster_ra(clusterID)
)

# robust standard errors
res_rob <- lm_robust(y ~ z, data = dat)
# tidy dataframes on command!
tidy(res_rob)
#>          term estimate std.error statistic p.value conf.low conf.high df
#> 1 (Intercept)    0.065      0.14      0.46    0.64    -0.21      0.34 98
#> 2           z   -0.067      0.21     -0.32    0.75    -0.48      0.35 98
#>   outcome
#> 1       y
#> 2       y

# cluster robust standard errors
res_cl <- lm_robust(y ~ z, data = dat, clusters = clusterID)
# standard summary view also available
summary(res_cl)
#> 
#> Call:
#> lm_robust(formula = y ~ z, data = dat, clusters = clusterID)
#> 
#> Standard error type:  CR2 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|) CI Lower CI Upper   DF
#> (Intercept)   0.0653      0.145   0.452    0.678   -0.358    0.489 3.53
#> z            -0.0670      0.202  -0.331    0.750   -0.544    0.410 7.05
#> 
#> Multiple R-squared:  0.00105 ,   Adjusted R-squared:  -0.00915 
#> F-statistic: 0.11 on 1 and 9 DF,  p-value: 0.748

# matched-pair design learned from blocks argument
data(sleep)
res_dim <- difference_in_means(extra ~ group, data = sleep, blocks = ID)
```

The [Getting Started Guide](/r/estimatr/articles/getting-started.html)
has more examples and uses, as do the reference pages. The [Mathematical
Notes](/r/estimatr/articles/mathematical-notes.html) provide more
information about what each estimator is doing under the hood.

## Fast to use

Getting estimates and robust standard errors is also faster than it used
to be. Compare our package to using `lm()` and the `sandwich` package to
get HC2 standard errors. More speed comparisons are available
[here](https://declaredesign.org/r/estimatr/articles/benchmarking-estimatr.html).
Furthermore, with many blocks (or fixed effects), users can use the
`fixed_effects` argument of `lm_robust` with HC1 standard errors to
greatly improve estimation speed. More on [fixed effects
here](https://declaredesign.org/r/estimatr/articles/absorbing-fixed-effects.html).

``` r
dat <- data.frame(X = matrix(rnorm(2000*50), 2000), y = rnorm(2000))

library(microbenchmark)
library(lmtest)
library(sandwich)
mb <- microbenchmark(
  `estimatr` = lm_robust(y ~ ., data = dat),
  `lm + sandwich` = {
    lo <- lm(y ~ ., data = dat)
    coeftest(lo, vcov = vcovHC(lo, type = 'HC2'))
  }
)
#> Warning in microbenchmark(estimatr = lm_robust(y ~ ., data = dat), `lm +
#> sandwich` = {: less accurate nanosecond times to avoid potential integer
#> overflows
```

| estimatr      | median run-time (ms) |
|:--------------|---------------------:|
| estimatr      |                  6.3 |
| lm + sandwich |                 20.1 |

------------------------------------------------------------------------

This project is generously supported by a grant from the [Laura and John
Arnold Foundation](http://www.arnoldfoundation.org) and seed funding
from [Evidence in Governance and Politics (EGAP)](http://egap.org).
