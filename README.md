
<!-- README.md is generated from README.Rmd. Please edit that file -->
estimatr: Fast Estimators for Design-Based Inference
====================================================

[![Travis-CI Build Status](https://travis-ci.org/DeclareDesign/estimatr.svg?branch=master)](https://travis-ci.org/DeclareDesign/estimatr) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/DeclareDesign/estimatr?branch=master&svg=true)](https://ci.appveyor.com/project/DeclareDesign/estimatr) [![Coverage Status](https://coveralls.io/repos/github/DeclareDesign/estimatr/badge.svg?branch=master)](https://coveralls.io/github/DeclareDesign/estimatr?branch=master)

Technical papers and textbooks demand complex estimation strategies that are often difficult to implement, even for scientists who are expert coders. The result is slow code copied and pasted from the internet, where the result is taken on faith.

**estimatr** provides a small set of commonly-used estimators (methods for estimating quantities of interest like treatment effects or regression parameters), using `C++` for speed, and implemented in `R` with simple, accessible syntax. We include two functions that implement means estimators, [`difference_in_means()`](http://estimatr.declaredesign.org/reference/difference_in_means.html) and [`horvitz_thompson()`](http://estimatr.declaredesign.org/reference/horvitz_thompson.html). In addition, we include two functions for linear regression estimators, [`lm_robust()`](http://estimatr.declaredesign.org/reference/lm_robust.html) and [`lm_lin()`](http://estimatr.declaredesign.org/reference/lm_lin.html). In each case, scientists can choose an estimator to reflect cluster-randomized, block-randomized, and block-and-cluster-randomized designs. The [Getting Started Guide](http://estimatr.declaredesign.org/articles/getting-started.html) describes each estimator provided by **estimatr** and how it can be used in your analysis.

Fast estimators also enable fast simulation of research designs to learn about their properties (see [DeclareDesign](http://declaredesign.org)).

### Installing estimatr

If you would like to use the current development release of **estimatr**, please ensure that you are running version 3.3 or later of R and run the following code:

``` r
install.packages("estimatr", dependencies = TRUE,
  repos = c("http://R.declaredesign.org", "https://cloud.r-project.org"))
```

### Easy to use

Once the package is installed, getting appropriate estimates and standard errors is now both fast and easy.

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

# cluster robust standard errors
res_cl <- lm_robust(y ~ z, data = dat, clusters = clusterID)
summary(res_cl)
#> 
#> Call:
#> lm_robust(formula = y ~ z, data = dat, clusters = clusterID)
#> 
#> Standard error type =  CR2 
#> 
#> Coefficients:
#>             Estimate Std. Error Pr(>|t|) CI Lower CI Upper   DF
#> (Intercept)   -0.180      0.200    0.428   -0.777    0.417 3.38
#> z              0.424      0.236    0.118   -0.142    0.990 6.54
#> 
#> Multiple R-squared:  0.0419 ,    Adjusted R-squared:  0.0321 
#> F-statistic: 4.28 on 1 and 98 DF,  p-value: 0.0411

# matched-pair design learned from blocks argument
data(sleep)
res_dim <- difference_in_means(extra ~ group, data = sleep, blocks = ID)
```

The [Getting Started Guide](http://estimatr.declaredesign.org/articles/getting-started.html) has more examples and uses, as do the reference pages. The [Technical Notes](http://estimatr.declaredesign.org/articles/technical-notes.html) provide more information about what each estimator is doing under the hood

### Fast to use

Getting estimates and robust standard errors is also faster than it used to be. Compare our package to using `lm()` and the `sandwich` package to get HC2 standard errors. More speed comparisons are available [here](http://estimatr.declaredesign.org/articles/benchmarking%20estimatr.html).

``` r
# example code
# estimatr
lm_robust(y ~ x1 + x2 + x3 + x4, data = dat)

# usual specification (lm + sandwich)
library(lmtest)
library(sandwich)
lm_out <- lm(y ~ x1 + x2 + x3 + x4, data = dat)
coeftest(lm_out, vcov = vcovHC(lm_out, type = 'HC2'))
```

![](vignettes/lm_speed.png)

![](vignettes/lm_speed_covars.png)

------------------------------------------------------------------------

This project is generously supported by a grant from the [Laura and John Arnold Foundation](http://www.arnoldfoundation.org) and seed funding from [Evidence in Governance and Politics (EGAP)](http://egap.org).
