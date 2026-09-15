estimatr: Fast Estimators for Design-Based Inference
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![CRAN status](https://www.r-pkg.org/badges/version/estimatr)](https://cran.r-project.org/package=estimatr) [![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/estimatr?color=green)](https://r-pkg.org/pkg/estimatr) [![Build status](https://github.com/DeclareDesign/estimatr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DeclareDesign/estimatr/actions/workflows/R-CMD-check.yaml) **estimatr** is an `R` package providing a range of commonly used linear estimators, designed for speed and for ease of use. Users can easily recover robust, cluster-robust, and other design-appropriate estimates. We include two functions that implement means estimators, `difference_in_means()` and `horvitz_thompson()`, and four linear regression estimators, `lm_robust()`, `lm_lin()`, `iv_robust()`, and `lh_robust()`. In each case, users can choose an estimator to reflect cluster-randomized, block-randomized, and block-and-cluster-randomized designs. The [Getting Started Guide](https://declaredesign.org/r/estimatr/articles/getting-started.html) describes each estimator provided by **estimatr** and how it can be used in your analysis.

You can also see the ways you can [get regression tables out of estimatr](https://declaredesign.org/r/estimatr/articles/regression-tables.html) using `texreg` and `modelsummary`. Fast estimators also enable fast simulation of research designs to learn about their properties (see [DeclareDesign](https://declaredesign.org)).

## Installing estimatr

To install the latest stable release of **estimatr**, please ensure that you are running version 3.6 or later of R and run the following code:

``` r
install.packages("estimatr")
```

## Easy to use

Once the package is installed, getting appropriate estimates and standard errors is both fast and easy.

``` r
library(estimatr)

# sample data from a cluster-randomized experiment
set.seed(343)
clusters <- sample(letters[1:10], size = 100, replace = TRUE)
treated_clusters <- sample(letters[1:10], 5)
dat <- data.frame(
  y = rnorm(100),
  clusterID = clusters,
  z = as.numeric(clusters %in% treated_clusters)
)

# robust standard errors
res_rob <- lm_robust(y ~ z, data = dat)
# tidy data frames on command
tidy(res_rob)
#> # A tibble: 2 × 9
#>   term        estimate std.error statistic  p.value conf.low conf.high    df outcome
#>   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl> <dbl> <chr>  
#> 1 (Intercept)   -0.440     0.153     -2.88 0.00490    -0.743    -0.137    98 y      
#> 2 z              0.760     0.194      3.92 0.000166    0.375     1.15     98 y

# cluster-robust standard errors
res_cl <- lm_robust(y ~ z, data = dat, clusters = clusterID)
# the standard summary view is also available
summary(res_cl)
#> 
#> Call:
#> lm_robust(formula = y ~ z, data = dat, clusters = clusterID)
#> 
#> Standard error type:  CR2 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|) CI Lower CI Upper    DF
#> (Intercept)  -0.4396     0.1256  -3.499 0.025593  -0.7908  -0.0885 3.935
#> z             0.7603     0.2129   3.571 0.007739   0.2661   1.2545 7.707
#> 
#> Multiple R-squared:  0.1388 ,    Adjusted R-squared:   0.13 
#> F-statistic: 12.75 on 1 and 9 DF,  p-value: 0.006011

# matched-pair design learned from the blocks argument
res_dim <- difference_in_means(extra ~ group, data = sleep, blocks = ID)
res_dim
#> Design:  Matched-pair 
#>        Estimate Std. Error  t value   Pr(>|t|)  CI Lower CI Upper DF
#> group2     1.58  0.3889587 4.062128 0.00283289 0.7001142 2.459886  9
```

The [Getting Started Guide](https://declaredesign.org/r/estimatr/articles/getting-started.html) has more examples and uses, as do the reference pages. The [Mathematical Notes](https://declaredesign.org/r/estimatr/articles/mathematical-notes.html) provide more information about what each estimator is doing under the hood.

## Fast to use

Getting estimates and robust standard errors is also faster than it used to be. Compare our package to using `lm()` and the `sandwich` package to get HC2 standard errors. More speed comparisons, including absorbed fixed effects and a head-to-head against estimatr 1.0.6, are in the [Performance](https://declaredesign.org/r/estimatr/articles/performance.html) article.

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
```

| estimatr      | median run-time (ms) |
|:--------------|---------------------:|
| estimatr      |                    3 |
| lm + sandwich |                   13 |

Measured on an Apple M4 under R 4.6.0 with estimatr 2.0.0, 200 replications. The chunk above is not run when this README is knit, so the table is typed in rather than generated.

## estimatr 2.0: a ground-up rewrite

estimatr 2.0 rewrites the package from the ground up for the DeclareDesign workflow: fit the same model thousands of times, as fast as possible, and get the same number every time. If you are coming from 1.x, read `vignette("estimatr2.0")` first. It covers what does not change, what changes and why, and how to port a 1.x script.

Three vignettes ship with the package, and three more (performance, the tidyverse, and regression tables) live on the [website](https://declaredesign.org/r/estimatr/) alone, because their numbers and their examples drift faster than a release cycle:

| vignette | what it is for |
|----|----|
| `estimatr2.0` | porting from 1.x: what changed, why, and the benchmarks |
| `getting-started` | the six estimators, on one worked example |
| `mathematical-notes` | every estimator defined, cited, and then validated against its own definition |

### What stays the same

The six estimators are `lm_robust`, `lm_lin`, `iv_robust`, `lh_robust`, `difference_in_means`, and `horvitz_thompson`, now with a Pashley and Miratrix (2021) blocked-variance estimator that 1.x does not have.

The six estimators keep their 1.0.6 signatures, except `horvitz_thompson()`, whose five probability arguments consolidate into `condition_prs`; the removals are below. Run side by side on one machine, 2.0 and 1.0.6 return the same estimates, standard errors, and degrees of freedom to 1e-12, across every supported standard error type, weighted and unweighted, clustered and unclustered, single and multivariate outcomes. The test suite holds that claim at 695 assertions against answers recorded from an installed 1.0.6, with `data-raw/make_estimatr_reference.R` producing the recording. Those run at 1e-9 rather than 1e-12, because a fixture recorded on one platform meets a different BLAS on another, and the floor there is the linear algebra rather than this package.

### What breaks

Two removals and one default, all deliberate, all covered in the vignette's porting section.

`horvitz_thompson()` takes one probability argument, `condition_prs`, in place of five. `blocks`, `clusters`, `simple`, `ra_declaration`, `condition_pr_mat`, `subset`, and `return_condition_pr_mat` are gone, along with `se_type = "constant"` and the three exported matrix builders that served them (`declaration_to_condition_pr_mat`, `gen_pr_matrix_cluster`, and `permutations_to_condition_pr_mat`). No design is lost: blocked, clustered, and custom designs reach the estimator through an `ra_declaration` passed as `condition_prs`, and that path uses exact design-aware joint probabilities rather than the conservative bound.

`commarobust()` and `starprep()` are removed. Both still exist as names that error and name their replacement, rather than failing with "could not find function".

**Fixed effects without clusters support HC2 and HC3 exactly, at any number of factors, and default to HC2, which is what 1.x returns for the same call**, so that code ports with no change in the numbers at all. The projection decomposes exactly, so no dummy matrix is built to get there. CR2 is the exception: it is built from cluster-level blocks of the hat matrix rather than from the leverage diagonal, so it still expands the dummies, and clustered `fixed_effects` therefore defaults to CR0 where 1.x defaulted to CR2. That is the only default in the release that moves, and it warns once per session rather than moving silently. Writing `se_type = "CR2"` gets the 1.x number back exactly.

A clustered block holding a single treated or single control cluster is refused rather than given a variance. That refusal is a correctness fix: 1.x returns a number there that is too small by roughly the block's cluster count.

`tidy()`, `glance()`, and `augment()` return tibbles, as broom's methods do, where 1.x returned plain data frames. The values are unchanged, and `$`, `[[`, and row indexing work as before. `tidy(fit)[, "estimate"]` now returns a one-column tibble rather than a vector, so write `tidy(fit)$estimate`; code that sets row names on the result should coerce with `as.data.frame()` first.

### Status

`R CMD check --as-cran`: 0 errors, 0 warnings, 1 NOTE (the maintainer change). Test suite 5,906 assertions under `R CMD check`, 0 failures, one skip (`blkvar`, worth 12 assertions, which CRAN does not serve).

Every open estimatr issue was read against this implementation: 26 are fixed here, 23 are feature requests, 7 are out of scope, 6 are not reproducible, 5 are superseded by the rewrite, and 4 remain open. `vignette("estimatr2.0")` names the four.

### How this was written

estimatr 2.0.0 was written with AI. `vignette("mathematical-notes")` says so in full and, in the same document, states each estimator's definition in mathematics and then, immediately underneath, measures estimatr against that definition to machine precision, computed when the vignette is built. Sixteen such checks, and the vignette refuses to build if any of them fails. That is the evidence that should decide whether you install it.

------------------------------------------------------------------------

This project is generously supported by a grant from the [Laura and John Arnold Foundation](https://www.arnoldventures.org) and seed funding from [Evidence in Governance and Politics (EGAP)](https://egap.org/).
