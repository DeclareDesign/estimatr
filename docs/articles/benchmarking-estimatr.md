---
title: 'Benchmarking estimatr'
author: "Luke Sonnet"
output:
  html_document:
    df_print: paged
link-citations: yes
vignette: |
  %\VignetteIndexEntry{Benchmarking estimatr} 
  %\VignetteEngine{knitr::knitr}
  \usepackage[utf8]{inputenc}
---

We built `estimatr` to provide accurate standard errors **quickly**.
This document benchmarks the speed of or linear regression estimator
against other estimators. Our performance is slightly better than base R
when using classical standard errors, but most of our improvements come
when estimating robust standard errors.

Furthermore, we provide an option in our `lm_robust()` and `lm_lin()`
estimators, `try_cholesky`, which users should set to `TRUE` if they are
concerned about speed and are certain their analysis does not suffer
from perfect multicollinearity (linear dependencies).

Linear regression
=================

I test our speed in estimating coefficients, standard errors, and doing
inference on four different datasets (500 and 5000 observations; 5 and
50 covariates) and across several different specifications. Below I
preview the results comparing `lm_robust()` to a commonly used competing
package, mostly the `sandwich` package and base R. In the two largest
datasets, our method is almost always faster, and is only slower by tiny
margins in the smallest datasets with classical standard errors, where
speed is hardly a factor. When it comes to the biggest gains, using
`lm_robust()` to get HC2 or Stata-like cluster-robust standard errors
will roughly halve your waiting time. If you want CR2 standard errors,
\`lm\_robust() can reduce your run time by a factor of 10!

<table>
<thead>
<tr class="header">
<th>N. Obs</th>
<th>N. Coefs</th>
<th>Estimator</th>
<th>Classical SEs</th>
<th>HC2 SEs</th>
<th>Stata clustered SEs</th>
<th>CR2 SEs</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>500</td>
<td>5</td>
<td><code>estimatr::lm_robust()</code></td>
<td>2.1</td>
<td><strong>2.4</strong></td>
<td><strong>2.2</strong></td>
<td><strong>3.9</strong></td>
</tr>
<tr class="even">
<td></td>
<td></td>
<td>base + sandwich/clubSandwich</td>
<td><strong>1.8</strong></td>
<td>5.1</td>
<td>4.8</td>
<td>44.4</td>
</tr>
<tr class="odd">
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td>5000</td>
<td>5</td>
<td><code>estimatr::lm_robust()</code></td>
<td><strong>5.2</strong></td>
<td><strong>9</strong></td>
<td><strong>7.9</strong></td>
<td><strong>150</strong></td>
</tr>
<tr class="odd">
<td></td>
<td></td>
<td>base + sandwich/clubSandwich</td>
<td><strong>5.2</strong></td>
<td>25</td>
<td>22.8</td>
<td>1882</td>
</tr>
<tr class="even">
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
</tr>
<tr class="odd">
<td>500</td>
<td>50</td>
<td><code>estimatr::lm_robust()</code></td>
<td><strong>6.2</strong></td>
<td><strong>10</strong></td>
<td><strong>8</strong></td>
<td><strong>56</strong></td>
</tr>
<tr class="even">
<td></td>
<td></td>
<td>base + sandwich/clubSandwich</td>
<td>7.5</td>
<td>22</td>
<td>29</td>
<td>162</td>
</tr>
<tr class="odd">
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td>5000</td>
<td>50</td>
<td><code>estimatr::lm_robust()</code></td>
<td><strong>30</strong></td>
<td><strong>44</strong></td>
<td><strong>44</strong></td>
<td><strong>2755</strong></td>
</tr>
<tr class="odd">
<td></td>
<td></td>
<td>base + sandwich/clubSandwich</td>
<td>38</td>
<td>181</td>
<td>216</td>
<td>11893</td>
</tr>
</tbody>
</table>

The times are milliseconds and are a median over 200 runs for all but
the CR2 case, which was taken on a sample of 50 runs, using the
`microbenchmark` package. This benchmarking was done on a 2017 MacBook
Air, with a 1.8 GHz Intel Core i5 CPU and 8 GB of memory.

To see the exact specifications specifications, see below.

    library(estimatr)
    library(microbenchmark)

    # Create some data sets of different sizes for testing below
    set.seed(42)
    data_size <- expand.grid(list(ns = c(500, 5000), ps = c(5, 50)))
    data_list <- lapply(
      1:nrow(data_size), 
      function(i) {
        n <- data_size$ns[i]
        p <- data_size$ps[i]
        y <- rnorm(n)
        X <- matrix(rnorm(n*p), n, p)
        return(data.frame(y, X))
      }
    )

First I compare to a couple other methods of the classical standard
errors. First, let's compare against base R, RcppEigen's `fastLm()`
function (from which we borrow much of our algorithm), and
RcppArmadillo's `fastLm()` function.

    library(RcppEigen)
    #> Warning: package 'RcppEigen' was built under R version 3.4.3
    library(RcppArmadillo)
    #> Warning: package 'RcppArmadillo' was built under R version 3.4.3
    #> 
    #> Attaching package: 'RcppArmadillo'
    #> The following objects are masked from 'package:RcppEigen':
    #> 
    #>     fastLm, fastLmPure

    test_base <- lapply(data_list, function(dat) {
      mbo <- summary(microbenchmark(
        'lm_robust' = lm_robust(y ~ ., data = dat, se_type = "classical"),
        'base' = summary(lm(y ~ ., data = dat)),
        'RcppEigen' = RcppEigen:::summary.fastLm(
          RcppEigen::fastLm(y ~ ., data = dat)
        ),
        "RcppArmadillo" = RcppArmadillo:::summary.fastLm(
          RcppArmadillo::fastLm(y ~ ., data = dat)
        ),
        times = 200L
      ),
      unit = "ms")
      return(mbo[, c("expr", "median")])
    })

The following table has the median time in milliseconds across 50 runs
of each estimator for each of the different data sets.

<table>
<thead>
<tr class="header">
<th align="left">Estimator</th>
<th align="right">N=500, P=5</th>
<th align="right">N=500, P=50</th>
<th align="right">N=5000, P=5</th>
<th align="right">N=500, P=50</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">lm_robust</td>
<td align="right">1.9</td>
<td align="right">5.5</td>
<td align="right">6.6</td>
<td align="right">28</td>
</tr>
<tr class="even">
<td align="left">base</td>
<td align="right">1.7</td>
<td align="right">5.1</td>
<td align="right">7.9</td>
<td align="right">34</td>
</tr>
<tr class="odd">
<td align="left">RcppEigen</td>
<td align="right">1.4</td>
<td align="right">5.8</td>
<td align="right">6.2</td>
<td align="right">35</td>
</tr>
<tr class="even">
<td align="left">RcppArmadillo</td>
<td align="right">1.8</td>
<td align="right">6.9</td>
<td align="right">11.3</td>
<td align="right">65</td>
</tr>
</tbody>
</table>
