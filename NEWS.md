# estimatr 0.12.0

* Fixed ambiguity about how interacted covariates were centered in `lm_lin`
* A series of fixes for bugs that occurred with multiple outcomes (multivariate regression):
  * Fixed bug pointed out by James Pustejovsky via the `sandwich` version 2.5-0 and off-diagonal blocks of multivariate regression vcov matrices
  * Fixed bugs in `lm_lin` preventing multivariate regression
  * Fixed bug that truncated degrees of freedom with "CR2"  standard errors
  * Fixed bug that returned incorrect R-squared for the second or later outcomes
* Fixed bug preventing integration with latest version of `margins`
* Fixed bug with `difference_in_means` when using `condition1` and `condition2` to subset a treatment vector with more than two treatment conditions. Previous estimates and standard errors were incorrect.

# estimatr 0.10.0

* Changed names of confidence interval columns in tidied data from `ci.lower` and `ci.upper` to `conf.low` and `conf.high` to be in line with other tidy methods
* Added support for `fixed_effects` that are just one block
* Added support for specifing `condition_prs` in `horvitz_thompson()` as a single number
* Added t- and z-statistics to output
* Limit unnecessary messaging in `horvitz_thompson()`

# estimatr 0.8.0

* Added support for absorbing fixed effects in `lm_robust` and `iv_robust`
* Added `commarobust` and `starprep` for stargazer integration
* Added `texreg` support for 2SLS IV models
* Fixed bugs for incorrect F-statistics with robust standard errors
* Refactor of main fitting engine for linear models

# estimatr 0.6.0

* Added support for multivariate linear models
* Added support for instrumental variables regression
* Major change to name of object output elements to mostly match with `broom::tidy`
  * old -> new
  * `coefficient_names` -> `term`
  * `se` -> `std.error`
  * `p` -> `p.values`
  * `ci_lower` -> `ci.lower`
  * `ci_upper` -> `ci.upper`
  * All of the above changes are also made to the column names on the output of `tidy`; furthermore for `tidy` objects one further name change from `coefficients` -> `estimate` has been made
* Fixed bug that caused variances, standard errors, and p-values to be wrong for weighted "CR2" variance estimation
* Fixed incorrect estimates when both weights and blocks were passed to `difference_in_means`
* Rewrite NSE handling to be done by `rlang`
* Rewrite `na.omit` handler in R
* Major refactor of C++ underlying regression estimators

# estimatr 0.4.0

* Changed suffix added to centered variables in `lm_lin()` from `_bar` to `_c`
* Added all vignettes to `.Rbuildignore`, only available on website now
* Fixed `lm_robust_helper.cpp` algorithm to not catch own exception and to deal with `valgrind` memory errors
* Bugfix where passing a formula as an object within a function would fail
* Simplified some tests for various CRAN test platforms

# estimatr 0.2.0

* First **CRAN** upload
