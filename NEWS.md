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
