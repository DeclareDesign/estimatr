# estimatr 0.5.0 (GitHub version)

* Fixed bug that caused variances, standard errors, and p-values to be wrong for weighted "CR2" variance estimation
* Added support for multivariate linear models
* Added support for instrumental variables regression
* Rewrite NSE handling to be done by `rlang`

# estimatr 0.4.0 (CRAN version)

* Changed suffix added to centered variables in `lm_lin()` from `_bar` to `_c`
* Added all vignettes to `.Rbuildignore`, only available on website now
* Fixed `lm_robust_helper.cpp` algorithm to not catch own exception and to deal with `valgrind` memory errors
* Bugfix where passing a formula as an object within a function would fail
* Simplified some tests for various CRAN test platforms

# estimatr 0.2.0

* First **CRAN** upload
