# clubSandwich (0.7.0)

* GitHub: <https://github.com/jepusto/clubSandwich>
* Email: <mailto:jepusto@gmail.com>
* GitHub mirror: <https://github.com/cran/clubSandwich>

Run `revdepcheck::revdep_details(, "clubSandwich")` for more info

## Newly broken

*   checking examples ... ERROR
     ```
     ...
     > library(estimatr)
     > lm_fit <- lm_robust(weight ~ Time + Diet:Time, data = ChickWeight)
     > vcovCR(lm_fit, cluster = ChickWeight$Chick, type = "CR2")
                 (Intercept)       Time Time:Diet2 Time:Diet3 Time:Diet4
     (Intercept)   3.9317840 -0.7851508 -0.2560762 -0.0750498  0.2695861
     Time         -0.7851508  0.4360726 -0.3314209 -0.3436398 -0.3665016
     Time:Diet2   -0.2560762 -0.3314209  1.3680862  0.3475327  0.3482061
     Time:Diet3   -0.0750498 -0.3436398  0.3475327  1.1195632  0.3483146
     Time:Diet4    0.2695861 -0.3665016  0.3482061  0.3483146  0.5416587
     > lm_fit_clust <- lm_robust(weight ~ Time + Diet:Time, data = ChickWeight, 
     +     clusters = Chick)
     > conf_int(lm_fit_clust, vcov = "CR2")
            Coef. Estimate    SE d.f. Lower 95% CI Upper 95% CI
      (Intercept)    27.86 1.983 48.7       23.874        31.84
             Time     7.05 0.660 31.4        5.703         8.40
       Time:Diet2     1.61 1.170 19.0       -0.837         4.06
       Time:Diet3     3.74 1.058 19.0        1.524         5.95
       Time:Diet4     2.86 0.736 18.4        1.317         4.41
     > lin_fit_clust <- lm_lin(weight ~ Diet, covariates = ~Time, data = ChickWeight, 
     +     clusters = Chick)
     > conf_int(lin_fit_clust, vcov = "CR2")
     Error in vcovCR.lm_robust(obj, type = vcov, ...) : 
       You must specify a clustering variable or `obj` must include a clustering variable.
     Calls: withAutoprint ... eval -> eval -> conf_int -> vcovCR -> vcovCR.lm_robust
     Execution halted
     ```

*   checking tests ...
     ```
       Running ‘testthat.R’
      ERROR
     Running the tests in ‘tests/testthat.R’ failed.
     Last 13 lines of output:
         1. ├─testthat::expect_equivalent(...) at test_lm_lin.R:444:3
         2. │ └─testthat::quasi_label(enquo(object), label)
         3. │   └─rlang::eval_bare(expr, quo_get_env(quo))
         4. ├─clubSandwich::vcovCR(...)
         5. └─clubSandwich:::vcovCR.lm_robust(...)
         6.   └─clubSandwich:::vcov_CR(...)
         7.     ├─clubSandwich:::model_matrix(obj)
         8.     └─clubSandwich:::model_matrix.lm_robust(obj)
         9.       ├─stats::model.matrix(frm, data = mf)
        10.       └─stats::model.matrix.default(frm, data = mf)
       
       [ FAIL 20 | WARN 1 | SKIP 30 | PASS 2191 ]
       Error:
       ! Test failures.
       Execution halted
     ```

