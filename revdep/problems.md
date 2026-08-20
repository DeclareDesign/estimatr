# clubSandwich (0.7.0)

* GitHub: <https://github.com/jepusto/clubSandwich>
* Email: <mailto:jepusto@gmail.com>
* GitHub mirror: <https://github.com/cran/clubSandwich>

Run `revdepcheck::revdep_details(, "clubSandwich")` for more info

## Newly broken

*   checking examples ... ERROR
     ```
     ...
     > x_exog <- rnorm(N)
     > x_endog <- 0.5 * z1 + 0.3 * z2 + rnorm(N)
     > y <- 1 + 2 * x_endog + 0.5 * x_exog + rnorm(N)
     > dat <- data.frame(y, x_endog, x_exog, z1, z2, cluster)
     > iv_fit <- iv_robust(y ~ x_endog + x_exog | x_exog + z1 + z2, data = dat, 
     +     clusters = cluster, se_type = "CR2")
     > vcovCR(iv_fit)
                   (Intercept)      x_endog       x_exog
     (Intercept)  0.0028574796 -0.001938895 0.0004459813
     x_endog     -0.0019388949  0.012634349 0.0020154301
     x_exog       0.0004459813  0.002015430 0.0059362835
     > conf_int(iv_fit, vcov = "CR2")
            Coef. Estimate     SE d.f. Lower 95% CI Upper 95% CI
      (Intercept)    1.006 0.0535 18.8        0.894        1.118
          x_endog    2.126 0.1124 16.1        1.888        2.364
           x_exog    0.458 0.0770 14.3        0.293        0.623
     > iv_fe <- iv_robust(y ~ x_endog + x_exog | x_exog + z1 + z2, fixed_effects = ~cluster, 
     +     data = dat, clusters = cluster, se_type = "CR2")
     Error in check_se_type(se_type, clustered, has_fe = (fe_rank > 0L), oneway_fe = !is.null(fe_leverage)) : 
       `se_type = "CR2"` requires hat values from the full [X | FE dummies] design matrix and cannot be used with `fixed_effects`.
     To get CR2 SEs, replace `fixed_effects` with explicit dummies:
       lm_robust(y ~ x + factor(fe_var), se_type = "CR2")
     With `fixed_effects`, available SE types are: "CR0", "stata", or "none".
     Calls: withAutoprint ... eval -> iv_robust -> lm_robust_fit -> check_se_type
     Execution halted
     ```

*   checking tests ...
     ```
       Running ‘testthat.R’
      ERROR
     Running the tests in ‘tests/testthat.R’ failed.
     Last 13 lines of output:
       ── Error ('test_lm_robust.R:590:7'): vcovCR works with se_type inherited from lm_robust(). ──
       Error in `check_se_type(se_type, clustered, has_fe = (fe_rank > 0L), oneway_fe = !is.null(fe_leverage))`: `se_type = "CR2"` requires hat values from the full [X | FE dummies] design matrix and cannot be used with `fixed_effects`.
       To get CR2 SEs, replace `fixed_effects` with explicit dummies:
         lm_robust(y ~ x + factor(fe_var), se_type = "CR2")
       With `fixed_effects`, available SE types are: "CR0", "stata", or "none".
       Backtrace:
           ▆
        1. └─estimatr::lm_robust(...) at test_lm_robust.R:590:7
        2.   └─estimatr::lm_robust_fit(...)
        3.     └─estimatr:::check_se_type(...)
       
       [ FAIL 3 | WARN 8 | SKIP 30 | PASS 2186 ]
       Error:
       ! Test failures.
       Execution halted
     ```

# eventstudyr (1.2.0)

* GitHub: <https://github.com/JMSLab/eventstudyr>
* Email: <mailto:santiago.hermo@monash.edu>
* GitHub mirror: <https://github.com/cran/eventstudyr>

Run `revdepcheck::revdep_details(, "eventstudyr")` for more info

## Newly broken

*   checking tests ...
     ```
       Running ‘testthat.R’
      ERROR
     Running the tests in ‘tests/testthat.R’ failed.
     Last 13 lines of output:
       `expected` is a logical vector (TRUE)
       
       ── Failure ('test-EventStudyOLS.R:278:5'): FE = TRUE,
                  TFE = FALSE,
                  cluster = FALSE works ──
       Expected `all.equal(reg$felevels$V1, as.character(unique(df_test_EventStudyOLS$id)))` to be TRUE.
       Differences:
       `actual` is a character vector ('target is NULL, current is character')
       `expected` is a logical vector (TRUE)
       
       
       [ FAIL 3 | WARN 103 | SKIP 0 | PASS 456 ]
       Error:
       ! Test failures.
       Execution halted
     ```

# RCT (1.2)

* Email: <mailto:isidoro.gu@gmail.com>
* GitHub mirror: <https://github.com/cran/RCT>

Run `revdepcheck::revdep_details(, "RCT")` for more info

## Newly broken

*   checking examples ... ERROR
     ```
     ...
     Error in `purrr::map()`:
     ℹ In index: 1.
     Caused by error in `clean_model_data()`:
     ! `fixed_effects` must be a one-sided formula, such as `~ blockID` or `~ block + year`. You passed an object of class numeric.
     estimatr 1.x accepted a bare grouping vector here; wrap it in `~`.
     Backtrace:
          ▆
       1. ├─RCT::impact_eval(...)
       2. │ └─purrr::map(...)
       3. │   └─purrr:::map_("list", .x, .f, ..., .progress = .progress)
       4. │     ├─purrr:::with_indexed_errors(...)
       5. │     │ └─base::withCallingHandlers(...)
       6. │     ├─purrr:::call_with_cleanup(...)
       7. │     └─RCT (local) .f(.x[[i]], ...)
       8. │       ├─... %>% ...
       9. │       └─estimatr::lm_robust(...)
      10. │         └─estimatr:::clean_model_data(data = data, datargs)
      11. │           └─base::stop(...)
      12. ├─dplyr::select(...)
      13. ├─broom::tidy(.)
      14. └─base::.handleSimpleError(...)
      15.   └─purrr (local) h(simpleError(msg, call))
      16.     └─cli::cli_abort(...)
      17.       └─rlang::abort(...)
     Execution halted
     ```

*   checking tests ...
     ```
       Running ‘testthat.R’
      ERROR
     Running the tests in ‘tests/testthat.R’ failed.
     Last 13 lines of output:
         8. │       ├─... %>% ...
         9. │       └─estimatr::lm_robust(...)
        10. │         └─estimatr:::clean_model_data(data = data, datargs)
        11. │           └─base::stop(...)
        12. ├─dplyr::select(...)
        13. ├─broom::tidy(.)
        14. └─base::.handleSimpleError(...)
        15.   └─purrr (local) h(simpleError(msg, call))
        16.     └─cli::cli_abort(...)
        17.       └─rlang::abort(...)
       
       [ FAIL 1 | WARN 1 | SKIP 0 | PASS 88 ]
       Error:
       ! Test failures.
       Execution halted
     ```

# statuser (0.3.1)

* Email: <mailto:urisohn@gmail.com>
* GitHub mirror: <https://github.com/cran/statuser>

Run `revdepcheck::revdep_details(, "statuser")` for more info

## Newly broken

*   checking examples ... ERROR
     ```
     ...
     > lm2(y ~ x, data = panel, fixed_effects = ~ firm_id)
     Call: lm2(formula = y ~ x, data = panel, fixed_effects = ~firm_id)
     
                 estimate  SE.robust SE.classical t.value p.value std.estimate
     intercept   [dropped]   --         --        --       --        --       
     x           .727**      .085       .082      8.527    <.0001    .191     
     firm_id     df=30       --         --        --       --        --       
                   mean missing red.flag
     intercept   --       --        --  
     x           -.066    0         --  
     firm_id     --       0         --  
     
     N = 150  | missing = 0  | df = 149  | R² = 0.951  | SE type: HC3 
     
     Note: to see explanations for lm2() output, run: `lm2_notes()`
     > 
     > # Two-way fixed effects (firm and year)
     > lm2(y ~ x, data = panel, fixed_effects = ~ firm_id + year)
     Error in check_se_type(se_type, clustered, has_fe = (fe_rank > 0L), oneway_fe = !is.null(fe_leverage)) : 
       `se_type = "HC3"` requires hat values from the full [X | FE dummies] design matrix and cannot be used with `fixed_effects`.
     To get HC3 SEs, replace `fixed_effects` with explicit dummies:
       lm_robust(y ~ x + factor(fe_var), se_type = "HC3")
     With `fixed_effects`, available SE types are: "HC0", "HC1", "stata", "classical", or "none".
     Calls: lm2 ... do.call -> <Anonymous> -> lm_robust_fit -> check_se_type
     Execution halted
     ```

*   checking Rd cross-references ... WARNING
     ```
     Missing link(s) in Rd file 'predict.lm2.Rd':
       ‘[estimatr]{predict.lm_robust}’
     
     See section 'Cross-references' in the 'Writing R Extensions' manual.
     ```

