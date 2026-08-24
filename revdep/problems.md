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

# projoint (1.1.3)

* GitHub: <https://github.com/yhoriuchi/projoint>
* Email: <mailto:yusaku.horiuchi@gmail.com>
* GitHub mirror: <https://github.com/cran/projoint>

Run `revdepcheck::revdep_details(, "projoint")` for more info

## Newly broken

*   checking running R code from vignettes ...
     ```
       ‘analyze.Rmd’ using ‘UTF-8’... failed
       ‘correct.Rmd’ using ‘UTF-8’... OK
       ‘design.Rmd’ using ‘UTF-8’... OK
       ‘explore.Rmd’ using ‘UTF-8’... OK
       ‘faq.Rmd’ using ‘UTF-8’... OK
       ‘read.Rmd’ using ‘UTF-8’... OK
       ‘structure.Rmd’ using ‘UTF-8’... OK
      ERROR
     Errors in running code in vignettes:
     when running code in ‘analyze.Rmd’
       ...
     
     > predicted_irr <- predict_tau(out1_arranged)
     
       When sourcing ‘analyze.R’:
     Error: ℹ In argument: `estimatr::tidy(...)`.
     ℹ In group 1: `x = 1`.
     Caused by error in `prep_data()`:
     ! `clusters` has only one level, so there is no between-cluster variation for a cluster-robust variance to estimate.
     Drop `clusters`, or use `se_type = "classical"` if the clustering is not the level you meant to inference over.
     Execution halted
     ```

