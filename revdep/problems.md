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

