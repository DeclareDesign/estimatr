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
       Differences:
       `actual` is a character vector ('target is NULL, current is character')
       `expected` is a logical vector (TRUE)
       
       ── Failure ('test-EventStudyPlot.R:230:5'): computed smoothest path for examples is within expectations ──
       Expected `p$data$smoothest_path` to equal `matrix(rep(0, nrow(p$data)))`.
       Differences:
       `dim(actual)` is absent
       `dim(expected)` is an integer vector (16, 1)
       
       
       [ FAIL 4 | WARN 112 | SKIP 0 | PASS 455 ]
       Error:
       ! Test failures.
       Execution halted
     ```

# hbal (1.2.15)

* GitHub: <https://github.com/xuyiqing/hbal>
* Email: <mailto:yiqingxu@stanford.edu>
* GitHub mirror: <https://github.com/cran/hbal>

Run `revdepcheck::revdep_details(, "hbal")` for more info

## Newly broken

*   checking examples ... WARNING
     ```
     Found the following significant warnings:
     
       Warning: Setting row names on a tibble is deprecated.
       Warning: Setting row names on a tibble is deprecated.
       Warning: Setting row names on a tibble is deprecated.
     Deprecated functions may be defunct as soon as of the next release of
     R.
     See ?Deprecated.
     ```

