# DesignLibrary (0.1.10)

* GitHub: <https://github.com/DeclareDesign/DesignLibrary>
* Email: <mailto:jjc2247@columbia.edu>
* GitHub mirror: <https://github.com/cran/DesignLibrary>

Run `revdepcheck::revdep_details(, "DesignLibrary")` for more info

## Newly broken

*   checking tests ...
     ```
       Running ‘testthat.R’
      ERROR
     Running the tests in ‘tests/testthat.R’ failed.
     Last 13 lines of output:
         8.         └─DeclareDesign:::future_lapply(...)
         9.           └─base::lapply(...)
        10.             └─DeclareDesign (local) FUN(X[[i]], ...)
        11.               ├─DeclareDesign:::run_design_internal(design)
        12.               └─DeclareDesign:::run_design_internal.design(design)
        13.                 └─DeclareDesign:::next_step(step, current_df, i)
        14.                   └─base::tryCatch(...)
        15.                     └─base (local) tryCatchList(expr, classes, parentenv, handlers)
        16.                       └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
        17.                         └─value[[3L]](cond)
       
       [ FAIL 3 | WARN 18 | SKIP 0 | PASS 323 ]
       Error:
       ! Test failures.
       Execution halted
     ```

## In both

*   checking whether package ‘DesignLibrary’ can be installed ... WARNING
     ```
     Found the following significant warnings:
       Warning: package ‘randomizr’ was built under R version 4.6.1
     See ‘/Users/alexandercoppock/git_projects/estimatr/revdep/checks.noindex/DesignLibrary/new/DesignLibrary.Rcheck/00install.out’ for details.
     ```

