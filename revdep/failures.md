# hbal (1.2.15)

* GitHub: <https://github.com/xuyiqing/hbal>
* Email: <mailto:yiqingxu@stanford.edu>
* GitHub mirror: <https://github.com/cran/hbal>

Run `revdepcheck::revdep_details(, "hbal")` for more info

## In both

*   checking whether package ‘hbal’ can be installed ... ERROR
     ```
     Installation failed.
     See ‘/Users/alexandercoppock/git_projects/estimatr/revdep/checks.noindex/hbal/new/hbal.Rcheck/00install.out’ for details.
     ```

## Installation

### Devel

```
* installing *source* package ‘hbal’ ...
** this is package ‘hbal’ version ‘1.2.15’
** package ‘hbal’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
using C++ compiler: ‘Apple clang version 21.0.0 (clang-2100.1.1.101)’
using SDK: ‘MacOSX26.5.sdk’
clang++ -arch arm64 -std=gnu++20 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/alexandercoppock/git_projects/estimatr/revdep/library.noindex/estimatr/new/Rcpp/include' -I'/Users/alexandercoppock/git_projects/estimatr/revdep/library.noindex/estimatr/new/RcppEigen/include' -I/opt/R/arm64/include     -fPIC  -falign-functions=64 -Wall -g -O2   -c RcppExports.cpp -o RcppExports.o
In file included from RcppExports.cpp:4:
In file included from /Users/alexandercoppock/git_projects/estimatr/revdep/library.noindex/estimatr/new/RcppEigen/include/RcppEigen.h:25:
...
      |       ^
5 warnings generated.
clang++ -arch arm64 -std=gnu++20 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -L/Library/Frameworks/R.framework/Resources/lib -L/opt/R/arm64/lib -o hbal.so RcppExports.o eigen.o -L/Library/Frameworks/R.framework/Resources/lib -lRlapack -L/Library/Frameworks/R.framework/Resources/lib -lRblas -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0 -L/opt/gfortran/lib -lemutls_w -lheapt_w -lgfortran -lquadmath -F/Library/Frameworks/R.framework/.. -framework R
ld: warning: search path '/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0' not found
ld: warning: search path '/opt/gfortran/lib' not found
ld: library 'emutls_w' not found
clang++: error: linker command failed with exit code 1 (use -v to see invocation)
make: *** [hbal.so] Error 1
ERROR: compilation failed for package ‘hbal’
* removing ‘/Users/alexandercoppock/git_projects/estimatr/revdep/checks.noindex/hbal/new/hbal.Rcheck/hbal’


```
### CRAN

```
* installing *source* package ‘hbal’ ...
** this is package ‘hbal’ version ‘1.2.15’
** package ‘hbal’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
using C++ compiler: ‘Apple clang version 21.0.0 (clang-2100.1.1.101)’
using SDK: ‘MacOSX26.5.sdk’
clang++ -arch arm64 -std=gnu++20 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/alexandercoppock/git_projects/estimatr/revdep/library.noindex/estimatr/old/Rcpp/include' -I'/Users/alexandercoppock/git_projects/estimatr/revdep/library.noindex/estimatr/old/RcppEigen/include' -I/opt/R/arm64/include     -fPIC  -falign-functions=64 -Wall -g -O2   -c RcppExports.cpp -o RcppExports.o
In file included from RcppExports.cpp:4:
In file included from /Users/alexandercoppock/git_projects/estimatr/revdep/library.noindex/estimatr/old/RcppEigen/include/RcppEigen.h:25:
...
      |       ^
5 warnings generated.
clang++ -arch arm64 -std=gnu++20 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -L/Library/Frameworks/R.framework/Resources/lib -L/opt/R/arm64/lib -o hbal.so RcppExports.o eigen.o -L/Library/Frameworks/R.framework/Resources/lib -lRlapack -L/Library/Frameworks/R.framework/Resources/lib -lRblas -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0 -L/opt/gfortran/lib -lemutls_w -lheapt_w -lgfortran -lquadmath -F/Library/Frameworks/R.framework/.. -framework R
ld: warning: search path '/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0' not found
ld: warning: search path '/opt/gfortran/lib' not found
ld: library 'emutls_w' not found
clang++: error: linker command failed with exit code 1 (use -v to see invocation)
make: *** [hbal.so] Error 1
ERROR: compilation failed for package ‘hbal’
* removing ‘/Users/alexandercoppock/git_projects/estimatr/revdep/checks.noindex/hbal/old/hbal.Rcheck/hbal’


```
