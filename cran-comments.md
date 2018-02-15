NOTE: This submission is identical to one roughly 20 hours ago, although this version fixes a small bug that we found in the interim. The main point of this submission remains the same as the one 20 hours ago.

The package is being resubmitted in response to several problems and errors found in the CRAN checks.

1) the package had tests ERROR on Solaris. These tests have been rewritten with greater tolerance and without reliance on small, fixed numeric values. Unfortunately, we were unable to check the new tests ourselves as the Solaris platform on r-hub was failing to compile RcppEigen.

2) there were tests that failed without long doubles enabled. These tests have been rewritten to only require equivalence and should pass.

3) Valgrind was erroring regarding some uninitialized values. All of the tests now pass without valgrind errors on Linux system running CentOS release 6.9 and R 3.4.2. Unfortunately, some of the tests were quite slow and running R CMD check on the package with "--use-valgrind" was not possible on r-hub without a timeout.

Other minor changes have been made, but this is largely just a patch for these issues which we were warned about. As we were told to resubmit by March 10, we wanted to try to resolve them ASAP and resubmit ASAP.

Also, as before, the NOTEs regarding misspellings in the DESCRIPTION file are false-positives.
