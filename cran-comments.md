NOTE: This submission is identical to one roughly 20 hours ago, although this version fixes a small bug that we found in the interim. The main point of this submission remains the same as the one 20 hours ago.

The package is being resubmitted to fix some bugs and add functionality, such as multivariate linear regression and two stage least squares instrumental variables regression.

The package pasts on Travis-CI, Appveyor, and winbuilder for the following configurations with only package size NOTEs:
- Ubuntu 14.04, R-release
- Ubuntu 14.04, R-devel
- Ubuntu 14.04, R-oldrel
- Mac OSX, R-release
- Mac OSX, R-oldrel
- Windows, R-release
- Windows, R-oldrel

There are, however, errors when the package is checked with `rchk` on R-hub. These errors appear to be internal to Rcpp or RcppEigen or in how they export Rcpp files and are largely out of our hands.
