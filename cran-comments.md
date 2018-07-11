
The package passes on Travis-CI, Appveyor, and winbuilder for the following configurations with only package size NOTEs:
- Ubuntu 14.04, R-release
- Ubuntu 14.04, R-devel
- Ubuntu 14.04, R-oldrel
- Mac OSX, R-release
- Mac OSX, R-oldrel
- Windows, R-release
- Windows, R-oldrel

There are, however, errors when the package is checked with `rchk` on R-hub. These errors appear to be internal to Rcpp or RcppEigen or in how they export Rcpp files and are largely out of our hands.
