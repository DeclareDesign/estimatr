This package 'Suggests' the fabricatr package which was recently added to CRAN. As a result, none of the binaries are available yet and causes errors in environments where packages must be installed from binary, such as R-release on win-builder.r-project.org.

Some of our other test environments have passed R CMD check with only one NOTE:

  * a local OS X install, R 3.4.3
  * ubuntu 14.04 (on travis-ci), R 3.4.2
  * ubuntu 14.04 (on travis-ci), R 3.3.3
  * win-builder (devel)

The only NOTE is:

```
checking installed package size ... NOTE
  installed size is 24.6Mb
  sub-directories of 1Mb or more:
    libs  24.3Mb
```
