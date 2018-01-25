This package 'Suggests' the fabricatr package which was recently added to CRAN. As a result, none of the binaries are available yet and causes errors in environments where packages must be installed from binary, such as R-release on win-builder.r-project.org.

One other NOTE may arise about possible misspellings in DESCRIPTION:

```
Possibly mis-spelled words in DESCRIPTION:
  Horvitz (12:285)
  pre (12:233)
```

Horvitz-Thompson is the commonly used name of an estimator that we implement, and pre-treatment, in its hyphenated form, is a commonly used phrase in the statistical literature.
