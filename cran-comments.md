## Submission

estimatr 2.0.1 is a patch release. It fixes wrong answers and regressions found in 2.0.0 after it was released, nearly all of them on rank-deficient or otherwise degenerate designs. The ones that moved a number: `horvitz_thompson()` returned `NA` for the standard error on every draw whenever a clustered or blocked-and-clustered declaration's `clusters` was a factor, which 1.0.6 answered correctly and which this release restores to fifteen digits; `try_cholesky = TRUE` returned a full set of coefficients at a design of deficient rank, because the rank guard compared a Gram-factor diagonal against a tolerance set for `dqrdc2`; absorbed fixed-effect estimates came back all `NA` on a rank-deficient fit; a rank-deficient fit reported its F statistic as `NA` unless the dropped column happened to be last; and which collinear column is dropped now follows `stats::lm()` rather than a pivoted QR's norm ordering, so a rank-deficient fit returns `NA` for the same coefficient `lm()` does. Separately, a fit no longer draws from the random number stream, so a seeded simulation that contains one reproduces across versions. Two of the fourteen will reach a user as a difference rather than as a fix, and both are deliberate: `HC2` and `HC3` now return `NA` for a coefficient that observations at leverage 1 alone identify, where 2.0 returned a number assembled entirely from rows carrying no information about it; and the notice that a collinear column was dropped is now a message rather than a warning, which matches `stats::lm()` on the one axis a caller can intercept, so code that silenced it with `suppressWarnings()` needs `suppressMessages()`. `NEWS.md` lists all fourteen changes.

This version was written by the maintainers working with AI assistance (Claude, from Anthropic), as 2.0.0 was. `vignette("estimatr2.0")` says so, and `vignette("mathematical-notes")` pairs each estimator's definition with a measurement of estimatr against that definition, to machine precision, computed when the vignette is built; the vignette ends in `stopifnot()`, so a broken identity fails the check rather than printing FALSE in a table. The evidence for this release: a suite of 7,177 assertions, 695 of them against answers recorded from an installed estimatr 1.0.6 and 820 against independent implementations (`sandwich`, `clubSandwich`, `ivreg`, Stata, `fixest`, `plm`, and `blkvar`), with the full returned surface of sixteen fit types pinned by test.

## Test environments

* local macOS 26.6 (aarch64, Apple M4), R 4.6.0: `Status: OK`, and `FAIL 0 | WARN 0 | SKIP 1 | PASS 7177`.
* GitHub Actions: ubuntu-latest (devel, release, oldrel-1), macOS-latest (release), windows-latest (release). All five green at `FAIL 0 | SKIP 1 | PASS 7177` (run 36290783813), read per job. R-devel alone reports two test warnings, both raised inside `clubSandwich::vcovCR()` on an `mlm` fit that a test uses as its reference value ("Replacing special names '.Dimnames' is deprecated"); they come from clubSandwich, not from estimatr.
* One test block skips by design, worth 18 assertions: the check that `blkvar`'s recorded values are still what `blkvar` and randomizr produce, that package being available only from GitHub. The 12 assertions comparing estimatr against the recording read the fixture and do not skip, on CRAN or anywhere else.
* win-builder, both arms, x86_64-w64-mingw32 under g++ 14.3.0 on Windows Server 2022, run 2026-09-27: `Status: OK` on each, 0 errors, 0 warnings, 0 notes, with CRAN incoming feasibility OK.
    * release, R 4.6.1 (2026-06-24 ucrt): `FAIL 0 | WARN 0 | SKIP 1 | PASS 7177`.
    * R-devel, R Under development (2026-09-25 r90590 ucrt): `FAIL 0 | WARN 2 | SKIP 1 | PASS 7177`. The two warnings are the clubSandwich ones described above, and they appear on R-devel alone here as they do on GitHub Actions.
    * The skip is the expected one on both, reported as `{blkvar} is not installed (1): 'test_blocked_variance.R:65:3'`. So 7,177 is what `R CMD check` prints on all seven platforms this release was checked on.

## R CMD check results

0 errors | 0 warnings | 0 notes

CRAN incoming feasibility passes with no note. The `2.0.0.9000` version-component note the 2.0.0 submission carried is gone with the version bump, and the maintainer change that release announced is already in place.

## Reverse dependencies

Re-run on 2026-09-27 against the commit this tarball is built from, with the baseline at CRAN 2.0.0: 38 checked, 37 OK, 1 new problem, 0 failed to check. The check library, the check directory, and the result database were deleted before the run, so no package's result was carried over from the earlier run of 2026-09-26, which had been measured against an earlier commit.

`DesignLibrary` 0.1.10 is the one new problem, and it is ours. `difference_in_means()` now returns a `call` field, which every other estimator in this package already had and which was added because `difference_in_means()` was the only one missing it. `DesignLibrary`'s `multi_arm_designer()` combines the fitted objects themselves with `rbind.data.frame()`, which requires every field of a fit to be of length 1, so it now fails with "invalid list argument: all variables should have the same length". Three of its test assertions fail and no estimate it reports changes.

We maintain `DesignLibrary`, and the fix is already written rather than planned: its 2.0.0 delegates that designer to an internal library call and contains no `rbind.data.frame()` anywhere. That release depends on the `fabricatr` and `DeclareDesign` 2.0 releases and is third in that sequence, so it will reach CRAN after this one. We are reporting the break rather than patching the 0.1.10 line because the successor removes the pattern entirely, and we did not want to ship a fix to a file that is about to be replaced.

`eventstudyr` reports an error under both the old and the new version, so it is not a new problem in this release. Its situation is the one described in the 2.0.0 submission: three assertions read a fitted object's `felevels` list by the element name `V1`, which 1.0.6 produced only in the case where it lost the term name, and a fourth checks the `dim` of a column assigned into `tidy()` output, which a tibble stores as a vector where a data frame stored a one-column matrix. Its maintainers were emailed on 2026-08-26 and on 2026-09-15.

`hbal`, whose examples raised a tibble row-names warning at 2.0.0, checks clean; its maintainer made the change they said they would.
