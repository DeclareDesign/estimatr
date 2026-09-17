## Submission

estimatr 2.0.0 is a rewrite of the package. The six estimators keep their signatures, with one exception: `horvitz_thompson()`, whose five probability arguments consolidate into `condition_prs`. Numerical results agree to 1e-12 wherever both versions answer. `tidy()`, `glance()`, and `augment()` now return tibbles, as broom's methods do. The removals and the breaking changes are listed in NEWS.md and in `vignette("estimatr2.0")`.

**This submission changes the maintainer** from Graeme Blair <graeme.blair@gmail.com> to Alexander Coppock <acoppock@gmail.com>. Graeme Blair has written to CRAN separately to confirm the transfer. He remains an author. One other change to `Authors@R`: Macartan Humphreys has asked to be listed as a contributor rather than an author, and his role moves from `aut` to `ctb`.

This version was written by the maintainers working with AI assistance (Claude, from Anthropic). `vignette("estimatr2.0")` says so, and `vignette("mathematical-notes")` pairs each estimator's definition with a measurement of estimatr against that definition, to machine precision, computed when the vignette is built; the vignette ends in `stopifnot()`, so a broken identity fails the check rather than printing FALSE in a table. The evidence for this release: a suite of 5,906 assertions, 695 of them against answers recorded from an installed estimatr 1.0.6 and 808 against independent implementations (`sandwich`, `clubSandwich`, `ivreg`, Stata, `fixest`, `plm`, and `blkvar`), with the full returned surface of sixteen fit types pinned by test. Every reverse dependency was checked.

## Test environments

* local macOS 26.6 (aarch64, Apple M4), R 4.6.0
* GitHub Actions: ubuntu-latest (devel, release, oldrel-1), macOS-latest (release), windows-latest (release). All five green at `FAIL 0 | SKIP 1 | PASS 5906` (run 35003141386). R-devel alone reports two test warnings, both raised inside `clubSandwich::vcovCR()` on an `mlm` fit that a test uses as its reference value ("Replacing special names '.Dimnames' is deprecated"); they come from clubSandwich, not from estimatr. Up to two test blocks skip by design: the `blkvar` comparison, that package being available only from GitHub, and one test that needs randomizr 2.0.1, which is not yet on CRAN.
* win-builder: R-devel (2026-09-14 r90539 ucrt) and R 4.6.1 (ucrt), Windows x86_64, `Status: 1 NOTE` on each, the maintainer change and nothing else. Tests OK and all three vignettes rebuilt on both.

## R CMD check results

0 errors | 0 warnings | 1 note

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Alexander Coppock <acoppock@gmail.com>'

New maintainer:
  Alexander Coppock <acoppock@gmail.com>
Old maintainer(s):
  Graeme Blair <graeme.blair@gmail.com>
```

The maintainer change is intentional and is covered by Graeme Blair's separate message to CRAN.

## Reverse dependencies

`revdepcheck` was run against the submitted code on 2026-09-14: 37 checked, 35 clean, 2 broken, none failed to check. Both breaks were expected, and each maintainer has been notified with the fix.

`eventstudyr` fails four test assertions, and no estimate changes. Three read a fitted object's `felevels` list by the element name `V1`. estimatr 1.0.6 named those elements after their terms except with a single fixed-effect factor on a model fitted with missing data, where it fell back to `V1`; this release names them consistently. Keeping 1.0.6's behaviour would mean keeping a bug that loses the term name, which also broke CR2 and `predict()` on those fits. The fourth assertion checks the `dim` of a column assigned into `tidy()` output: a data frame stores a one-column matrix as a matrix column, and a tibble stores it as a vector, with identical values. Its maintainers were emailed on 2026-08-26 about the first and on 2026-09-15 about the second.

`hbal`'s examples raise a warning. `att()` sets row names on `tidy()` output, which tibble deprecates. The one-line fix on their side is to coerce with `as.data.frame()` before setting row names, and every value in its table is unchanged. Its maintainer was emailed on 2026-09-15 and replied the same day that they will make the change.

`projoint`, which an earlier run of this release broke on a single-cluster fit, checks clean at version 1.1.4. An earlier run also broke `clubSandwich`, `RCT`, and `statuser`; those were regressions in the rewrite, they are fixed, and all three check clean.

One default changes, away from an estimator that would require expanding the fixed effects into dummies. Absorbed fixed effects with clusters default to CR0 where 1.0.6 defaulted to CR2, and emit a warning, once per session, that names the 1.0.6 default and names the `se_type` that accepts the new one. Naming `se_type = "CR2"` returns the 1.0.6 number exactly. Unclustered absorbed fixed effects keep 1.0.6's HC2 default at any number of factors.

One further difference is deliberate and affects only fixed-effect designs that are rank deficient, where one absorbed factor is spanned by the others (a nested factor, or a disconnected design). 1.0.6 sized the rank correction from the nominal level count, so its absorbed fit disagreed with its own explicit-dummy fit on those designs. This release takes the exact rank and the two now agree, which is also what `lm()` and `plm` report for the same data. Designs of full rank are unaffected.
