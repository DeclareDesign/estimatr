<!-- NOT READY TO SUBMIT. Delete this block once every line below is true.
     Outstanding as of 2026-08-20:
       - Graeme Blair's maintainer-transfer email to CRAN has not been sent.
       - win-builder has not been run.
       - The eventstudyr maintainer has not actually been emailed.
       - CI is green on all five platforms as of 245296d, but has NOT run
         against the fixed-effects leverage work that follows it.
       - The reverse-dependency counts below were measured against the final
         code on 2026-08-20 and are current.
     Every claim in the sections below is written for the submitted state, so
     do not send this file while this block is still here. -->

## Submission

estimatr 2.0.0 is a rewrite of the package. The six estimators keep their
signatures, with one exception: `horvitz_thompson()`, whose five probability
arguments consolidate into `condition_prs`. Numerical results agree to 1e-12
wherever both versions answer. The removals and the breaking changes are listed
in NEWS.md and in `vignette("estimatr2.0")`.

**This submission changes the maintainer** from Graeme Blair
<graeme.blair@gmail.com> to Alexander Coppock <acoppock@gmail.com>. Graeme Blair
has written to CRAN separately to confirm the transfer. He remains an author.

This version was written by the maintainer working with Claude (Anthropic).
`NEWS.md` says so and sets out the evidence: signatures unchanged except
`horvitz_thompson()`'s probability arguments, with the numbers agreeing to 1e-12
same-machine, a suite of 1,873 assertions
green on five platforms, 336 of them against answers recorded from an installed
1.0.6, 343 against independent implementations (`sandwich`, `clubSandwich`,
`ivreg`, Stata, `fixest`, `plm`, `blkvar`), and the full returned surface of sixteen fit
types pinned by test. All 35 reverse dependencies were checked.

## Test environments

* local macOS 15 (aarch64), R 4.6.0
* GitHub Actions: ubuntu-latest (devel, release, oldrel-1), macOS-latest
  (release), windows-latest (release). All five green, each reporting the same
  counts; the single skip is a test guarded on `blkvar`, which is not on CRAN.
* win-builder (devel and release)

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

The maintainer change is intentional and is covered by Graeme Blair's separate
message to CRAN.

## Reverse dependencies

estimatr has 17 strong reverse dependencies and 18 that suggest or enhance it.
All 35 were checked with `revdepcheck` against the submitted code; none failed
to check, and 34 are clean. (`hbal` needs one local correction to build at all
on the test machine, for a toolchain reason unrelated to this package: the
machine has no `/opt/gfortran`, which is the path R's `FLIBS` hardcodes. With
`FLIBS` pointed at the compiler that is installed, `hbal` installs and checks
`OK`.)

An earlier run of this release broke four. Three of those breaks were
regressions rather than intended changes: 1.0.6 computed HC2, HC3 and CR2 under
`fixed_effects` by expanding the fixed effects into dummies, and accepted a bare
grouping vector for that argument, and this release had dropped both. Both are
restored, and `clubSandwich`, `RCT` and `statuser` now check clean.

One remains. `eventstudyr` fails on three assertions that read a fitted object's
`felevels` list by the element name `V1`. estimatr 1.0.6 named those elements
after their terms except with a single fixed-effect factor on a model fitted
with missing data, where it fell back to `V1`; this release names them
consistently. Keeping 1.0.6's behaviour would mean keeping a bug that loses the
term name. No estimate changes: eventstudyr's one-way branch sizes its
degrees-of-freedom correction as `1 + rank` without consulting the field, and
its comparisons against recorded Stata output are unaffected. Its maintainer has
been notified.

One default changes, away from an estimator that would require expanding the
fixed effects into dummies. Absorbed fixed effects with clusters default to CR0
where 1.0.6 defaulted to CR2, and emit a warning, once per session, that names
the 1.0.6 default and names the `se_type` that accepts the new one. Naming
`se_type = "CR2"` returns the 1.0.6 number exactly. Unclustered absorbed fixed
effects keep 1.0.6's HC2 default at any number of factors: HC2 and HC3 no longer
need the dummies, so there is nothing to trade away.

One further difference is deliberate and affects only fixed-effect designs that
are rank deficient, where one absorbed factor is spanned by the others (a nested
factor, or a disconnected design). 1.0.6 sized the rank correction from the
nominal level count, so its absorbed fit disagreed with its own explicit-dummy
fit on those designs. This release takes the exact rank and the two now agree,
which is also what `lm()` and `plm` report for the same data. Designs of full
rank are unaffected.
