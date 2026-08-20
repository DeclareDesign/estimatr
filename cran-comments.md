<!-- NOT READY TO SUBMIT. Delete this block once every line below is true.
     Outstanding as of 2026-08-20:
       - Graeme Blair's maintainer-transfer email to CRAN has not been sent.
       - win-builder has not been run.
       - The eventstudyr maintainer has not actually been emailed.
       - CI has not been re-run since the revdep fixes; the counts below are
         from the local suite only.
     Every claim in the sections below is written for the submitted state, so
     do not send this file while this block is still here. -->

## Submission

estimatr 2.0.0 is a rewrite of the package. Interfaces are unchanged from 1.0.6
and the numerical results agree to 1e-12 wherever both versions answer; the
breaking changes are listed in NEWS.md and in `vignette("estimatr2.0")`.

**This submission changes the maintainer** from Graeme Blair
<graeme.blair@gmail.com> to Alexander Coppock <acoppock@gmail.com>. Graeme Blair
has written to CRAN separately to confirm the transfer. He remains an author.

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
All 35 were checked. 34 are clean, `hbal` among them. `hbal` initially failed
to install on the test machine for a local toolchain reason unrelated to this
package: the machine has no `/opt/gfortran`, which is the path R's `FLIBS`
hardcodes. With that corrected it installs and checks `OK`.

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

Two defaults do change, in both cases away from an estimator that would require
expanding the fixed effects into dummies. Absorbed fixed effects with clusters
default to CR0 where 1.0.6 defaulted to CR2, and two or more absorbed factors
default to HC1 where 1.0.6 defaulted to HC2. Both emit a warning, once per session, that names the
1.0.6 default and names the `se_type` that accepts the new one. Naming `se_type`
explicitly returns the 1.0.6 number exactly.
