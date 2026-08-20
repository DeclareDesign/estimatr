<!-- NOT READY TO SUBMIT. Delete this block once every line below is true.
     Outstanding as of 2026-08-19:
       - Graeme Blair's maintainer-transfer email to CRAN has not been sent.
       - The GitHub Actions matrix has not been run against the renamed package.
       - win-builder has not been run.
       - revdepcheck has not been run.
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
  (release), windows-latest (release)
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
All 35 were checked with `revdepcheck::revdep_check()`.
