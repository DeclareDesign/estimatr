# estimatr 2.0

**This branch is the 2.0.0 release candidate.** It is a ground-up rewrite of estimatr aimed at the DeclareDesign workflow: fit the same model thousands of times, as fast as possible, and get the same number every time. `main` still holds 1.0.6, which is what CRAN serves.

```r
remotes::install_github("DeclareDesign/estimatr@rewrite", build_vignettes = TRUE)
vignette("estimatr2.0")
```

**Installing this replaces your CRAN estimatr**, since the package now carries its own name. Reinstall from CRAN to go back.

The vignette is the document to read first. It covers what does not change, what changes and why, how to port a 1.x script, and the benchmarks.

## What it is

`lm_robust`, `lm_lin`, `iv_robust`, `lh_robust`, `difference_in_means` and `horvitz_thompson`, in roughly half the lines, with a Pashley and Miratrix (2021) blocked-variance estimator that 1.x does not have.

Interfaces match 1.0.6 and, run side by side on one machine, the numbers match to 1e-12 across every supported standard error type, weighted and unweighted, clustered and unclustered, single and multivariate outcomes. The suite is 1,729 assertions, of which 336 run against answers recorded from an installed estimatr 1.0.6 (`data-raw/make_estimatr_reference.R` produces the recording) at a tolerance of 1e-9, because a fixture recorded on one platform meets a different BLAS on another and the floor there is the linear algebra rather than this package.

## What breaks

Two removals and one default, all deliberate, all covered in the vignette's porting section.

`horvitz_thompson()` takes one probability argument, `condition_prs`, in place of five. `blocks`, `clusters`, `simple`, `ra_declaration`, `condition_pr_mat`, `subset` and `return_condition_pr_mat` are gone, along with `se_type = "constant"` and the three exported matrix builders that served them (`declaration_to_condition_pr_mat`, `gen_pr_matrix_cluster`, `permutations_to_condition_pr_mat`). No design is lost: blocked, clustered and custom designs reach the estimator through an `ra_declaration` passed as `condition_prs`, and that path uses exact design-aware joint probabilities rather than the conservative bound.

`commarobust()` and `starprep()` are removed. Both still exist as names that error and name their replacement, rather than failing with "could not find function".

**Fixed effects without clusters support HC2 and HC3 exactly, at any number of factors, and default to HC2, which is what 1.x returns for the same call**, so that code ports with no change in the numbers at all. The projection decomposes exactly, so no dummy matrix is built to get there. CR2 is the exception: it is built from cluster-level blocks of the hat matrix rather than from the leverage diagonal, so it still expands the dummies, and clustered `fixed_effects` therefore defaults to CR0 where 1.x defaulted to CR2. That is the only default in the release that moves, and it warns once per session rather than moving silently. Writing `se_type = "CR2"` gets the 1.x number back exactly.

A clustered block holding a single treated or single control cluster is refused rather than given a variance. That is a correctness fix: 1.x returns a number there that is too small by roughly the block's cluster count.

## Status

`R CMD check --as-cran` clean, 1 NOTE (the maintainer change). Test suite 1,729 assertions, 0 failures, green on five platforms with identical counts on each.

A ledger of all 71 open estimatr issues against this implementation is in `notes/` (private): 26 fixed here, 7 out of scope, 6 not reproducible, 5 superseded by the rewrite, 4 still open.

Sibling branches: `DeclareDesign/fabricatr@rewrite`, `DeclareDesign/DeclareDesign@rewrite` and `DeclareDesign/randomizr@rewrite`.

## How this was written

estimatr 2.0.0 was written by Alexander Coppock working with Claude (Anthropic), across design, implementation, tests, benchmarks, and documentation. The evidence that should decide whether you install it is in `vignette("estimatr2.0")` under "How this was checked", and in `NEWS.md`: 336 assertions against a recording of an installed 1.0.6, and 305 more against implementations that share no lineage with estimatr, including `sandwich`, `clubSandwich`, `ivreg`, Stata, `fixest`, `plm` and `blkvar`.
