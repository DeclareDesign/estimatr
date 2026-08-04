# estimatr, rewritten

**This branch is not the CRAN package and it is not `main`.** It holds a lean rewrite of estimatr targeting the DeclareDesign workflow. The released estimatr is unaffected by anything here.

The package on this branch is still named `estimatrZero`, so installing it leaves your CRAN estimatr in place and both can be loaded in the same session. That matters more here than for the sibling branches: the test suite checks this implementation against estimatr's to 1e-12, so the two must coexist.

```r
remotes::install_github("DeclareDesign/estimatr@rewrite", build_vignettes = TRUE)
vignette("estimatr2.0")
```

The vignette is the document to read first. It covers what does not change (with the live identity check against estimatr), what changes and why, how to port an existing script, and the benchmarks.

## What it is

4,007 lines of R against estimatr's 5,964. `lm_robust`, `lm_lin`, `iv_robust`, `lh_robust`, `difference_in_means` and `horvitz_thompson`, with a Pashley and Miratrix (2021) blocked-variance estimator that estimatr does not have.

**It is not a drop-in replacement, and the differences are the reason this branch is not on the same schedule as the other two.** Seven of estimatr's exports are absent: `commarobust`, `starprep`, `extract.lm_robust`, `extract.iv_robust`, `declaration_to_condition_pr_mat`, `gen_pr_matrix_cluster`, `permutations_to_condition_pr_mat`. Three behaviors move numbers or refuse where estimatr answered: the default under `fixed_effects` is HC1 rather than HC2, HC2/HC3/CR2 with `fixed_effects` is an error rather than a slow computation, and a clustered block with a singleton arm is refused rather than given a variance. That last one is a correctness fix; estimatr returns a number there that is too small by roughly the block's cluster count.

## Status

As of 2026-07-30: `R CMD check` 0 errors / 0 warnings / 0 notes, identity to estimatr holding at 1e-12.

A ledger of all 71 open estimatr issues against this implementation is in `notes/` (private) and mirrored to Alex's `claude_control/logs/estimatrZero_issue_ledger.md`: 26 fixed here, 7 out of scope, 6 not reproducible, 5 superseded by the rewrite.

The near-term plan is to back-port the fixes into estimatr 1.x as ordinary patches rather than to renumber this as 2.0.0. Nothing on this branch asserts a version: the DESCRIPTION still reads `estimatrZero 0.1.0`.

Sibling branches: `DeclareDesign/fabricatr@rewrite` and `DeclareDesign/DeclareDesign@rewrite`.
