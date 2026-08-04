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

**It is not yet a drop-in replacement.** Measured against estimatr 1.0.4, four things are missing and three behaviours differ.

Absent entirely: `commarobust` and `starprep`. Written but never registered: `extract.lm_robust` and `extract.iv_robust`, whose implementation is sitting in `R/S3_extract.R` as `extract.robust_default`. Absent along with the argument they serve: `declaration_to_condition_pr_mat`, `gen_pr_matrix_cluster` and `permutations_to_condition_pr_mat`, which build the `condition_pr_mat` that `horvitz_thompson()` no longer accepts. And on a `horvitz_thompson` fit, `confint()`, `vcov()` and `glance()` error where estimatr answers, which is a defect rather than a design choice.

`horvitz_thompson()` also has a narrower signature: `blocks`, `clusters`, `simple`, `ra_declaration`, `subset` and `return_condition_pr_mat` are gone, and `se_type` no longer takes `"constant"`. The capability is not lost, since blocked, clustered and custom designs all reach the estimator through an `ra_declaration` passed as `condition_prs`, and that path uses exact design-aware joint probabilities. It is a deliberate narrowing of the interface rather than an omission, and it is the one item here that is a question for the group rather than a gap to fill.

Three behaviours move numbers or refuse where estimatr answered: the default under `fixed_effects` is HC1 rather than HC2, HC2/HC3/CR2 with `fixed_effects` is an error rather than a slow computation, and a clustered block with a singleton arm is refused rather than given a variance. That last one is a correctness fix; estimatr returns a number there that is too small by roughly the block's cluster count.

## Status

As of 2026-07-30: `R CMD check` 0 errors / 0 warnings / 0 notes, identity to estimatr holding at 1e-12.

A ledger of all 71 open estimatr issues against this implementation is in `notes/` (private) and mirrored to Alex's `claude_control/logs/estimatrZero_issue_ledger.md`: 26 fixed here, 7 out of scope, 6 not reproducible, 5 superseded by the rewrite.

**This branch is versioned 2.0.0**, against estimatr 1.0.6 on CRAN, and heads to CRAN with its three siblings. It is still *named* `estimatrZero`, and that is deliberate: the version says where it is going, the name is what keeps it installable alongside the released package. That matters more here than anywhere else, since the test suite checks this implementation against estimatr's to 1e-12 and both must load at once. The rename is the last step before release.

Sibling branches: `DeclareDesign/fabricatr@rewrite` and `DeclareDesign/DeclareDesign@rewrite`.
