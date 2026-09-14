# estimatr 2.0

A ground-up rewrite of estimatr aimed at the DeclareDesign workflow: fit the same model thousands of times, as fast as possible, and get the same number every time.

```r
install.packages("estimatr")
vignette("estimatr2.0")
```

The second line opens the vignette in the help pane or a browser. It is the document to read first if you are coming from 1.x: what does not change, what changes and why, and how to port a 1.x script.

Three vignettes ship with the package, and three more (performance, the tidyverse, and regression tables) live on the [website](https://declaredesign.org/r/estimatr/) alone, because their numbers and their examples drift faster than a release cycle:

| vignette | what it is for |
|---|---|
| `estimatr2.0` | porting from 1.x: what changed, why, and the benchmarks |
| `getting-started` | the six estimators, on one worked example |
| `mathematical-notes` | every estimator defined, cited, and then validated against its own definition |

## What it is

`lm_robust`, `lm_lin`, `iv_robust`, `lh_robust`, `difference_in_means` and `horvitz_thompson`, with a Pashley and Miratrix (2021) blocked-variance estimator that 1.x does not have.

The six estimators keep their 1.0.6 signatures, except `horvitz_thompson()`, whose five probability arguments consolidate into `condition_prs`; the removals are below. Run side by side on one machine, 2.0 and 1.0.6 return the same estimates, standard errors, and degrees of freedom to 1e-12, across every supported standard error type, weighted and unweighted, clustered and unclustered, single and multivariate outcomes. The test suite holds that claim at 695 assertions against answers recorded from an installed 1.0.6, with `data-raw/make_estimatr_reference.R` producing the recording. Those run at 1e-9 rather than 1e-12, because a fixture recorded on one platform meets a different BLAS on another, and the floor there is the linear algebra rather than this package.

## What breaks

Two removals and one default, all deliberate, all covered in the vignette's porting section.

`horvitz_thompson()` takes one probability argument, `condition_prs`, in place of five. `blocks`, `clusters`, `simple`, `ra_declaration`, `condition_pr_mat`, `subset` and `return_condition_pr_mat` are gone, along with `se_type = "constant"` and the three exported matrix builders that served them (`declaration_to_condition_pr_mat`, `gen_pr_matrix_cluster`, `permutations_to_condition_pr_mat`). No design is lost: blocked, clustered and custom designs reach the estimator through an `ra_declaration` passed as `condition_prs`, and that path uses exact design-aware joint probabilities rather than the conservative bound.

`commarobust()` and `starprep()` are removed. Both still exist as names that error and name their replacement, rather than failing with "could not find function".

**Fixed effects without clusters support HC2 and HC3 exactly, at any number of factors, and default to HC2, which is what 1.x returns for the same call**, so that code ports with no change in the numbers at all. The projection decomposes exactly, so no dummy matrix is built to get there. CR2 is the exception: it is built from cluster-level blocks of the hat matrix rather than from the leverage diagonal, so it still expands the dummies, and clustered `fixed_effects` therefore defaults to CR0 where 1.x defaulted to CR2. That is the only default in the release that moves, and it warns once per session rather than moving silently. Writing `se_type = "CR2"` gets the 1.x number back exactly.

A clustered block holding a single treated or single control cluster is refused rather than given a variance. That is a correctness fix: 1.x returns a number there that is too small by roughly the block's cluster count.

## Status

`R CMD check --as-cran`: 0 errors, 0 warnings, 1 NOTE (the maintainer change). Test suite 5,852 assertions under `R CMD check`, 0 failures, one skip (`blkvar`, worth 12 assertions, which CRAN does not serve).

Every open estimatr issue was read against this implementation: 26 are fixed here, 23 are feature requests, 7 are out of scope, 6 are not reproducible, 5 are superseded by the rewrite, and 4 remain open. `vignette("estimatr2.0")` names the four.

## How this was written

estimatr 2.0.0 was written with AI. `vignette("mathematical-notes")` says so in full and, in the same document, states each estimator's definition in mathematics and then, immediately underneath, measures estimatr against that definition to machine precision, computed when the vignette is built. Sixteen such checks, and the vignette refuses to build if any of them fails. That is the evidence that should decide whether you install it.
