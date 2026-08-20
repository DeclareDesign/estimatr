# estimatr 2.0.0

estimatr 2.0.0 is a ground-up rewrite aimed at the DeclareDesign use case: OLS, Lin-adjusted OLS, 2SLS IV, difference-in-means, Horvitz-Thompson, and linear hypothesis tests with heteroskedasticity- and cluster-robust standard errors. It fixes several long-standing correctness bugs, improves performance on the critical path, adds feols-style fixed effects absorption, and replaces the O(N²) Horvitz-Thompson variance with a design-aware O(1) computation. Interfaces are unchanged from 1.0.6 and the numbers agree to 1e-12 wherever both versions answer.

See `vignette("estimatr2.0")` for a user-facing tour of what changes and what does not.

---

## Found by `revdepcheck`, and fixed here

Running the 35 CRAN reverse dependencies turned up six fields that had gone
missing from fitted objects, and one that was returned with a wrong value. None
of it was visible to a NAMESPACE diff, since no export changed, nor to the
estimator comparisons, which name the fields they check, nor to a grep of the
reverse dependencies' sources, since the calls are unchanged. It surfaced
because `eventstudyr` reads `felevels` on a fixed-effects fit and got `NULL`.

Restored: `felevels` and the unprojected `tss` on every fixed-effects fit;
`fixed_effects` (the absorbed group effects) and `proj_fstatistic` on
`iv_robust()` fixed-effects fits; `treatment_levels` on `lm_lin()` fits, which
holds every level including the baseline and is not the same thing as the
`treatment_vals` this version added; and `nclusters`, reported as `0`, on an
unclustered blocked `difference_in_means()`.

**`proj_fstatistic` was wrong, not merely missing.** With absorbed fixed
effects the intercept is already out of the design, and the numerator degrees
of freedom subtracted one anyway, so the statistic tested a hypothesis one
restriction short of the intended one. A two-regressor fixed-effects model
reported F on 1 numerator degree of freedom rather than 2, with a
correspondingly wrong value. A single-regressor fixed-effects model took the
count to zero and dropped the statistic altogether, which is how `iv_robust()`
came to have no `proj_fstatistic` at all.

`felevels` also came back in the wrong order. The fixed-effects matrix is
coerced to character before demeaning, which re-derives the levels by string
sort, so a factor with levels 1 to 30 reported them as 1, 10, 11. The names are
now taken before that coercion.

**`felevels` also had to keep listing the levels that are in the fit rather
than the levels that are in the data.** The first attempt at the ordering fix
above read the level names off the model frame before rows were dropped for
missingness, which counted levels that no observation survives into. That is
not cosmetic: downstream code sizes a degrees-of-freedom correction off
`length(fit$felevels[[term]])`. On `eventstudyr`'s example data, 1,700 complete
rows of 2,000 leave 34 of 40 time periods in the estimation sample, and
counting all 40 moved its standard errors in the fourth significant figure.
`felevels` now carries the estimation-sample levels, in the order the data
gave them.

One difference from 1.0.6 does remain, and it moves no estimate. 1.0.6 named
each element after its term, `~ get(idvar)` giving `felevels[["get(idvar)"]]`,
except with a single fixed-effect factor on a model fitted with missing data,
where it fell back to `V1`. 2.0 names it consistently. Code that reads
`felevels$V1` will need the term name instead.

`tests/testthat/test_return_surface.R` now pins the whole returned surface of
sixteen fit types against answers recorded from 1.0.6, so a field cannot go
missing again without a test failing. A missing field is the failure mode worth
guarding: `fit$felevels` returns `NULL` rather than erroring, so the reader gets
a wrong answer instead of a stop.

---

## Two more breaking changes, which had not been written down

Both were deliberate and neither was listed. Found the same way.

**`fixed_effects` must be a one-sided formula.** 1.x also accepted a bare
column name or an already-evaluated grouping vector. Refusing those is the
resolution of issue #304, which asked for the argument to be enforced, and the
error now says so. `RCT` relies on the old form.

**HC2, HC3 and CR2 with `fixed_effects`**: see below. `clubSandwich` and
`statuser` both call for CR2 or HC3 on an absorbed model.

### What the reverse-dependency run leaves broken

Of estimatr's 35 CRAN reverse dependencies, 30 check clean against 2.0.0 and
one (`hbal`) fails to install on the test machine under both versions, for a
local toolchain reason unrelated to this package.

Four break. Three need a one-line fix at the call site:

* **`clubSandwich`** calls `iv_robust()` and `lm_robust()` with `fixed_effects`
  and `se_type = "CR2"`, in one example and three tests.
* **`statuser`** documents and defaults to HC3, and its examples absorb two
  fixed-effect factors.
* **`RCT`** passes an evaluated grouping vector to `fixed_effects` on its
  single-factor path.

* **`eventstudyr`** has one assertion that reads `felevels$V1` on a one-way
  fixed-effects fit and wants the 1.0.6 fallback name. Nothing else in its
  suite changes, and no estimate moves: its one-way branch sizes its
  degrees-of-freedom correction as `1 + rank` without consulting `felevels`.

---

## What is kept

All public functions from 1.x that are relevant to the DeclareDesign workflow are present with identical interfaces:

- `lm_robust()` — OLS with HC0/HC1/HC2/HC3/CR0/CR2/stata SEs
- `lm_lin()` — Lin (2013) covariate-adjusted estimator
- `iv_robust()` — two-stage least squares with the same SE menu
- `difference_in_means()` — Neyman variance, paired, blocked, and clustered designs
- `horvitz_thompson()` — inverse probability weighting with Young's inequality variance
- `lh_robust()` — linear hypothesis tests via `car::linearHypothesis` with robust variance
- S3 methods: `tidy`, `glance`, `summary`, `print`, `predict`, `coef`, `confint`, `vcov`, `nobs`, `update`

Return objects are structurally compatible with estimatr 1.0.6: field names, classes, and method dispatch are unchanged. Cross-package identity tests confirm coefficient and standard error agreement to 1e-12 across all supported SE types.

`horvitz_thompson()` fits now carry the same post-estimation methods as every other estimator here. `vcov()`, `confint()`, `glance()`, `summary()` and `nobs()` were missing, so the first three errored and the last two fell through to the base defaults and returned something that resembled output without being the estimate. All five now agree with estimatr 1.0.6 exactly, including `confint(level = )` at a non-default level, where the interval is normal rather than t because the estimator has no degrees of freedom to spend. The `glance()` gap was the one worth finding: without it `modelsummary()` did not fail, it silently dropped the goodness-of-fit rows and printed a coefficient table that looked complete.

---

## HC2 and HC3 with fixed effects

Absorbed fixed effects used to refuse `se_type = "HC2"` and `"HC3"`, on the grounds that they need hat values from the full `[dummies | X]` design while absorption leaves only the demeaned ones. That is true for two or more FE factors and false for one, where the hat value splits exactly:

```
h_ii = h_ii(demeaned X) + w_i / sum(w over i's group)
```

So with a single FE factor both are now exact, and cost one vector rather than a dummy hat matrix. Verified to machine precision against writing the dummies out, and against 1.x, weighted and unweighted, balanced and unbalanced. At n = 40,000 with 2,000 groups, HC2 takes 17.8 ms here against 41.3 s in 1.0.6, with standard errors bit-identical.

**The default with one-way fixed effects therefore returns to `"HC2"`**, which is the package default everywhere else and is also what 1.x returns for the same call. It had fallen back to `"stata"` only because HC2 was unaffordable.

Two or more FE factors keep the restriction: the analogous sum is wrong in the third decimal place, so `"HC2"` and `"HC3"` still error there and name the dummy formulation. `"CR2"` is still refused with any `fixed_effects`, since its adjustment is built from cluster-level blocks of the hat matrix rather than from `h_ii`, and the identity says nothing about it. `iv_robust()` keeps the restriction too, since the identity has not been derived for a second stage running on fitted regressors.

Fixing this exposed a latent crash: `r_fe` was left at `r + fe_rank` when the meat falls back to the plain `XtX_inv`, so `X.leftCols(r_fe)` overran and aborted inside Eigen. Unreachable while HC2 was refused with fixed effects, reachable now.

---

## What is dropped

**Horvitz-Thompson auxiliary arguments.** `horvitz_thompson()` is included, but the probability specification is consolidated into the single `condition_prs` argument. The `blocks`, `clusters`, `simple`, `condition_pr_mat`, `subset`, and `return_condition_pr_mat` arguments and `se_type = "constant"` are dropped. An `ra_declaration` already carries the block structure, cluster structure, per-unit probabilities, and simple-versus-complete flag, so the separate arguments restated information the declaration holds, and gave the estimator two sources of truth that could disagree. Fully custom designs go through `declare_ra(permutation_matrix = perm)`, which replaces `condition_pr_mat`.

**`starprep()` and `commarobust()`.** Both are removed, but they remain as names that error with a pointer to the replacement rather than being deleted, so a script written against 1.x reports what happened instead of failing with `could not find function`. `commarobust()` recomputed robust standard errors on a fitted `lm`, which is what `lm_robust()` does. `starprep()` prepared fits for stargazer, which has not been maintained in years; table building now goes through modelsummary, which reads `tidy()` and `glance()` and therefore works on every estimator here with no adapter.

The condition-probability-matrix helpers `declaration_to_condition_pr_mat()`, `gen_pr_matrix_cluster()`, and `permutations_to_condition_pr_mat()` go with the Horvitz-Thompson arguments above: the new variance computation does not build that matrix.

**`extract.lm_robust()` and `extract.iv_robust()` are kept**, for \pkg{texreg}, which is their only consumer. They are exported as plain functions rather than registered as S3 methods because texreg looks them up by name.

**HAC (Newey-West) standard errors.** Not implemented.

**CR3 standard errors.** Not implemented.

---

## Bug fixes

### Weighted R² formula: investigated, no change (affects `lm_robust`, `lm_lin`, `iv_robust`)

An early version of 2.0 "corrected" the weighted TSS from `weights_internal²` to `weights_internal`. That change was wrong and was reverted. estimatr 1.x's weighted R² formula is correct, and 2.0 reproduces it exactly. The reasoning is recorded here so it is not re-litigated.

The standard weighted TSS is:

```
TSS = Σᵢ wᵢ (yᵢ - ȳ_w)²
```

where ȳ_w is the standard weighted mean (weights wᵢ). Inside `lm_robust_fit`, raw weights are rescaled to `weights_internal = sqrt(w / μ_w)`, so `weights_internal² = w / μ_w`. Using `weights_internal²` in the TSS formula recovers the correct standard weighted mean as the centering point. Using `weights_internal` instead (with sqrt-weights) centers on ȳ_{sqrt(w)} ≠ ȳ_w, which breaks the WLS guarantee that RSS ≤ TSS and allows R² < 0.

2.0 uses `weights_internal²` throughout, matching estimatr 1.x's formula exactly. For an intercept-only WLS model, R² = 0 exactly. For any WLS model, R² ∈ [0, 1].

### iv_robust diagnostic tests (affects `iv_robust` with `diagnostics = TRUE`)

Two bugs in the diagnostic computations were identified and fixed:

**First-stage F-test missing p-value.** The returned vector from `first_stage_ftest()` lacked a `p.value` entry, causing `glance()` on `iv_robust` objects with diagnostics to silently drop the p-value or error downstream.

**Wu-Hausman endogeneity test wrong numerator df.** `wu_hausman_reg_ftest()` computed `numdf` as the number of columns in the first-stage residual matrix. This over-counted by one because the intercept residuals are collinear and are automatically dropped by `lm_robust_fit`, reducing the effective rank by one relative to the column count. The fix computes `numdf = lm_resids$rank - lm_noresids$rank`, which reflects the actual rank increase from adding the residuals.

### Ordered factor cluster variable crashes (#421)

Passing an ordered factor as the `clusters` argument caused the error `'length = 2' in coercion to 'logical(1)'`. The class check `class(x) %in% c("factor", "integer")` returns a length-2 logical vector for ordered factors (whose class is `c("ordered", "factor")`), which R's `if()` rejects. Fixed by replacing the class check with `is.factor()` and `is.integer()`.

### `formula()` object passed to `fixed_effects` crashes (#348)

When a pre-built formula object is stored in a variable and passed as `fixed_effects = fe_formula`, the rlang quosure captures the variable name (a symbol), not the formula value. Passing the raw quosure to `stats::model.frame.default()` produced the error `invalid type (language)`. Fixed by calling `eval_tidy()` on the quosure before handing it to `model.frame.default()`.

### Intercept-only model with fixed effects crashes (#303)

`lm_robust(y ~ 1, fixed_effects = ~block)` crashed with `length of 'dimnames' [2] not equal to array extent`. After FE demeaning the intercept column is dropped, leaving a zero-column design matrix. The alternating-projections loop in `demean_fes()` ran 50 iterations on the empty matrix (the convergence check `max(abs(empty_matrix))` returns `-Inf`, which never satisfies the criterion), accumulating 50+ warnings before a downstream dimension error.

The fix has two parts. First, `demean_fes()` skips the demeaning loop entirely when there are no non-intercept columns. Second, `lm_robust()` detects the zero-column case after demeaning and returns a well-formed result directly: empty coefficient vector, correct `df.residual = n - fe_rank`, full-model R², and fitted values equal to the within-group means.

### `lh_robust` CIs inconsistent with `lm_robust` when clusters present (#405)

The individual hypothesis CIs and p-values from `lh_robust()` did not match those from `lm_robust()` when clusters were specified. The root cause: `lh_robust()` used `lmr$df.residual` (= n − k, e.g., 97 for n = 100, k = 3) as the degrees of freedom for t-tests and confidence intervals. For clustered models, the correct df is the cluster-adjusted value stored per-coefficient in `lmr$df` (e.g., 9 for a CR2 model with 10 clusters, or 8–9 via Satterthwaite approximation).

Using the residual df of 97 instead of the cluster df of ~9 produces confidence intervals that are far too narrow. The fix parses the hypothesis expression (e.g., `"x=0"`) to extract the coefficient name (`"x"`), looks up `lmr$df["x"]`, and uses that as the df. For composite hypotheses involving multiple coefficients (e.g., `"x - z=0"`), the conservative minimum over all per-coefficient dfs is used.

### `lh_robust` joint hypothesis test absent (#320)

`lh_robust()` returned only individual hypothesis tests with no joint F-test. The fix adds a `joint_hypothesis` element to the return object containing a Wald F-statistic:

```
F = (Lβ̂)ᵀ (L Vcov Lᵀ)⁻¹ (Lβ̂) / m
```

where m is the number of restrictions, tested against F(m, df) using the same cluster-adjusted df as the individual tests. The vcov matrix `L Vcov Lᵀ` is obtained directly from `car::linearHypothesis()` via the `vcov` attribute of its return value, so the cluster-robust variance is correctly propagated.

The same work closes #390, where `lh_robust(Y ~ X1 + X2 + X3, linear_hypothesis = c("X1", "X2"))` should return one joint test and 1.x instead errors: "lh_robust currently implements tests for hypotheses involving linear combinations of variables but not joint hypotheses."

### Residuals are returned (#345)

Neither `lm_robust()` nor `iv_robust()` nor `lm_lin()` returned residuals: `residuals(fit)` was `NULL`, and `fit$residuals` did not exist, though `fitted.values` did. `lm_return()` dropped the field.

2.0 returns `residuals` on the scale of the data. Weighted fits return the unweighted residuals, so `residuals + fitted.values` recovers the outcome. Clustered fits are returned in the original row order, undoing the internal sort by cluster. Fixed-effects models return full-model residuals, which by Frisch-Waugh-Lovell equal the residuals of the demeaned regression. `iv_robust()` returns structural residuals, y − Xβ̂ on the original endogenous regressors rather than their first-stage projections. Because `residuals.default()` reads `object$residuals`, `residuals(fit)` now works with no new method.

### Dropped collinear regressors warn (#411)

A regressor collinear with the others is dropped and returned as an NA coefficient, silently. The reporting user compared `lm()` and `lm_robust()` output, found coefficients that differed, and had no way to see that a term had been dropped and the rest were therefore conditional on a different set of regressors.

2.0 warns once, naming the dropped terms. Note that an under-identified `iv_robust()` call now raises two warnings, "More regressors than instruments" and the collinearity warning for the intercept it drops; both are accurate.

### `augment()` for use with broom-aware packages (#377)

`augment()` returns the model frame with `.fitted` and `.resid` appended, which is what downstream packages expect before they will accept a model object. It was not implementable before this release, since no estimator returned residuals (#345). `augment(newdata = )` predicts on new data instead. Multivariate outcomes are refused. Methods are provided for `lm_robust` (covering `lm_lin`) and `iv_robust`.

### Closed by the rewrite rather than by a fix

- **#233, Horvitz-Thompson messages.** estimatr 1.x prints `` `simple` = FALSE, using complete cluster randomization `` and `` `condition_prs` not found, estimating probability of treatment to be constant... `` on every call, which in a simulation means every iteration. 2.0 requires `condition_prs`, so it never guesses a probability and never prints either message.
- **#334, faster fixed effects.** The issue asks for `fixest`'s algorithm. 2.0 uses it: alternating projections in C++, 200x to 4,000x faster on the designs measured above.
- **#365, replace AER tests with ivreg**, and **#402, skip tests relying on suggested packages.** Both are true by construction here: the suite never used AER, and every test needing `randomizr`, `carData` or `blkvar` is guarded with `skip_if_not_installed()`. The comparisons against 1.0.6 need no guard at all, because they read recorded answers rather than calling the older version.
- **#260, bad defaulting to matched pairs.** The reporter has 17 blocks, one of which holds a single treated and single control unit, and estimatr 1.x collapses the whole design to matched pairs on 16 degrees of freedom. With two or more such blocks 2.0 now estimates the design properly; with exactly one, as here, it errors naming the block, because one small block leaves nothing to estimate its contribution against.

`R CMD check` runs clean: 0 errors, 0 warnings, 0 notes.

### Blocked designs with blocks of different sizes (#336)

1.x handles two kinds of blocked design and refuses everything between them. If every block has at least two treated and two control units, each block carries its own Neyman variance. If every block is a matched pair, the variance comes from the variation across pairs. A design with both, or with a block holding one treated and three control units, either errors ("every block must have at least two treated and control units") or silently applies the matched-pairs estimator to every block, big ones included, after a warning. Such designs are not exotic: coarsened exact matching, full matching, and multisite trials with one or two sites per stratum all produce them.

2.0 implements the estimators of Pashley and Miratrix (2021). Blocks are classified by how many units each **arm** holds, not by how large the block is, which is the substantive correction: a block of eight units with one of them treated has no more estimable within-block variance than a matched pair does.

- Blocks with at least two treated and two control units contribute their own Neyman variance, their equation 4.
- Blocks with a singleton treated or control unit contribute through the variation across such blocks, their equation 8, the "unified" estimator that requires no two blocks to share a size. With equal sizes this reduces to the usual matched-pairs estimator, their equation 5, which is used directly because equation 8 is undefined at two equal-sized blocks.
- A design with both kinds combines the parts by squared share of the sample, their section 3.3.

Degrees of freedom are not treated in the paper. 2.0 combines the two components by Welch-Satterthwaite, which reduces to `n - 2K` for an all-big design and to `K - 1` for an all-small one, matching what each literature uses on its own.

The `design` element now reports `"Blocked"`, `"Matched-pair"`, `"Small blocks"`, or `"Hybrid blocked"`. Standard errors agree to 1e-10 with `blkvar::block_estimator(method = "hybrid_p")`, the authors' own implementation, across all-big, all-small, matched-pair, and hybrid designs.

Two designs are refused, because the variance genuinely cannot be estimated rather than because the software is unwilling: exactly one block with a singleton arm, which leaves nothing to compare it against and would contribute a variance of zero, and a set of different-sized such blocks in which one holds half or more of their units, which is the condition equation 8 needs to stay defined and conservative.

Unchanged: all-big designs and matched pairs return exactly what they did before.

### Blocks of clusters with a singleton cluster in an arm are now refused

Pashley and Miratrix cover treatment assigned within blocks, not blocks of clusters: clusters appear once in their paper, to be set aside ("we are focusing on treatment assigned within cluster", section 2). The hybrid above therefore does not extend to cluster-randomized blocks, and blocked designs with `clusters` keep the previous estimators.

Checking that boundary turned up a separate defect, present in 1.x too. A block with one treated or one control **cluster** has no estimable within-block variance, because the singleton arm contributes nothing to the between-cluster variability, yet both packages returned a number for it with no warning. Exhaustive enumeration of every assignment, with cluster-level potential outcomes held fixed:

| block | ratio of mean estimated variance to true variance |
|---|---:|
| 6 clusters, 3 treated | 1.00 |
| 6 clusters, 2 treated | 1.00 |
| 4 clusters, 2 treated | 1.00 |
| 4 clusters, 1 treated | 0.25 |
| 6 clusters, 1 treated | 0.17 |
| 8 clusters, 1 treated | 0.12 |

With two or more clusters per arm the estimator is exact. With one it is too small by roughly the block's cluster count, and it worsens as the block grows. End to end, blocks of 6, 6, 6, and 4 clusters where the last has a single treated cluster gave 88% coverage of a nominal 95% interval.

2.0 now errors, naming the offending blocks and suggesting either merging them so every block has two clusters per arm, or `lm_robust()` with block fixed effects and `clusters`. Matched-pair clustered designs are exempt, since every block there has a singleton arm by construction and the variance is estimated across blocks rather than within them.

### Rank detection matches `lm()` (#351, #395)

`lm_robust(y ~ x)` with `x` constant returned coefficients of order 1e11 instead of NA. The design matrix is exactly rank deficient, but Eigen's default `ColPivHouseholderQR` threshold is roughly `epsilon * ncol` relative to the largest pivot, tight enough that the collinear column survives as a pivot of order 1e-14. `stats::lm()` uses LINPACK `dqrdc2` with `tol = 1e-7`; 2.0 now sets the same threshold, so rank detection agrees with `lm()` on degenerate designs. Cross-package agreement with estimatr 1.x is unaffected on any full-rank design.

A near-saturated design also produces observations with leverage exactly 1, which makes HC2, HC3, and CR2 divide by zero and return `NaN`. 2.0 now says so, naming the SE type and suggesting `"HC1"` or `"classical"`, rather than returning a bare `NaN`.

### `predict()` with fixed effects (#403, #404)

Three separate failures. `predict()` with no `newdata` errored instead of returning the in-sample fit. `predict()` on a fixed-effects model with a factor regressor failed with "non-conformable arguments", because the terms object still carries the intercept that fixed-effects absorption removed from the coefficients. And 2.0 did not store the absorbed group effects at all, so no fixed-effects prediction was possible.

`lm_robust()` now stores the absorbed group effects for a single-outcome model with one set of fixed effects, and `predict()` adds them back, reproducing the dummy-variable regression exactly. New levels in `newdata` and multi-way fixed effects both error with an explanation: with several sets of fixed effects the absorbed effects are identified in sum but not separately, so a new observation cannot be assigned its share.

### `glance()` on `iv_robust` with several endogenous regressors (#389)

The first-stage F test carries one entry per endogenous regressor, named `"<var>:value"` once there is more than one, so indexing it by `"value"` produced an NA with an NA name and the `data.frame()` call failed on the row name. `glance()` must return one row, so it now reports the weakest of the per-regressor first stages, which is the quantity a weak-instrument diagnostic asks about. The per-regressor tests remain in `diagnostic_first_stage_fstatistic`.

### `model.frame()` on `iv_robust` (#397)

`model.frame()` returned a two-column frame whose second column was named `"x + hp | am + hp"`, because `model.frame.default()` reads the `|` of a two-part formula as a logical operator. A `model.frame.iv_robust()` method now rebuilds a one-part formula naming every variable and defers to the default method.

### Smaller fixes

- **`lm_robust_fit()` with an unnamed or integer matrix (#269).** The exported fitting function failed on an integer design matrix ("Wrong R type for mapped matrix") and produced an object `tidy()` could not handle when `y` had no column names. Both inputs are now coerced and named.
- **`variable.names()` (#123).** Added for `lm_robust` and `iv_robust`.
- **`fixed_effects` must be a formula (#304).** Passing a bare column name produced "invalid formula" from deep inside `model.frame`. The error now names the argument and shows the expected form.
- **`lh_robust()` with multiple outcomes (#297).** Coefficients of a multivariate fit are named `"<outcome>:<term>"`, which a hypothesis such as `"cyl = 2"` cannot refer to, and `car::linearHypothesis()` reported it as malformed. 2.0 errors with that explanation and tells the user to fit one outcome at a time.

### All-missing outcome with weights (#370)

`lm_robust(X ~ 1, weights = w, data = df)` where `X` is entirely `NA` errors in estimatr 1.x with "'x' must be an array of at least two dimensions", while the same model without weights returns an NA coefficient. The asymmetry matters when one model is fitted across many subgroups and some subgroup has no observed outcome. 2.0 returns the NA coefficient in both cases.

---

## Performance improvements

### CR2 degrees-of-freedom computation: O(J²) → O(r²·J)

CR2 standard errors require, for each cluster j, a scalar adjustment derived from the hat matrix block H_jj. The original 1.x implementation formed a J×J matrix `P_array` and computed a trace and Frobenius norm over it — O(J²) memory and compute when J is large.

2.0 uses an algebraic identity to reduce this to five r×r Gram matrix products, where r is the number of model parameters (typically ≪ J). For a model with 3 parameters and 500 clusters the allocation drops from a 500×500 matrix to five 3×3 matrices. Measured speedup on a 1,000-observation / 100-cluster CR2 model: approximately 2× end-to-end.

### Formula parsing bypass for non-IV calls

`iv_robust()` requires `Formula::as.Formula()` to parse the two-sided formula `y ~ x | z`. `lm_robust()` and `lm_lin()` do not, but the original code called `as.Formula()` unconditionally in `clean_model_data()`. This path is now skipped for non-IV calls, saving approximately 80 µs per `lm_robust()` invocation — meaningful in simulations that call it thousands of times.

### Horvitz-Thompson variance: O(N²) → O(1)

1.x constructs the full 2N×2N joint inclusion probability matrix for every Horvitz-Thompson call. The Young's inequality variance for any complete-randomization design (including within blocks, and at the cluster level for clustered designs) reduces to six scalar sufficient statistics per block, so the matrix is never needed. The custom-design pathway, where an explicit permutation matrix is the only description of the design, is the one case that still forms a matrix, and does so in one `tcrossprod()` call in place of 1.x's 538-line `helper_condition_pr_matrix.R`.

Measured against 1.x on complete randomization: 10× at N = 200, 73× at N = 1,000, 837× at N = 3,000. Standard errors are identical to floating-point precision for simple, complete, blocked, clustered, and blocked-and-clustered designs.

### Multi-arm Horvitz-Thompson

estimatr 1.x refuses an `ra_declaration` with more than two arms. 2.0 supports the contrast between any two arms of a multi-arm design: `condition1` and `condition2` pick the contrast, the estimand stays the average treatment effect over all N units, and the variance uses the joint assignment probabilities implied by the arm sizes.

Two generalizations were needed. The denominator of the estimator is N, the size of the design, and not the number of units landing in the two selected conditions; those coincide only when the design has two arms. The Young's inequality coefficients come from the pairwise joint probabilities `m2(m2-1)/(n(n-1))`, `m1(m1-1)/(n(n-1))`, and `m2 m1/(n(n-1))`, which do not assume `m1 = n - m2`. The cross coefficient collapses to `1/n` for every complete design. Two-arm contrasts continue to run through 1.x's `gen_joint_pr_complete` path, including its floor/ceiling mixture for non-integer implied counts, so every two-arm result is numerically unchanged.

Tested by exhaustive enumeration: over all 90 assignments of six units to three arms of two, the average estimate equals the true ATE exactly, and the closed-form variance equals the one computed from the exact joint probability matrix of the same 90 permutations.

Because probabilities are matched to units by row position, `data` must contain one row per unit of the declaration. Passing a declaration whose size does not match `nrow(data)` is now an error rather than a silent misalignment.

### TSS computation vectorized

The total sum of squares calculation replaced an `apply(y, 1, '-', colMeans(y))` column-mean sweep with `sweep(y, 2, colMeans(y))`, which is handled by a single optimized C routine rather than an R-level loop. Speedup on a 1,000-observation multivariate model: approximately 280 µs.

---

## Fixed effects

Fixed effects absorption is implemented using the Frisch-Waugh-Lovell theorem with alternating projections (the algorithm used by `fixest::feols`). All FE variables are demeaned out of both the outcome and the design matrix before estimation; the intercept is absorbed into the group means and dropped.

**SE types with fixed effects.** With a **single** FE factor there is no restriction and no cost: the hat value of the full `[dummies | X]` design splits as `h_ii = h_ii(demeaned X) + w_i / sum(w over i's group)`, so HC2 and HC3 are exact from the demeaned fit plus one vector. The default is therefore `"HC2"`, the same as with no FE, and the same as 1.x returns.

```r
lm_robust(y ~ z, data = dat, fixed_effects = ~block)
# se_type = "HC2", bit-identical to 1.x, 17.8 ms where 1.0.6 takes 41.3 s
# at n = 40,000 with 2,000 blocks

lm_robust(y ~ z, data = dat, fixed_effects = ~block, clusters = cl)
# se_type = "CR0"

lm_robust(y ~ z, data = dat, fixed_effects = ~ block + year, se_type = "HC2")
# Error: `se_type = "HC2"` requires hat values from the full [X | FE dummies]
# design matrix and cannot be used with `fixed_effects`.
# To get HC2 SEs, replace `fixed_effects` with explicit dummies:
#   lm_robust(y ~ x + factor(fe_var), se_type = "HC2")
```

With **two or more** FE factors the identity fails, by roughly 3e-3 rather than by rounding, so HC2 and HC3 still error there and name the escape hatch, and the default falls back to `"stata"` (= HC1). `"CR2"` is refused with any `fixed_effects`: its adjustment is built from cluster-level blocks of the hat matrix rather than from `h_ii`. The dummy-variable formulation returns exactly the 1.x `fixed_effects` result in every case; that equivalence is tested.

**Multi-way FE.** Two-way (and higher) FE are supported via alternating projections that iterate to convergence (tolerance 1e-8, maximum 50 iterations).

**Return object additions.** FE models add:
- `fes`: logical flag indicating whether fixed effects were used
- `proj_r.squared`: within-group R² from the demeaned regression
- `r.squared`: full-model R² using the original (un-demeaned) outcome
- `fitted.values`: reconstructed to the original scale via FWL (projected residuals equal full residuals)

**Degrees of freedom.** For a model with one FE variable having B levels, `fe_rank = B` degrees of freedom are consumed (B − 1 dummy contrasts plus the absorbed intercept), so `df.residual = n − k − B` where k is the number of non-FE regressors.
