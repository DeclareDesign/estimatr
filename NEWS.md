# estimatr 2.0.0

estimatr 2.0.0 is a ground-up rewrite aimed at the DeclareDesign use case: OLS, Lin-adjusted OLS, 2SLS IV, difference-in-means, Horvitz-Thompson, and linear hypothesis tests with heteroskedasticity- and cluster-robust standard errors. It fixes several long-standing correctness bugs, improves performance on the critical path, adds feols-style fixed effects absorption, and replaces the O(N²) Horvitz-Thompson variance with a design-aware O(1) computation. Interfaces are unchanged from 1.0.6 and, run side by side on one machine, the numbers agree to 1e-12 wherever both versions answer.

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

### A note on the comparison tolerance

The suite compares against a fixture of answers recorded from an installed
1.0.6 rather than calling that version live, which is what lets the package
carry its own name. One consequence was not obvious until CI ran: **the
recording is a fixed set of numbers made on one machine, so comparing it to
values computed on another machine has a floor set by BLAS and LAPACK rather
than by this package.** Before the freeze both sides were computed in the same
process, so several comparisons could reasonably ask for exact equality; they
cannot any more. The first CI run after the rename failed 11 assertions on
Ubuntu and Windows and none on macOS, where the fixture was recorded, and every
difference was too small for waldo to render.

Recorded comparisons therefore run at `REF_TOL = 1e-9`, set from the worst case
rather than by taste: the tightest recorded quantity is the weighted
`adj.r.squared` at -7.1e-4, whose one ulp is a relative difference of 3.1e-13,
leaving about 3,200 ulps of headroom, and nothing else in the fixture comes
within 100 ulps of the tolerance. Comparisons between two things computed in
the same session keep their tight tolerances, since those have no platform
floor.

`tests/testthat/test_return_surface.R` now pins the whole returned surface of
sixteen fit types against answers recorded from 1.0.6, so a field cannot go
missing again without a test failing. A missing field is the failure mode worth
guarding: `fit$felevels` returns `NULL` rather than erroring, so the reader gets
a wrong answer instead of a stop.

---

## Two capabilities that had been dropped, and are now restored

Both were removed on the understanding that 1.x had restricted them too. It had
not: 1.0.6 answered every one of these calls, and answered them correctly. The
reverse-dependency run is what surfaced the mistake.

**HC2, HC3 and CR2 work with `fixed_effects` again.** 1.0.6 built the full
dummy matrix and took the hat values off the `[X | FE dummies]` design whenever
the requested `se_type` needed them. 2.0 does the same. One fixed-effect factor
still takes the cheap route through the leverage identity below; everything
else expands the dummies. Every case is exact and equal to the explicit-dummy
fit, and equal to 1.0.6, which the fixture now records. The one combination
1.0.6 also refused, weighted CR2 with `fixed_effects`, is still refused.

**A bare grouping vector passed to `fixed_effects` warns instead of erroring.**
1.x accepted a bare column name or an already-evaluated vector. Issue #304
asked for the formula to be enforced; a warning that names the argument and
shows the expected form does that without breaking working code, so the vector
is still accepted and still names the term the way 1.0.6 named it. Anything
that is neither a formula nor a grouping vector is an error, as before.

### One default that moves, and says so

`fixed_effects` with `clusters` defaults to `"CR0"`; 1.0.6 defaulted to
`"CR2"`. CR2 is the one estimator here still built from cluster-level blocks of
the hat matrix rather than from `h_ii`, so it is the one case that still has to
expand the dummies, and a default that quietly costs O(g^3) in the number of
levels would undo the reason for absorbing fixed effects at all. The warning
fires **once per session**, since absorbed fixed effects are usually fitted in
a loop or a simulation and a per-call warning would be noise. Writing
`se_type = "CR0"` accepts the default and removes the warning; writing
`se_type = "CR2"` returns the 1.0.6 number exactly, also without warning.

**Unclustered `fixed_effects` keeps the `"HC2"` default at any number of
factors**, which is what 1.0.6 returned. An earlier draft of 2.0 fell back to
`"HC1"` for two or more factors and warned, on the understanding that HC2 there
meant expanding the dummies. It does not; see below.

### What the reverse-dependency run leaves broken

Of estimatr's 35 CRAN reverse dependencies, 34 check clean against 2.0.0.
(`hbal` at first would not install on the test machine at all, for a local
gfortran path problem unrelated to this package; corrected, it checks clean
too.)

`clubSandwich`, `RCT` and `statuser` were broken by the two removals above and
now check clean.

One remains:

* **`eventstudyr`** has three assertions that read `felevels$V1` on a one-way
  fixed-effects fit and want the 1.0.6 fallback name, in
  `test-EventStudyFHS.R` and `test-EventStudyOLS.R`. Nothing else in its suite
  changes and no estimate moves: its one-way branches size their
  degrees-of-freedom correction as `1 + rank` without consulting `felevels`,
  and its comparisons against recorded Stata output all hold.

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

Return objects are structurally compatible with estimatr 1.0.6: field names, classes, and method dispatch are unchanged. Cross-version identity tests confirm coefficient and standard error agreement across all supported SE types, to 1e-12 when both versions are run on the same machine and to 1e-9 against the recorded fixture the suite ships.

`horvitz_thompson()` fits now carry the same post-estimation methods as every other estimator here. `vcov()`, `confint()`, `glance()`, `summary()` and `nobs()` were missing, so the first three errored and the last two fell through to the base defaults and returned something that resembled output without being the estimate. All five now agree with estimatr 1.0.6 exactly, including `confint(level = )` at a non-default level, where the interval is normal rather than t because the estimator has no degrees of freedom to spend. The `glance()` gap was the one worth finding: without it `modelsummary()` did not fail, it silently dropped the goodness-of-fit rows and printed a coefficient table that looked complete.

---

## HC2 and HC3 with fixed effects

Absorbed fixed effects used to refuse `se_type = "HC2"` and `"HC3"`, on the grounds that they need hat values from the full `[dummies | X]` design while absorption leaves only the demeaned ones. The hat values in fact decompose exactly, for any number of FE factors:

```
P_[X | D] = P_D + P_{M_D X}
```

so `h_ii` is the demeaned-X hat value plus `diag(P_D)`. With one factor `P_D` is diagonal and the second term is just `w_i / sum(w over i's group)`. With several it is not diagonal, but writing `D` with its widest factor in full dummies and the rest contrast-coded leaves `D'WD` with a diagonal leading block, so `diag(P_D)` costs a factorisation of order `sum_{k>1}(g_k - 1)` -- the narrowest the design allows -- and no dummy matrix is built at any point.

Verified to machine precision against writing the dummies out, weighted and unweighted, balanced and unbalanced, from one factor to five, and against 1.x. At n = 40,000 with 2,000 groups, one-way HC2 takes 17.8 ms here against 41.3 s in 1.0.6. At n = 50,000 across 1,000 x 30 groups, two-way HC2 takes 36 ms and peaks at 174 MB, against 13.8 s and 970 MB before.

**The default with fixed effects therefore returns to `"HC2"`**, at any number of factors, which is the package default everywhere else and is what 1.x returns for the same call. It had fallen back to `"stata"` only because HC2 was unaffordable.

The identity generalises to any number of factors. What fails beyond one factor is the *sum of each factor's own within-group share*, which ignores the cross terms and is wrong in the third decimal place; that is what an earlier draft tested and rejected. The projection itself decomposes exactly:

```
P_[X | D] = P_D + P_{M_D X}
```

so `h_ii` is the demeaned-X hat value plus `diag(P_D)`, whatever the number of factors. Writing the FE design with its widest factor in full dummies and the rest contrast-coded makes `D'WD` block-diagonal in its leading block, so `diag(P_D)` costs a factorisation of order `sum_{k>1}(g_k - 1)` -- as small as the design allows -- rather than an n-by-g dummy hat matrix. Inverting that block through its eigendecomposition rather than a solve means a **disconnected or nested** design, where `D` is rank deficient, gets the pseudo-inverse and the right projection anyway.

`iv_robust()` gets the same treatment: its second stage runs on fitted regressors, but those are demeaned by the same fixed effects, so the decomposition applies unchanged. It never used the shortcut before, even for one factor.

`"CR2"` is the exception. Its adjustment is built from cluster-level *blocks* of the hat matrix rather than from `h_ii`, and blocks do not decompose the way the diagonal does, so CR2 with `fixed_effects` still expands the dummies exactly as 1.0.6 did. It is available and exact; it simply costs what 1.0.6 charged for it.

**What this is worth.** On 50,000 observations across 1,000 x 30 groups, `se_type = "HC2"` with two-way `fixed_effects` took 13.8 seconds and peaked at 970 MB; it now takes 36 milliseconds and peaks at 174 MB. The dummy matrix is never built, so the memory that used to make wide fixed effects impractical is simply not allocated.

**A rank-deficiency bug, fixed.** When one FE factor is spanned by the others -- a nested factor, or a disconnected design -- the FE design is rank deficient, and the nominal count `sum(levels) - K + 1` overstates its rank. 1.0.6 used that nominal count for the residual degrees of freedom, and for HC2 and HC3 it expanded the dummies, let a pivoted QR drop the redundant columns, and then read the hat values off the padded design. Its absorbed answer therefore disagreed with its own explicit-dummy fit, and disagreed differently depending on the `se_type`.

2.0 takes the exact rank from the same eigendecomposition that produces the leverage, and uses it for **every** `se_type`, so an absorbed fit and the same fit with the dummies written out now agree exactly on standard errors, degrees of freedom, p-values and intervals however degenerate the design is. Absorbing fixed effects is a computational choice and should not move the answer.

This is a deliberate departure from 1.0.6, and it only bites where the FE design is genuinely rank deficient; on a design of full rank the nominal and exact counts agree and nothing moves. Canvassed against the ecosystem on a nested-factor design and a disconnected two-way design: base `lm()` with explicit dummies and `plm` (twoways within) both report the exact rank, and estimatr now matches them. `fixest::feols` reports the nominal one and issues no note, so it differs from its own dummy regression on both designs; that difference is `fixest`'s to make, and is recorded here only so the discrepancy is not mistaken for a bug in either package.

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
- **`fixed_effects` given a bare column name (#304).** It produced "invalid formula" from deep inside `model.frame`. It now warns, names the argument, shows the expected form, and goes on to fit the model the way 1.x did.
- **`lh_robust()` with multiple outcomes (#297).** Coefficients of a multivariate fit are named `"<outcome>:<term>"`, which a hypothesis such as `"cyl = 2"` cannot refer to, and `car::linearHypothesis()` reported it as malformed. 2.0 errors with that explanation and tells the user to fit one outcome at a time.

### A single fixed-effects factor lost its name when rows were dropped

With one factor in `fixed_effects` and any observation dropped for missingness, the model frame's one-column matrix of factor codes degrades to a vector, and rebuilding it as a matrix leaves it without a column name. Three things read that name, and all three failed. `se_type = "CR2"` errored outright with "'.' in formula and no 'data' argument", because CR2 is the one estimator that still expands the fixed effects and it does so through `model.matrix(~ 0 + .)` on a data frame that had no names to expand. The absorbed group effects in `fixed_effects` came back labelled `1`, `2`, `3` instead of `g1`, `g2`, `g3`. And `predict()` rejected its own data as containing new levels, since the levels it was matching against had lost the term.

estimatr 1.0.6 lost the name in the same place and fell back to `V1`, so it returned `V11`, `V12`, `V13` and its `predict()` failed too, but it did answer the CR2 call. The name is now taken from the level names captured before the model frame is built, which are the same ones `felevels` reports, so the two cannot disagree and complete and incomplete data give the same labels. The restored CR2 fit agrees with 1.0.6 and with the explicit-dummy fit to twelve digits.

Every fixed-effects test until now used complete data, which is why nothing reported this.

### `res_var` and `weights` were named by accident on weighted fixed-effects fits

Both fields took whatever dimnames survived `data[["y"]] - fitted.values`, and R resolves that in favour of the first operand carrying any. The outcome was carrying the model frame's row names, except on the fixed-effects path, where it has been demeaned into a fresh matrix. So a weighted fixed-effects fit disagreed with a weighted fit without fixed effects, and with 1.0.6, in two fields at once: `res_var` (and the `proj_` R-squareds derived from it) came back named after the outcome column when 1.0.6 leaves it unnamed, and `weights` came back unnamed when 1.0.6 names it. Now that the row names are carried explicitly rather than inherited, all three cases agree with 1.0.6. Verified against an installed 1.0.6 and pinned in `test_return_surface.R`.

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

The total sum of squares calculation replaced an `apply(y, 1, '-', colMeans(y))` column-mean sweep with direct centring and `colSums()`. 1.x ran an R-level loop over rows; an intermediate version of 2.0 used `sweep()`, which builds its subtrahend with `array()` and then transposes it with `aperm()`, so it materialised two extra n-by-k copies to do what recycling does in one. Speedup on a 1,000-observation multivariate model: approximately 280 µs, and two fewer allocations the size of the outcome on every fit.

### Front-end overhead

At n = 100,000 the numerical core of a fit is a small fraction of the call. Most of the rest was the model frame's row names.

`stats::model.frame()` turns a data frame's compact automatic row names into a character vector as long as the data, and from there every object derived from the design matrix carries them: the fitted values from `X %*% beta`, the residuals formed from those, the row subset the cluster sort takes, and the `drop()` that turns each n-by-1 matrix into a vector. Only `fitted.values` keeps them in the returned object; estimatr returns residuals unnamed, as 1.0.6 did, so on that branch they were built only to be discarded.

They are now taken off the design matrix and the outcome once, in `clean_model_data()`, where both are still unshared and stripping them is free, and put back on `fitted.values` once in `lm_return()`. Nothing in between has to carry them. `weights` is named to match, as before.

The rest are local. Gaps in fixed-effect codes are detected with `tabulate()` rather than by hashing every observation with `sort(unique())`. Unweighted group sizes in the leverage computation are `tabulate()` rather than `rowsum()`, which rediscovers the groups it was already given, and the unit weight vector is no longer materialised on the one-way path that does not index it. The sums of squares behind the R-squared of a fixed-effects fit go through `crossprod()` instead of squaring the vector first. The absorbed group effects find a representative row per group by scatter assignment over the codes rather than by hashing all n of them. `sweep()` is gone from the total sum of squares, and the design matrix is no longer copied to itself when no column was dropped for collinearity.

### The fixed-effects kernels

Two things dominated a multi-way fixed-effects fit, and both were doing work they already had the answer to.

`fe_leverage()` cross-tabulates each pair of factors. It built a composite index and handed it to `rowsum()`, which hashes all n observations to rediscover groups that are known to run 1..g, and then labels its result rows with those group values as strings, so the caller converted them back with `as.integer(rownames(.))`. At n = 100,000 over 1,000 by 30 groups that pair cost 7.3 ms, against 0.04 ms for a single counting pass, and there is one per factor pair. A dense weighted cross-tabulation in C++ (`xtab_cpp`) replaces it and serves the weighted case too, so the unit weight vector is no longer built for an unweighted fit. The masks and column offsets the leverage vector needs are now built only when a leverage vector is wanted, which is HC2 and HC3; every other `se_type` returns the rank without them.

`demean_cpp()` rebuilt the group weight sums on every sweep of the alternating projections, though they depend only on the codes and the weights and are identical each time, and then made a second full pass over the matrix to find its largest absolute value for the convergence test. The sums are now computed once, the maximum is accumulated by the subtraction pass that already touches every cell, and an unweighted fit no longer multiplies each element by a 1.0 it had just read out of a vector of ones. Every result is bit-identical; there is simply less of it.

Measured at n = 100,000 with two regressors, against 2.0.0 before this work: a plain fit 15.8 to 3.2 ms, a clustered fit 18.5 to 5.1 ms, a one-way fixed-effects fit 14.7 to 5.8 ms, a two-way fixed-effects fit 32.7 to 9.3 ms. On the cases that carry the most fixed-effects machinery: two-way HC2 36.7 to 12.8 ms, two-way weighted 35.0 to 12.6 ms.

All of it is in `clean_model_data()`, `lm_robust_fit()`, `lm_return()` and the fixed-effects helpers, so `iv_robust()` and `lm_lin()` get it for nothing: 2SLS 20.0 to 5.5 ms, clustered 2SLS 30.0 to 7.6 ms, 2SLS with two-way fixed effects 37.8 to 14.4 ms, `lm_lin()` 26.0 to 6.1 ms.

Every returned value is unchanged. Across 56 specifications covering each `se_type`, weights, clusters, multivariate outcomes, IV, and one to two fixed-effect factors, the largest relative difference from the previous version is 9e-15, and the names on every field are identical apart from the three naming bugs fixed below.

---

## Fixed effects

Fixed effects absorption is implemented using the Frisch-Waugh-Lovell theorem with alternating projections (the algorithm used by `fixest::feols`). All FE variables are demeaned out of both the outcome and the design matrix before estimation; the intercept is absorbed into the group means and dropped.

**SE types with fixed effects.** With a **single** FE factor there is no restriction and no cost: the hat value of the full `[dummies | X]` design splits as `h_ii = h_ii(demeaned X) + w_i / sum(w over i's group)`, so HC2 and HC3 are exact from the demeaned fit plus one vector. The default is therefore `"HC2"`, the same as with no FE, and the same as 1.x returns.

```r
lm_robust(y ~ z, data = dat, fixed_effects = ~block)
# se_type = "HC2", bit-identical to 1.x, 17.8 ms where 1.0.6 takes 41.3 s
# at n = 40,000 with 2,000 blocks

lm_robust(y ~ z, data = dat, fixed_effects = ~block, clusters = cl)
# se_type = "CR0", with a warning: 1.0.6 defaulted to "CR2" here

lm_robust(y ~ z, data = dat, fixed_effects = ~ block + year)
# se_type = "HC2", as in 1.0.6, and no warning: the identity covers any
# number of factors, so there is nothing to trade away

lm_robust(y ~ z, data = dat, fixed_effects = ~ block + year, se_type = "CR2",
          clusters = cl)
# se_type = "CR2", exact, by expanding the dummies. No warning: you asked.
```

With **two or more** FE factors the naive version of the identity fails, by roughly 3e-3 rather than by rounding: summing each factor's own within-group share drops the cross terms. The projection decomposition `P_[X | D] = P_D + P_{M_D X}` is exact for any number of factors, and that is what is used. Only `"CR2"` with `fixed_effects` still expands the dummies, its adjustment being built from cluster-level blocks of the hat matrix rather than from `h_ii`. Every answer equals the explicit-dummy fit and equals 1.0.6; that equivalence is tested and recorded in the fixture.

What changes is only what a **default** will reach for, and only in the clustered case. Expanding the dummies is roughly O(g^3) in the number of levels, so no default expands them; `"CR2"` is the only estimator that still needs the expansion, so `fixed_effects` with `clusters` defaults to `"CR0"` with a warning, shown once per session, that names what 1.0.6 returned and names the `se_type` that accepts the new default. Writing `se_type = "CR0"` removes it. Unclustered fits keep the `"HC2"` default at any number of factors, as in 1.0.6.

**Multi-way FE.** Two-way (and higher) FE are supported via alternating projections that iterate to convergence (tolerance 1e-8, maximum 50 iterations).

**Return object additions.** FE models add:
- `fes`: logical flag indicating whether fixed effects were used
- `proj_r.squared`: within-group R² from the demeaned regression
- `r.squared`: full-model R² using the original (un-demeaned) outcome
- `fitted.values`: reconstructed to the original scale via FWL (projected residuals equal full residuals)

**Degrees of freedom.** For a model with one FE variable having B levels, `fe_rank = B` degrees of freedom are consumed (B − 1 dummy contrasts plus the absorbed intercept), so `df.residual = n − k − B` where k is the number of non-FE regressors.
