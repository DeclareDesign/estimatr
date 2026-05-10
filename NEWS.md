# estimatrZero 0.1.0

estimatrZero is a ground-up rewrite of [estimatr](https://declaredesign.org/r/estimatr/) targeting the DeclareDesign use case: OLS, Lin-adjusted OLS, 2SLS IV, difference-in-means, and linear hypothesis tests with heteroskedasticity- and cluster-robust standard errors. The rewrite fixes several long-standing correctness bugs, improves performance on the critical path, and adds feols-style fixed effects absorption.

---

## What is kept

All public functions from estimatr that are relevant to the DeclareDesign workflow are present with identical interfaces:

- `lm_robust()` — OLS with HC0/HC1/HC2/HC3/CR0/CR2/stata SEs
- `lm_lin()` — Lin (2013) covariate-adjusted estimator
- `iv_robust()` — two-stage least squares with the same SE menu
- `difference_in_means()` — Neyman variance, paired, blocked, and clustered designs
- `lh_robust()` — linear hypothesis tests via `car::linearHypothesis` with robust variance
- S3 methods: `tidy`, `glance`, `summary`, `print`, `predict`, `coef`, `confint`, `vcov`, `nobs`, `update`

Return objects are structurally compatible with estimatr 1.0.6: field names, classes, and method dispatch are unchanged. Cross-package identity tests confirm coefficient and standard error agreement to 1e-12 across all supported SE types.

---

## What is dropped

**Horvitz-Thompson estimator.** `horvitz_thompson()` is not included. The HT estimator requires a known randomization probability matrix, which in the DeclareDesign workflow would have to be recomputed for every simulation draw. Supplying it adds complexity and latency that outweigh the benefits for that use case; `difference_in_means()` covers the common cases.

**`starprep()` and `commarobust()`.** Stargazer integration utilities not relevant to the DeclareDesign workflow.

**HAC (Newey-West) standard errors.** Not implemented.

**CR3 standard errors.** Not implemented.

---

## Bug fixes

### Weighted R² formula (affects `lm_robust`, `lm_lin`, `iv_robust`)

estimatr's weighted R² formula is correct. The standard weighted TSS is:

```
TSS = Σᵢ wᵢ (yᵢ - ȳ_w)²
```

where ȳ_w is the standard weighted mean (weights wᵢ). Inside `lm_robust_fit`, raw weights are rescaled to `weights_internal = sqrt(w / μ_w)`, so `weights_internal² = w / μ_w`. Using `weights_internal²` in the TSS formula recovers the correct standard weighted mean as the centering point. Using `weights_internal` instead (with sqrt-weights) centers on ȳ_{sqrt(w)} ≠ ȳ_w, which breaks the WLS guarantee that RSS ≤ TSS and allows R² < 0.

estimatrZero uses `weights_internal²` throughout, matching estimatr's formula exactly. For an intercept-only WLS model, R² = 0 exactly. For any WLS model, R² ∈ [0, 1].

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

---

## Performance improvements

### CR2 degrees-of-freedom computation: O(J²) → O(r²·J)

CR2 standard errors require, for each cluster j, a scalar adjustment derived from the hat matrix block H_jj. The original estimatr implementation formed a J×J matrix `P_array` and computed a trace and Frobenius norm over it — O(J²) memory and compute when J is large.

estimatrZero uses an algebraic identity to reduce this to five r×r Gram matrix products, where r is the number of model parameters (typically ≪ J). For a model with 3 parameters and 500 clusters the allocation drops from a 500×500 matrix to five 3×3 matrices. Measured speedup on a 1,000-observation / 100-cluster CR2 model: approximately 2× end-to-end.

### Formula parsing bypass for non-IV calls

`iv_robust()` requires `Formula::as.Formula()` to parse the two-sided formula `y ~ x | z`. `lm_robust()` and `lm_lin()` do not, but the original code called `as.Formula()` unconditionally in `clean_model_data()`. This path is now skipped for non-IV calls, saving approximately 80 µs per `lm_robust()` invocation — meaningful in simulations that call it thousands of times.

### TSS computation vectorized

The total sum of squares calculation replaced an `apply(y, 1, '-', colMeans(y))` column-mean sweep with `sweep(y, 2, colMeans(y))`, which is handled by a single optimized C routine rather than an R-level loop. Speedup on a 1,000-observation multivariate model: approximately 280 µs.

---

## Fixed effects

Fixed effects absorption is implemented using the Frisch-Waugh-Lovell theorem with alternating projections (the algorithm used by `fixest::feols`). All FE variables are demeaned out of both the outcome and the design matrix before estimation; the intercept is absorbed into the group means and dropped.

**SE type defaults with fixed effects.** HC2 and HC3 require the full n×n hat matrix over the augmented design matrix including all FE dummies, which eliminates the computational advantage of FE demeaning. estimatrZero instead defaults to HC1 (no clusters) or `stata` CR0-scaled (with clusters) when `fixed_effects` is specified, and warns if HC2, HC3, or CR2 are requested:

```r
lm_robust(y ~ z, data = dat, fixed_effects = ~block)
# se_type = "HC1"

lm_robust(y ~ z, data = dat, fixed_effects = ~block, clusters = cl)
# se_type = "stata"

lm_robust(y ~ z, data = dat, fixed_effects = ~block, se_type = "HC2")
# Warning: 'HC2' is not supported with fixed_effects in estimatrZero...
# se_type = "HC1"
```

**Multi-way FE.** Two-way (and higher) FE are supported via alternating projections that iterate to convergence (tolerance 1e-8, maximum 50 iterations).

**Return object additions.** FE models add:
- `fes`: logical flag indicating whether fixed effects were used
- `proj_r.squared`: within-group R² from the demeaned regression
- `r.squared`: full-model R² using the original (un-demeaned) outcome
- `fitted.values`: reconstructed to the original scale via FWL (projected residuals equal full residuals)

**Degrees of freedom.** For a model with one FE variable having B levels, `fe_rank = B` degrees of freedom are consumed (B − 1 dummy contrasts plus the absorbed intercept), so `df.residual = n − k − B` where k is the number of non-FE regressors.
