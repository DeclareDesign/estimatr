// ===========================================================================
// estimatr: which leverage does Stata use for HC2/HC3 with 2SLS?
//
// WHAT THIS IS FOR
//
// estimatr's iv_robust() and R's sandwich package disagree about HC2 and HC3
// standard errors for two-stage least squares. They agree exactly on HC0, HC1
// and cluster-robust; the disagreement is only about which "leverage" value
// belongs in the HC2/HC3 correction, and it is worth about 1e-4 relative.
//
// Writing Xhat for the instrumented regressors, the two candidates are
//
//   (A) h_i = xhat_i' (Xhat'Xhat)^-1 xhat_i    "second-stage leverage" <- estimatr
//   (B) h_i = xhat_i' (Xhat'Xhat)^-1 x_i       "influence"             <- sandwich
//
// Neither is obviously right: for 2SLS the hat matrix is idempotent but NOT
// symmetric, so the usual HC2 unbiasedness argument does not pick one out.
// We are looking for an outside opinion.
//
// We have since largely answered this ourselves: (A) is the diagonal of an
// orthogonal projection and so is bounded in [0,1], while (B) is the diagonal
// of an idempotent but non-symmetric matrix and is not bounded at all. On a
// weak first stage (B) routinely exceeds 1, at which point sandwich's HC2 is
// undefined and returns NaN, while (A) is always computable. Where both are
// defined the two agree closely. So we are no longer asking which is correct;
// we would like to know, for the record, whether Stata implements this at all
// and if so which convention it picked.
//
// As far as we can tell, official `ivregress` does not accept vce(hc2) (its
// manual lists only unadjusted, robust, cluster, bootstrap, jackknife, hac),
// and `ivreg2` does not either. This file finds out for certain on YOUR Stata,
// and in either case collects the one thing Stata CAN produce that bears on
// the question: the delete-one jackknife, which HC3 is meant to approximate.
//
// The second dataset, weak-iv.csv, is built so the two conventions cannot be
// confused: there (A) peaks at 0.18 and (B) peaks at 49.9. If your Stata does
// accept vce(hc2) for ivregress, whether it returns a number or a missing
// value on that file identifies the convention on its own.
//
// WHAT TO DO
//
//   1. Put this file, mtcars.csv and weak-iv.csv in the same folder.
//   2. Open Stata, cd to that folder, and run:   do run-stata-iv-hc2.do
//   3. Send back the file  stata-iv-hc2.log  that it creates.
//
// It takes a few seconds, changes nothing on your system, installs nothing,
// and needs no packages. Every line we need is written into the log, so the
// log on its own is enough. Nothing outside this folder is read or written.
//
// If something fails, send the log anyway: a failure is a result here, and the
// file is written so that a refused option is recorded rather than fatal.
// ===========================================================================

clear all
set more off
capture log close _all
log using "stata-iv-hc2.log", text replace

display "RESULT|meta.stata_version|`c(stata_version)'"
display "RESULT|meta.version|`c(version)'"
display "RESULT|meta.flavor|`c(flavor)'"
display "RESULT|meta.os|`c(os)'"
display "RESULT|meta.date|`c(current_date)'"

import delimited "mtcars.csv", clear case(preserve)

// Fingerprint, so we can be certain your data is the data we have. The
// signature is Stata's own; the sums are a second check that does not depend
// on Stata's hashing.
datasignature
local sig "`r(datasignature)'"
display "RESULT|data.signature|`sig'"
display "RESULT|data.N|`=_N'"
foreach v of varlist mpg hp am wt gear cyl drat carb {
    quietly summarize `v'
    display "RESULT|data.sum.`v'|" %24.17e r(sum)
}

// A weight variable, defined exactly as in the 2019 do-files.
generate double w = drat / 5

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

capture program drop dumpfit
program define dumpfit
    args tag
    tempname b_ V_
    matrix `b_' = e(b)
    matrix `V_' = e(V)
    local k = colsof(`b_')
    local nm : colnames `b_'
    display "RESULT|`tag'.k|`k'"
    display "RESULT|`tag'.names|`nm'"
    display "RESULT|`tag'.N|" e(N)
    display "RESULT|`tag'.df_r|" e(df_r)
    forvalues i = 1/`k' {
        display "RESULT|`tag'.b_`i'|" %24.17e `b_'[1,`i']
        forvalues j = 1/`k' {
            display "RESULT|`tag'.V_`i'_`j'|" %24.17e `V_'[`i',`j']
        }
    }
end

// Run a command; record whether Stata accepted it, and dump the fit if so.
capture program drop tryfit
program define tryfit
    args tag cmdline
    capture noisily `cmdline'
    local rc = _rc
    display "RESULT|`tag'.rc|`rc'"
    if `rc' == 0 {
        dumpfit "`tag'"
    }
    else {
        display "RESULT|`tag'.refused|1"
    }
end

// ---------------------------------------------------------------------------
// 1. Does this Stata accept vce(hc2) / vce(hc3) for ivregress?
//
//    This is the question. rc == 0 means yes and the matrix follows; any other
//    return code means the option was refused, which is itself the answer.
//    Both the just-identified and the over-identified case are tried, because
//    the two leverage conventions separate further when the model is
//    over-identified.
// ---------------------------------------------------------------------------

tryfit iv_just_hc2   "ivregress 2sls mpg (hp am = wt gear), small vce(hc2)"
tryfit iv_just_hc3   "ivregress 2sls mpg (hp am = wt gear), small vce(hc3)"
tryfit iv_over_hc2   "ivregress 2sls mpg (hp = wt gear drat), small vce(hc2)"
tryfit iv_over_hc3   "ivregress 2sls mpg (hp = wt gear drat), small vce(hc3)"
tryfit iv_just_hc2_w "ivregress 2sls mpg (hp am = wt gear) [aweight = w], small vce(hc2)"

// ---------------------------------------------------------------------------
// 2. What Stata will certainly give us.
//
//    vce(jackknife) is the useful one: HC3 exists as an approximation to the
//    delete-one jackknife, so whichever convention tracks Stata's jackknife
//    more closely has the better claim. The others are anchors, and two of
//    them reproduce rows we already hold from 2019, which lets us confirm your
//    Stata behaves like the Stata those came from before we trust anything new.
// ---------------------------------------------------------------------------

tryfit iv_just_unadj "ivregress 2sls mpg (hp am = wt gear), small"
tryfit iv_just_rob   "ivregress 2sls mpg (hp am = wt gear), small rob"
tryfit iv_just_cl    "ivregress 2sls mpg (hp am = wt gear), small vce(cluster cyl)"
tryfit iv_just_jack  "ivregress 2sls mpg (hp am = wt gear), small vce(jackknife)"
tryfit iv_over_unadj "ivregress 2sls mpg (hp = wt gear drat), small"
tryfit iv_over_rob   "ivregress 2sls mpg (hp = wt gear drat), small rob"
tryfit iv_over_jack  "ivregress 2sls mpg (hp = wt gear drat), small vce(jackknife)"

// OLS hc2/hc3 for contrast. These we already match exactly, so if any of them
// disagrees with what we hold, the problem is the data or the Stata version
// rather than anything to do with 2SLS.
tryfit ols_hc2 "regress mpg hp, vce(hc2)"
tryfit ols_hc3 "regress mpg hp, vce(hc3)"

// ---------------------------------------------------------------------------
// 3. ivreg2, only if it happens to be installed. Not required.
// ---------------------------------------------------------------------------

capture which ivreg2
local has_ivreg2 = (_rc == 0)
display "RESULT|ivreg2.installed|`has_ivreg2'"
if `has_ivreg2' {
    tryfit ivreg2_rob "ivreg2 mpg (hp am = wt gear), robust small"
}

// ---------------------------------------------------------------------------
// 4. The two candidate leverage vectors, computed in Mata.
//
//    This is not Stata's opinion, it is our arithmetic done in Stata's linear
//    algebra rather than R's. It matters for two reasons: if vce(hc2) above
//    was accepted, comparing it against these two says which convention Stata
//    chose; and if it was refused, these still confirm that the two candidates
//    are what we think they are, computed by an independent implementation.
//
//    Column order is (hp, am, _cons), matching how ivregress reports.
// ---------------------------------------------------------------------------

mata:
void dumpiv(string scalar yv, string rowvector xv, string rowvector zv,
            string scalar tag)
{
    real colvector  y, u, ha, hb, hx, b, sel, yi, e
    real matrix     X, Z, PZX, Bread, H, IH, V, bj, D, Vj, Xi, Zi, Pi, dens
    real scalar     n, k, i, j, t
    real rowvector  bbar
    string rowvector tags

    y = st_data(., yv)
    n = rows(y)
    X = st_data(., xv), J(n, 1, 1)
    Z = st_data(., zv), J(n, 1, 1)
    k = cols(X)

    PZX   = Z * (invsym(quadcross(Z, Z)) * quadcross(Z, X))
    Bread = invsym(quadcross(PZX, PZX))
    b     = Bread * quadcross(PZX, y)
    u     = y - X * b

    ha = rowsum((PZX * Bread) :* PZX)   // second-stage leverage (estimatr)
    hb = rowsum((PZX * Bread) :* X)     // influence             (sandwich)

    H  = X * Bread * PZX'
    IH = I(n) - H
    hx = diagonal(IH * IH')             // exact E[e^2]/sigma^2

    printf("RESULT|%s.mata.n|%g\n", tag, n)
    printf("RESULT|%s.mata.k|%g\n", tag, k)
    printf("RESULT|%s.mata.sum_ha|%24.17e\n", tag, sum(ha))
    printf("RESULT|%s.mata.sum_hb|%24.17e\n", tag, sum(hb))
    printf("RESULT|%s.mata.max_ha|%24.17e\n", tag, max(ha))
    printf("RESULT|%s.mata.max_hb|%24.17e\n", tag, max(hb))
    printf("RESULT|%s.mata.max_asym|%24.17e\n", tag, max(abs(H - H')))
    printf("RESULT|%s.mata.max_idem|%24.17e\n", tag, max(abs(H * H - H)))

    for (i = 1; i <= n; i++) {
        printf("RESULT|%s.mata.ha_%g|%24.17e\n", tag, i, ha[i])
        printf("RESULT|%s.mata.hb_%g|%24.17e\n", tag, i, hb[i])
        printf("RESULT|%s.mata.hx_%g|%24.17e\n", tag, i, hx[i])
        printf("RESULT|%s.mata.u_%g|%24.17e\n",  tag, i, u[i])
    }
    for (i = 1; i <= k; i++)
        printf("RESULT|%s.mata.b_%g|%24.17e\n", tag, i, b[i])

    // sqrt of a negative entry gives a missing value rather than an error,
    // which is exactly the behaviour we want recorded for convention (B).
    dens = (sqrt(1 :- ha), sqrt(1 :- hb), sqrt(hx),
            (1 :- ha), (1 :- hb), hx, J(n, 1, 1))
    tags = ("hc2_a", "hc2_b", "hc2_x", "hc3_a", "hc3_b", "hc3_x", "hc0")
    for (t = 1; t <= 7; t++) {
        e = u :/ dens[., t]
        V = Bread * quadcross(PZX :* e, PZX :* e) * Bread
        for (i = 1; i <= k; i++)
            for (j = 1; j <= k; j++)
                printf("RESULT|%s.mata.%s.V_%g_%g|%24.17e\n",
                       tag, tags[t], i, j, V[i, j])
    }

    // Exact delete-one jackknife, to cross-check Stata's vce(jackknife).
    bj = J(n, k, .)
    for (i = 1; i <= n; i++) {
        sel = (1::n) :!= i
        yi = select(y, sel); Xi = select(X, sel); Zi = select(Z, sel)
        Pi = Zi * (invsym(quadcross(Zi, Zi)) * quadcross(Zi, Xi))
        bj[i, .] = (invsym(quadcross(Pi, Pi)) * quadcross(Pi, yi))'
    }
    bbar = mean(bj)
    D    = bj :- bbar
    Vj   = ((n - 1) / n) * quadcross(D, D)
    for (i = 1; i <= k; i++)
        for (j = 1; j <= k; j++)
            printf("RESULT|%s.mata.jack.V_%g_%g|%24.17e\n", tag, i, j, Vj[i, j])
}
end

// Column order is (hp, am, _cons), matching how ivregress reports.
mata: dumpiv("mpg", ("hp", "am"), ("wt", "gear"), "mtcars")

// ---------------------------------------------------------------------------
// 5. The weak-instrument file, where the two conventions cannot be confused.
// ---------------------------------------------------------------------------

import delimited "weak-iv.csv", clear case(preserve)
datasignature
local sig2 "`r(datasignature)'"
display "RESULT|weak.data.signature|`sig2'"
display "RESULT|weak.data.N|`=_N'"
foreach v of varlist y en x inst {
    quietly summarize `v'
    display "RESULT|weak.data.sum.`v'|" %24.17e r(sum)
}

tryfit weak_hc2   "ivregress 2sls y (en = inst) x, small vce(hc2)"
tryfit weak_hc3   "ivregress 2sls y (en = inst) x, small vce(hc3)"
tryfit weak_unadj "ivregress 2sls y (en = inst) x, small"
tryfit weak_rob   "ivregress 2sls y (en = inst) x, small rob"
tryfit weak_jack  "ivregress 2sls y (en = inst) x, small vce(jackknife)"
quietly regress en inst x
display "RESULT|weak.first_stage_r2|" %24.17e e(r2)

mata: dumpiv("y", ("en", "x"), ("inst", "x"), "weak")

display "RESULT|meta.done|1"

log close
