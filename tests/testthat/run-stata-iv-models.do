// This file fits many models in stata and outputs the estimates for comparison with estimatr

clear all
import delimited mtcars.csv

gen w = drat / 5

file open outf using stata-iv-ests.txt, write r

ivregress 2sls mpg (hp am = wt gear), small
mat V=e(V)
file write outf _n "classical" _tab (V[1,1])  _tab (V[2,2])  _tab (V[3,3]) _tab (e(F)) _tab (e(r2)) _tab (e(r2_a)) _tab (e(rmse))

ivregress 2sls mpg (hp am = wt gear), small rob
mat V=e(V)
file write outf _n "rob" _tab (V[1,1])  _tab (V[2,2])  _tab (V[3,3]) _tab (e(F)) _tab (e(r2)) _tab (e(r2_a)) _tab (e(rmse))

ivregress 2sls mpg (hp am = wt gear), small vce(cluster cyl)
mat V=e(V)
file write outf _n "cl" _tab (V[1,1])  _tab (V[2,2])  _tab (V[3,3]) _tab (e(F)) _tab (e(r2)) _tab (e(r2_a)) _tab (e(rmse))

ivregress 2sls mpg (hp am = wt gear) [aweight = w], small
mat V=e(V)
file write outf _n "classical_w" _tab (V[1,1])  _tab (V[2,2])  _tab (V[3,3]) _tab (e(F)) _tab (e(r2)) _tab (e(r2_a)) _tab (e(rmse))

ivregress 2sls mpg (hp am = wt gear) [aweight = w], small rob
mat V=e(V)
file write outf _n "rob_w" _tab (V[1,1])  _tab (V[2,2])  _tab (V[3,3]) _tab (e(F)) _tab (e(r2)) _tab (e(r2_a)) _tab (e(rmse))

ivregress 2sls mpg (hp am = wt gear) [aweight = w], small vce(cluster cyl)
mat V=e(V)
file write outf _n "cl_w" _tab (V[1,1])  _tab (V[2,2])  _tab (V[3,3]) _tab (e(F)) _tab (e(r2)) _tab (e(r2_a)) _tab (e(rmse))

file close outf


cap file close outfdiag
file open outfdiag using stata-iv-diagnostics.txt, write r

#delimit ;
local formulae = `" 
"(hp = wt)"
"(hp am = wt gear)"
"gear (hp = wt)"
"gear (hp = wt am)"
"' ;
local options = `" 
"small" 
"rob" 
"cluster(cyl)" 
"small noconstant" 
"rob noconstant" 
"cluster(cyl) noconstant" 
"' ;
#delimit cr

foreach f in `formulae' {
	display "`f'"
	foreach opt in `options' {
		ivregress 2sls mpg `f', `opt'
		estat firststage, all
		mat singleresults=r(singleresults)
		local rows = rowsof(singleresults)
		forvalues i=1/`rows' {
			cap file write outfdiag "`f'" ";" "`opt'" ";" "weak`i'" ";" (singleresults[`i',5]) ";" (singleresults[`i',6]) ";" (singleresults[`i',4]) ";" (singleresults[`i',7]) _n 
		}
		estat endogenous
		if strpos("`opt'", "rob") > 0 | strpos("`opt'", "cluster") > 0 {
			file write outfdiag "`f'" ";" "`opt'" ";" "endog" ";" (r(regFdf_n)) ";" (r(regFdf_d)) ";" (r(regF)) ";" (r(p_regF)) _n 
			cap estat overid
			file write outfdiag "`f'" ";" "`opt'" ";" "overid" ";" (r(df)) ";.;" (r(score)) ";" (r(p_score)) _n 
		} 
		else {
			file write outfdiag "`f'" ";" "`opt'" ";" "endog" ";" (r(df)) ";" (r(wudf_r)) ";" (r(wu)) ";" (r(p_wu)) _n 
			cap estat overid
			file write outfdiag "`f'" ";" "`opt'" ";" "overid" ";" (r(df)) ";.;" (r(sargan)) ";" (r(p_sargan)) _n 
		}
	}

}

file close outfdiag


ivregress 2sls mpg gear (hp = wt am), rob
