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



    mpg ~ hp | wt,
    mpg ~ 0 + hp | 0 + wt,
    mpg ~ hp + am | wt + gear,
    mpg ~ 0 + hp + am | 0 + wt + gear,
    mpg ~ hp + gear | wt + gear,
    mpg ~ 0 + hp + gear | 0 + wt + gear,
    mpg ~ hp + gear | wt + gear + am,
    mpg ~ 0 + hp + gear | 0 + wt + gear + am
	
	
cap file close outfdiag
file open outfdiag using stata-iv-diagnostics.txt, write r

#delimit ;
local formulae = `" 
"(hp = wt)"
"(hp am = wt gear)"
"gear (hp = wt)"
"gear (hp = wt am)"
"' ;
#delimit cr
foreach f in `formulae' {
	display "`f'"
	ivregress 2sls mpg `f', small
	estat firststage, all
	mat singleresults=r(singleresults)
	local rows = rowsof(singleresults)
	forvalues i=1/`rows' {
		cap file write outfdiag "`f'" _tab "weak`i'" _tab (singleresults[`i',5]) _tab (singleresults[`i',6]) _tab (singleresults[`i',4]) _tab (singleresults[`i',7]) _n 
	}
	estat endogenous
	file write outfdiag "`f'" _tab "wu-hausman" _tab (r(df)) _tab (r(wudf_r)) _tab (r(wu)) _tab (r(p_wu)) _n 
	cap estat overid
	file write outfdiag "`f'" _tab "overid" _tab (r(df)) _tab (r(sargan)) _tab (r(p_sargan)) _n 
}

file close outfdiag


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
