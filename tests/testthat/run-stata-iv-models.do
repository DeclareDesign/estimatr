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
