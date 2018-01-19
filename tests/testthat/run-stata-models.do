// This file fits many models in stata and outputs the estimates for comparison with estimatr

cd ~/Dropbox/coding/dd/estimatr/tests/testthat/

clear all
import delimited mtcars.csv

gen w = drat / 5

file open outf using stata-ests.txt, write r

reg mpg hp
mat V=e(V)
file write outf _n "classical" _tab (V[1,1]) _tab (V[2,2]) _tab (e(df_r))

reg mpg hp, vce(robust)
mat V=e(V)
file write outf _n "HC1" _tab (V[1,1]) _tab (V[2,2]) _tab (e(df_r))

reg mpg hp, vce(hc2)
mat V=e(V)
file write outf _n "HC2" _tab (V[1,1]) _tab (V[2,2]) _tab (e(df_r))

reg mpg hp, vce(hc3)
mat V=e(V)
file write outf _n "HC3" _tab (V[1,1]) _tab (V[2,2]) _tab (e(df_r))

reg mpg hp, vce(cluster cyl)
mat V=e(V)
file write outf _n "stata_cl" _tab (V[1,1]) _tab (V[2,2]) _tab (e(df_r))

reg mpg hp [aweight=w]
mat V=e(V)
file write outf _n "classicalw" _tab (V[1,1]) _tab (V[2,2]) _tab (e(df_r))

reg mpg hp [aweight=w], vce(robust)
mat V=e(V)
file write outf _n "HC1w" _tab (V[1,1]) _tab (V[2,2]) _tab (e(df_r))

reg mpg hp [aweight=w], vce(hc2)
mat V=e(V)
file write outf _n "HC2w" _tab (V[1,1]) _tab (V[2,2]) _tab (e(df_r))

reg mpg hp [aweight=w], vce(hc3)
mat V=e(V)
file write outf _n "HC3w" _tab (V[1,1]) _tab (V[2,2]) _tab (e(df_r))

reg mpg hp [aweight=w], vce(cluster cyl)
mat V=e(V)
file write outf _n "stata_clw" _tab (V[1,1]) _tab (V[2,2]) _tab (e(df_r))

file close outf
