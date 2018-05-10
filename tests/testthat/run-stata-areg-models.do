// This file fits many models in stata and outputs the estimates for comparison with estimatr

clear all
import delimited mtcars.csv

gen w = drat / 5

file open outf using stata-fe-ests.txt, write r

// xtset carb
// xtreg mpg hp, fe

areg mpg hp, absorb(carb)
mat V=e(V)
file write outf _n "classical" _tab (V[1,1]) _tab (e(F))

areg mpg hp, absorb(carb) vce(robust)
mat V=e(V)
file write outf _n "HC1" _tab (V[1,1]) _tab (e(F))

areg mpg hp, absorb(carb) vce(cluster cyl)
mat V=e(V)
file write outf _n "stata_cl" _tab (V[1,1]) _tab (e(F))

areg mpg hp [aweight=w], absorb(carb)
predict hii, hat
mat V=e(V)
file write outf _n "classicalw" _tab (V[1,1]) _tab (e(F))

areg mpg hp [aweight=w], absorb(carb) vce(robust)
mat V=e(V)
file write outf _n "HC1w" _tab (V[1,1]) _tab (e(F))

areg mpg hp [aweight=w], absorb(carb) vce(cluster cyl)
mat V=e(V)
file write outf _n "stata_clw" _tab (V[1,1]) _tab (e(F))

file close outf
