// This file fits many models in stata and outputs the estimates for comparison with estimatr

clear all
import delimited mtcars.csv

gen w = drat / 5

reg mpg hp
reg mpg hp, vce(robust)
reg mpg hp, vce(hc2)
reg mpg hp, vce(hc3)
reg mpg hp, vce(cluster cyl)

reg mpg hp [aweight=w]
reg mpg hp [aweight=w], vce(robust)
reg mpg hp [aweight=w], vce(hc2)
reg mpg hp [aweight=w], vce(hc3)
reg mpg hp [aweight=w], vce(cluster cyl)
