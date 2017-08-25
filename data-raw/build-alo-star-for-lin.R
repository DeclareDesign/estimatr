# This script takes the dta file from the ALO replication data at this url:
# https://www.aeaweb.org/articles?id=10.1257/app.1.1.136
# and turns it in to the .rda we export with the package
# Full citation:
# Angrist, Joshua, Daniel Lang, and Philip Oreopoulos. 2009. “Incentives and
# Services for College Achievement: Evidence from a Randomized Trial.” American
# Economic Journal: Applied Economics 1(1): 136–63.

# wd should be package root

dat <-
  foreign::read.dta(
    'data-raw/STAR_public_use.dta'
  )

alo_star_men <- dat[dat$sex == 'M' &
                      !is.na(dat$GPA_year1) &
                      (dat$sfsp == 1 | dat$ssp == 1) &
                      dat$noshow == 0,
                    c('gpa0', 'sfsp', 'ssp', 'sfp', 'GPA_year1', 'GPA_year2')
                    ]

devtools::use_data(alo_star_men,
                   internal = F,
                   overwrite = T)
