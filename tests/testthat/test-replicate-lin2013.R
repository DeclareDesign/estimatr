context("lm lin replicates lin 2013")
# Lin paper available here: www.stat.berkeley.edu/~winston/agnostic.pdf
# Citation:
# Lin, Winston. 2013. "Agnostic notes on regression adjustments to experimental
# data: Reexamining Freedman’s critique." The Annals of Applied Statistics.
# Stat. 7(1): 295-318. doi:10.1214/12-AOAS583.
# https://projecteuclid.org/euclid.aoas/1365527200.

test_that("lm_lin recreates Lin 2013 Table 2", {

  data("alo_star_men")

  ## Table 2
  # Lin uses "classic sandwich," or in our package, HC0

  # unadjusted, Lin est = -0.036, se = 0.158
  expect_equivalent(
    round(
      lm_robust(GPA_year1 ~ sfsp,
                          data = alo_star_men,
                          se_type = 'HC0',
                          coefficient_name = 'sfsp')[, c('est', 'se')],
      3
    ),
    c(-0.036, 0.158)
  )


  # usual adjusted for HS gpa, Lin est = -0.083, se = 0.146
  expect_equivalent(
    unlist(round(
      lm_robust(GPA_year1 ~ sfsp + gpa0,
                          data = alo_star_men,
                          se_type = 'HC0',
                          coefficient_name = 'sfsp')[, c('est', 'se')],
      3
    )),
    c(-0.083, 0.146)
  )

  # interaction adjusted, Lin est = -0.081, se = 0.146
  expect_equivalent(
    unlist(round(
      lm_lin(GPA_year1 ~ sfsp,
                       covariates = ~ gpa0,
                       data = alo_star_men,
                       se_type = 'HC0',
                       coefficient_name = 'sfsp')[, c('est', 'se')],
      3
    )),
    c(-0.081, 0.146)
  )

})


## Table 3 too long to run
rep_table_3 <- FALSE

if (rep_table_3) {

  data("alo_star_men")

  ## Table 3
  samp_dat <- alo_star_men
  its <- 250000
  set.seed(161235)
  check_cover <- function(obj, point = 0) {
    return(obj$ci_lower < point & obj$ci_upper > point)
  }
  ci_dist <- function(obj) {
    return(obj$ci_upper - obj$ci_lower)
  }
  ci_custom <- function(obj) {
    return(list(ci_upper = obj$est + obj$se * 1.96,
                ci_lower = obj$est - obj$se * 1.96))
  }

  ses <- c('HC0', 'HC1', 'HC2', 'HC3')

  ests <- matrix(NA,
                 nrow = its,
                 ncol = 3)
  sd_mats <- cover_mats <- width_mats <-
    array(NA,
          dim = c(its, length(ses), 3))
  for (i in 1:its) {
    samp_dat$sfsp <- sample(samp_dat$sfsp)
    sd_mat <- cover_mat <- width_mat <-
      matrix(NA,
             nrow = length(ses),
             ncol = 3)
    for(j in 1:length(ses)) {
      unadj <- lm_robust(GPA_year1 ~ sfsp,
                         data = samp_dat,
                         se_type = ses[j],
                         coefficient_name = 'sfsp')
      tradadj <- lm_robust(GPA_year1 ~ sfsp + gpa0,
                           data = samp_dat,
                           se_type = ses[j],
                           coefficient_name = 'sfsp')
      intadj <- lm_lin(GPA_year1 ~ sfsp,
                       covariates = ~ gpa0,
                       data = samp_dat,
                       se_type = ses[j],
                       coefficient_name = 'sfsp')

      sd_mat[j, ] <- c(unadj$se, tradadj$se, intadj$se)
      cover_mat[j, ] <- c(check_cover(ci_custom(unadj)),
                          check_cover(ci_custom(tradadj)),
                          check_cover(ci_custom(intadj)))
      width_mat[j, ] <- c(ci_dist(ci_custom(unadj)),
                          ci_dist(ci_custom(tradadj)),
                          ci_dist(ci_custom(intadj)))

    }

    ests[i, ] <- c(unadj$est, tradadj$est, intadj$est)
    sd_mats[i, , ] <- sd_mat
    cover_mats[i, , ] <- cover_mat
    width_mats[i, , ] <- width_mat


    if(i %% 1000 == 0) print(i)
  }


  # Panel A
  colMeans(ests)
  # Panel B
  apply(sd_mats, c(2, 3), mean) - apply(ests, 2, sd)
  # Panel C
  apply(sd_mats, c(2, 3), sd)
  # Panel D, not replicated because he uses normal dist. while we use t dist, all slightly larger
  apply(cover_mats, c(2, 3), mean)
  # Panel E, not replicated because he uses normal dist. while we use t dist, all slightly larger
  apply(width_mats, c(2, 3), mean)
}
