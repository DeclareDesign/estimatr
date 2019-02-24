context("Estimator - lh_robust")

test_that("lh_robust works with classical", {

     for (se_type in se_types) {
       lhro <- tidy(lh_robust(mpg ~ cyl + disp, mtcars,  linearHypothesis = "cyl + disp = 0", se_type =  se_type))
       lmro <- lm_robust(mpg ~ cyl + disp, mtcars, se_type = se_type)
       linHyp <-  car::linearHypothesis(lmro,  hypothesis.matrix =  "cyl + disp = 0")

       expect_equivalent(
         lhro$std.error[ lhro$term == "cyl + disp = 0"],
         sqrt(as.numeric(attr(linHyp , "vcov")))
       )

     }


})




#
# test_that("lh_robust with clusters works for all se_types ", {
#   set.seed(42)
#   N <- 20
#   dat <- data.frame(
#     Y = rnorm(N),
#     Y2 = rnorm(N),
#     Z = rbinom(N, 1, .5),
#     X = rnorm(N),
#     B = factor(rep(1:2, times = c(8, 12))),
#     cl = sample(1:4, size = N, replace = T),
#     w = runif(N)
#   )
#
#
#   for (se_type in cr_se_types) {
#     lhro <- tidy(lh_robust(Y ~ Z*X, data = dat, clusters =  cl, linearHypothesis = "Z + Z:X = 0", se_type = se_type))
#     lmro <- lm_robust(Y ~ Z*X,  data = dat, se_type = se_type, clusters =  cl)
#     linHyp <-  car::linearHypothesis(lmro,  hypothesis.matrix =  "Z + Z:X = 0")
#     expect_equal(
#       lhro$std.error[ lhro$term == "Z + Z:X = 0"],
#       sqrt(as.numeric(attr(linHyp , "vcov")))
#     )
#   }
# })
#
#


test_that("lh_robust matches lm_robust with fixed effects", {
N <- 40
  dat <- data.frame(
    Y = rnorm(N),
    Y2 = rnorm(N),
    Z = rbinom(N, 1, .5),
    X = rnorm(N),
    B = factor(rep(1:2, times = c(8, 12))),
    cl = sample(1:4, size = N, replace = T),
    w = runif(N)
  )

  lhro <- lh_robust(Y ~ Z*X, data = dat, fixed_effects = ~ B, linearHypothesis = c("Z + Z:X = 0"))
  lmro_lh <- tidy(lhro$lm_robust)
  lmro <- tidy(lm_robust(Y ~ Z*X, data = dat, fixed_effects = ~ B))
  linHyp <-  car::linearHypothesis(lmro,  hypothesis.matrix =  "Z + Z:X = 0")

  expect_equal(lmro,  lmro_lh)
  expect_equal(
    lhro$std.error[ lhro$term == "Z + Z:X = 0"],
    sqrt(as.numeric(attr(linHyp , "vcov")))
  )

})


test_that("lh_robust matches lm_robust with weights", {

  set.seed(42)
  N <- 20
  dat <- data.frame(
    Y = rnorm(N),
    Y2 = rnorm(N),
    Z = rbinom(N, 1, .5),
    X = rnorm(N),
    B = factor(rep(1:2, times = c(8, 12))),
    cl = sample(1:4, size = N, replace = T),
    w = runif(N)
  )

  lhro <- lh_robust(Y ~ Z*X, data = dat,  weights = w, linearHypothesis = c("Z + Z:X = 0"))
  lmro_lh <- tidy(lhro$lm_robust)
  lmro <- tidy(lm_robust(Y ~ Z*X, data = dat,  weights = w))
  linHyp <-  car::linearHypothesis(lmro,  hypothesis.matrix =  "Z + Z:X = 0")

  expect_equal(lmro,  lmro_lh)
  expect_equal(
    lhro$std.error[ lhro$term == "Z + Z:X = 0"],
    sqrt(as.numeric(attr(linHyp , "vcov")))
  )

})

test_that("lh_robust matches lm_robust with subsetted data.frame", {
  set.seed(42)
  N <- 20
  dat <- data.frame(
    Y = rnorm(N),
    Y2 = rnorm(N),
    Z = rbinom(N, 1, .5),
    X = rnorm(N),
    B = factor(rep(1:2, times = c(8, 12))),
    cl = sample(1:4, size = N, replace = T),
    w = runif(N)
  )
  lhro <- lh_robust(Y ~ Z*X, data = dat,  subset =  B == 1, linearHypothesis = c("Z + Z:X = 0"))
  lmro_lh <- tidy(lhro$lm_robust)
  lmro <- tidy(lm_robust(Y ~ Z*X, data = dat,  subset =  B == 1))
  linHyp <-  car::linearHypothesis(lmro,  hypothesis.matrix =  "Z + Z:X = 0")

  expect_equal(lmro,  lmro_lh)
  expect_equal(
    lhro$std.error[ lhro$term == "Z + Z:X = 0"],
    sqrt(as.numeric(attr(linHyp , "vcov")))
  )

})

