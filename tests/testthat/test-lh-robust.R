context("Estimator - lh_robust")
set.seed(40)
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

# se tests
test_that("lh_robust works with all se types", {
  for (se_type in se_types) {
    lhro <-
      tidy(
        lh_robust(
          mpg ~ cyl + disp,
          data = mtcars,
          linear_hypothesis = "cyl + disp = 0",
          se_type =  se_type
        )
      )
    lmro <-
      lm_robust(mpg ~ cyl + disp, data = mtcars, se_type = se_type)
    linHyp <-
      car::linearHypothesis(lmro,  hypothesis.matrix =  "cyl + disp = 0")

    expect_equal(lhro$std.error[lhro$term == "cyl + disp = 0"],
                 sqrt(as.numeric(attr(linHyp , "vcov"))))
  }
})

test_that("lh_robust with clusters works for all se_types ", {
  for (se_type in cr_se_types) {
    lhro <-
      tidy(
        lh_robust(
          Y ~ Z * X,
          data = dat,
          clusters = cl,
          linear_hypothesis = "Z + Z:X = 0",
          se_type = se_type
        )
      )
    lmro <-
      lm_robust(Y ~ Z * X,
                data = dat,
                se_type = se_type,
                clusters =  cl)
    linHyp <-
      car::linearHypothesis(lmro,  hypothesis.matrix =  "Z + Z:X = 0")
    expect_equal(lhro$std.error[lhro$term == "Z + Z:X = 0"],
                 sqrt(as.numeric(attr(linHyp , "vcov"))))
  }
})

test_that("lh_robust matches lm_robust with fixed effects", {
  lhro <-
    lh_robust(
      Y ~ Z * X,
      data = dat,
      fixed_effects = ~ B,
      linear_hypothesis = c("Z + Z:X = 0")
    )
  lmro <- lm_robust(Y ~ Z * X, data = dat, fixed_effects = ~ B)
  linHyp <-
    car::linearHypothesis(lmro,  hypothesis.matrix =  "Z + Z:X = 0")
  tidy_lhro <- tidy(lhro)

  expect_equal(tidy_lhro$std.error[tidy_lhro$term == "Z + Z:X = 0"],
               sqrt(as.numeric(attr(linHyp , "vcov"))))

})

test_that("lh_robust matches lm_robust with weights", {
  lhro <-
    lh_robust(
      Y ~ Z * X,
      data = dat,
      weights = w,
      linear_hypothesis = c("Z + Z:X = 0")
    )
  tidy_lhro <- tidy(lhro)
  lmro <- lm_robust(Y ~ Z * X, data = dat,  weights = w)
  linHyp <-
    car::linearHypothesis(lmro,  hypothesis.matrix =  "Z + Z:X = 0")


  expect_equal(tidy_lhro$std.error[tidy_lhro$term == "Z + Z:X = 0"],
               sqrt(as.numeric(attr(linHyp , "vcov"))))

})

test_that("lh_robust matches lm_robust with subsetted data.frame", {
  lhro <-
    lh_robust(Y ~ Z * X,
              data = dat,
              subset = B == 1,
              linear_hypothesis = c("Z + Z:X = 0"))
  tidy_lhro <- tidy(lhro)
  lmro <- lm_robust(Y ~ Z * X, data = dat,  subset =  B == 1)
  linHyp <-
    car::linearHypothesis(lmro,  hypothesis.matrix =  "Z + Z:X = 0")


  expect_equal(tidy_lhro$std.error[tidy_lhro$term == "Z + Z:X = 0"],
               sqrt(as.numeric(attr(linHyp , "vcov"))))

})

test_that("lh_robust matches lm_robust with subsetted data.frame", {
  lhro <-
    lh_robust(Y ~ Z * X,
              data = dat,
              subset =  B == 1,
              linear_hypothesis = c("Z + Z:X = 0"))
  tidy_lhro <- tidy(lhro)
  lmro <- lm_robust(Y ~ Z * X, data = dat,  subset =  B == 1)
  linHyp <-
    car::linearHypothesis(lmro,  hypothesis.matrix =  "Z + Z:X = 0")


  expect_equal(tidy_lhro$std.error[tidy_lhro$term == "Z + Z:X = 0"],
               sqrt(as.numeric(attr(linHyp , "vcov"))))

})

# lh test
test_that("returns error when no linear hypothesis is specified", {
  expect_error(lh_robust(Y ~ Z * X, data = dat))
})
