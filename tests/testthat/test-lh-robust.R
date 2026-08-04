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
  skip_if_not_installed("car")
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

test_that("lh_robust with clusters works for all se_types (except CR2)", {
  skip_if_not_installed("car")
  for (se_type in cr_se_types_lh) {
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
  skip_if_not_installed("car")

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
  skip_if_not_installed("car")

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
  skip_if_not_installed("car")

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
  skip_if_not_installed("car")

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

# Consistency of cluster results on single coefficient
test_that("lh single coefficient consistency", {
  skip_if_not_installed("car")
  for (se_type in cr_se_types_lh) {
    lhro <-
      tidy(
        lh_robust(
          Y ~ Z * X,
          data = dat,
          clusters = cl,
          linear_hypothesis = "X = 0",
          se_type = se_type
        )
      )
    expect_equal(lhro[3,5], lhro[5,5])
    }
})


# lh test
test_that("returns error when no linear hypothesis is specified", {
  expect_error(lh_robust(Y ~ Z * X, data = dat))
})

# lh test
test_that("returns error when joint  hypothesis is specified", {
  expect_error(lh_robust(Y ~ Z * X, linear_hypothesis = c("Z", "X"), data = dat))
})

# lh test
test_that("returns error when CR2 is specified", {
  expect_error(lh_robust(Y ~ Z * X, se_type = "CR2", linear_hypothesis = c("Z"), data = dat))
})

# lh test
test_that("returns error when default CR2 are called", {
  expect_error(lh_robust(Y ~ Z * X, linear_hypothesis = c("Z"), data = dat, clusters = cl))
})

test_that("issue #405 fixed", {

  library(estimatr)
  nSize = 12
  dat = data.frame(
    x = rnorm(nSize),
    e = rnorm(nSize),
    # Irrelevant clusters for errors
    eg = sample(2, nSize, replace=T)
  )
  dat$z = dat$x + dat$e

  # already worked
  fit <- lh_robust(z~x, data=dat, se_type='HC2', linear_hypothesis='x=0')

  expect_equivalent(fit$lm_robust$conf.low[2],
               fit$lh$conf.low[1])

  # was broken
  fit_2 <- lh_robust(z~x, data=dat, clusters=eg, se_type='stata', linear_hypothesis='x=0')

  expect_equivalent(fit_2$lm_robust$conf.low[2],
                    fit_2$lh$conf.low[1])

})
