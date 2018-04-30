context("Estimator - lm_robust, fixed effects")

test_that("test matches lm_robust with dummies", {
  set.seed(41)
  N <- 10
  dat <- data.frame(
    Y = rnorm(N),
    Z = rbinom(N, 1, .5),
    X = rnorm(N),
    B = factor(rep(1:2, times = c(4, 6))),
    B2 = factor(rep(1:3, times = c(3, 3, 4))),
    cl = sample(1:4, size = N, replace = T)
  )

  ## One FE, one covar

  ## Classical
  ro <- tidy(lm_robust(Y ~ Z + factor(B), data = dat, se_type = "classical"))
  rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~ B, data = dat, se_type = "classical"))

  ## TODO f stats and stuff

  expect_equivalent(
    ro[ro$term %in% c("Z"), ],
    rfo[rfo$term %in% c("Z"), ]
  )

  ## HC0
  ro <- tidy(lm_robust(Y ~ Z + factor(B), data = dat, se_type = "HC0"))
  rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~ B, data = dat, se_type = "HC0"))

  expect_equivalent(
    ro[ro$term %in% c("Z"), ],
    rfo[rfo$term %in% c("Z"), ]
  )

  ## HC1
  ro <- tidy(lm_robust(Y ~ Z + factor(B), data = dat, se_type = "HC1"))
  rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~ B, data = dat, se_type = "HC1"))

  expect_equivalent(
    ro[ro$term %in% c("Z"), ],
    rfo[rfo$term %in% c("Z"), ]
  )

  ## HC2
  ro <- tidy(lm_robust(Y ~ Z + factor(B), data = dat, se_type = "HC2"))
  rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~ B, data = dat, se_type = "HC2"))

  # Not equivalent
  expect_equivalent(
    ro[ro$term %in% c("Z"), ],
    rfo[rfo$term %in% c("Z"), ]
  )

  ## HC3
  ro <- tidy(lm_robust(Y ~ Z + factor(B), data = dat, se_type = "HC3"))
  rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~ B, data = dat, se_type = "HC3"))

  # Not equivalent
  expect_equivalent(
    ro[ro$term %in% c("Z"), ],
    rfo[rfo$term %in% c("Z"), ]
  )

  ## CR0
  ro <- tidy(lm_robust(Y ~ Z + factor(B), clusters = cl, data = dat, se_type = "CR0"))
  rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR0"))

  expect_equivalent(
    ro[ro$term %in% c("Z"), ],
    rfo[rfo$term %in% c("Z"), ]
  )

  ## CR stata
  ro <- tidy(lm_robust(Y ~ Z + factor(B), clusters = cl, data = dat, se_type = "stata"))
  rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "stata"))

  expect_equivalent(
    ro[ro$term %in% c("Z"), ],
    rfo[rfo$term %in% c("Z"), ]
  )

  ## CR2
  ro <- tidy(lm_robust(Y ~ Z + factor(B), clusters = cl, data = dat, se_type = "CR2"))
  rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR2"))

  # Not equivalent, only works if clusters are factors
  expect_equivalent(
    ro[ro$term %in% c("Z"), ],
    rfo[rfo$term %in% c("Z"), ]
  )

  ## Multiple FEs, multiple covars

  ## Classical
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "classical"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "classical"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC0
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "HC0"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC0"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC1
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "HC1"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC1"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC2
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "HC2"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC2"))

  # Not equivalent
  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC3
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "HC3"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC3"))

  # Not equivalent
  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR0
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR0"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR0"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR stata
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "stata"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "stata"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR2
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR2"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR2"))

  # Not equivalent, only works if clusters are factors
  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

})


test_that("test matches stata absorb", {

  rfo <- tidy(lm_robust(mpg ~ hp, mtcars, fixed_effects = ~ carb, se_type = "classical")) # areg mpg hp, absorb(carb)
  rfo
  rfo <- tidy(lm_robust(mpg ~ hp, mtcars, fixed_effects = ~ carb, se_type = "HC1")) # areg mpg hp, absorb(carb) rob
  rfo
  rfo <- tidy(lm_robust(mpg ~ hp, mtcars, fixed_effects = ~ carb, clusters = cyl, se_type = "stata")) # areg mpg hp, absorb(carb) cl(cyl)
  rfo

  rfo <- tidy(lm_robust(mpg ~ hp, mtcars, fixed_effects = ~ carb, weights = wt, se_type = "classical")) # areg mpg hp [aweight=w], absorb(carb)
  rfo
  rfo <- tidy(lm_robust(mpg ~ hp, mtcars, fixed_effects = ~ carb, se_type = "HC1")) # areg mpg hp [aweight=w], absorb(carb) rob
  rfo
  rfo <- tidy(lm_robust(mpg ~ hp, mtcars, fixed_effects = ~ carb, clusters = cyl, se_type = "stata")) # areg mpg hp [aweight=w], absorb(carb) cl(cyl)
  rfo


})
