context("Estimator - lm_robust, fixed effects")

test_that("test matches lm_robust with dummies", {
  set.seed(43)
  N <- 20
  dat <- data.frame(
    Y = rnorm(N),
    Y2 = rnorm(N),
    Z = rbinom(N, 1, .5),
    X = rnorm(N),
    B = factor(rep(1:2, times = c(8, 12))),
    B2 = factor(rep(1:4, times = c(3, 3, 4, 10))),
    cl = sample(1:4, size = N, replace = T)
  )
  dat$Xdup <- dat$X
  dat$Bdup <- dat$B

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

  expect_equivalent(
    ro[ro$term %in% c("Z"), ],
    rfo[rfo$term %in% c("Z"), ]
  )

  ## HC3
  ro <- tidy(lm_robust(Y ~ Z + factor(B), data = dat, se_type = "HC3"))
  rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~ B, data = dat, se_type = "HC3"))

  expect_equivalent(
    ro[ro$term %in% c("Z"), ],
    rfo[rfo$term %in% c("Z"), ]
  )

  ## CR0
  ro <- tidy(lm_robust(Y ~ Z + factor(B), clusters = cl, data = dat, se_type = "CR0"))
  rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~ B, clusters = cl, data = dat, se_type = "CR0"))

  expect_equivalent(
    ro[ro$term %in% c("Z"), ],
    rfo[rfo$term %in% c("Z"), ]
  )

  ## CR stata
  ro <- tidy(lm_robust(Y ~ Z + factor(B), clusters = cl, data = dat, se_type = "stata"))
  rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~ B, clusters = cl, data = dat, se_type = "stata"))

  expect_equivalent(
    ro[ro$term %in% c("Z"), ],
    rfo[rfo$term %in% c("Z"), ]
  )

  ## CR2
  ro <- tidy(lm_robust(Y ~ Z + factor(B), clusters = cl, data = dat, se_type = "CR2"))
  rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~ B, clusters = cl, data = dat, se_type = "CR2"))

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

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC3

  # Uncomment for perfect fits
  # set.seed(41)
  # N <- 10
  # dat <- data.frame(
  #   Y = rnorm(N),
  #   Z = rbinom(N, 1, .5),
  #   X = rnorm(N),
  #   B = factor(rep(1:2, times = c(4, 6))),
  #   B2 = factor(rep(1:3, times = c(3, 3, 4))),
  #   cl = sample(1:4, size = N, replace = T)
  # )
  # lmo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = dat)
  # summary(lmo)
  # sandwich::vcovHC(lmo, "HC3")

  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "HC3"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC3"))

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

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## Multiple Outcomes

  ## Classical
  ro <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "classical"))
  rfo <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "classical"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC0
  ro <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "HC0"))
  rfo <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC0"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC1
  ro <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "HC1"))
  rfo <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC1"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC2
  ro <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "HC2"))
  rfo <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC2"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC3

  ro <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "HC3"))
  rfo <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC3"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR0
  ro <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR0"))
  rfo <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR0"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR stata
  ro <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "stata"))
  rfo <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "stata"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR2
  ro <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR2"))
  rfo <- tidy(lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR2"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## Collinear factors
  ## Only works for HC2, HC3 for now

  ## Classical
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), data = dat, se_type = "classical"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, data = dat, se_type = "classical"))
  # lo <- felm(Y ~ Z + X|B + B2 + Bdup, data = dat)
  #
  # mtcars$cyl2 <- mtcars$cyl
  # # Doesn't properly deal if collinearity not in first two
  # summary(felm(mpg ~ hp | cyl + am + cyl2, data = mtcars))$coefficients
  # library(estimatr)
  # mtcars$cyl2 <- mtcars$cyl
  # tidy(lm_robust(mpg ~ hp, fixed_effects = ~ cyl + cyl2 + am, data = mtcars))
  #
  # tidy(lm_robust(mpg ~ hp, fixed_effects = ~ cyl + am, data = mtcars))
  # tidy(lm_robust(mpg ~ hp + factor(cyl) + factor(am), data = mtcars))[2,]
  #
  # tidy(lm_robust(mpg ~ hp, fixed_effects = ~ cyl + cyl2 + am, data = mtcars))
  # tidy(lm_robust(mpg ~ hp + factor(cyl) + factor(cyl2)  + factor(am), data = mtcars))[2,]
  # tidy(lm_robust(mpg ~ hp + factor(cyl) + factor(am), data = mtcars))[2,]
  #
  # # LFE does the right thing if the dependency is in the first two
  # summary(felm(mpg ~ hp | cyl + cyl2 + am, data = mtcars))$coefficients
  # tidy(lm_robust(mpg ~ hp + factor(cyl) + factor(cyl3) + factor(am), data = mtcars, se_type = "classical"))[2,]

  # Not the same (denom is wrong)
  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC0
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), data = dat, se_type = "HC0"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, data = dat, se_type = "HC0"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC1
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), data = dat, se_type = "HC1"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, data = dat, se_type = "HC1"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC2
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), data = dat, se_type = "HC2"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, data = dat, se_type = "HC2"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )
  #
  # ## HC3
  #
  # ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), data = dat, se_type = "HC3"))
  # rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, data = dat, se_type = "HC3"))
  #
  # expect_equivalent(
  #   ro[ro$term %in% c("Z", "X"), ],
  #   rfo[rfo$term %in% c("Z", "X"), ]
  # )
  #
  # ## CR0
  # ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), clusters = cl, data = dat, se_type = "CR0"))
  # rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, clusters = cl, data = dat, se_type = "CR0"))
  #
  # expect_equivalent(
  #   ro[ro$term %in% c("Z", "X"), ],
  #   rfo[rfo$term %in% c("Z", "X"), ]
  # )
  #
  # ## CR stata
  # ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), clusters = cl, data = dat, se_type = "stata"))
  # rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, clusters = cl, data = dat, se_type = "stata"))
  #
  # expect_equivalent(
  #   ro[ro$term %in% c("Z", "X"), ],
  #   rfo[rfo$term %in% c("Z", "X"), ]
  # )
  #
  # ## CR2
  # ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), clusters = cl, data = dat, se_type = "CR2"))
  # rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, clusters = cl, data = dat, se_type = "CR2"))
  #
  # expect_equivalent(
  #   ro[ro$term %in% c("Z", "X"), ],
  #   rfo[rfo$term %in% c("Z", "X"), ]
  # )

  ## Collinear covariates

  ## Classical
  ro <- tidy(lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = dat, se_type = "classical"))
  rfo <- tidy(lm_robust(Y ~ Z + X + Xdup, fixed_effects = ~ B + B2, data = dat, se_type = "classical"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X", "Xdup"), ],
    rfo[rfo$term %in% c("Z", "X", "Xdup"), ]
  )

  ## HC0
  ro <- tidy(lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = dat, se_type = "HC0"))
  rfo <- tidy(lm_robust(Y ~ Z + X + Xdup, fixed_effects = ~ B + B2, data = dat, se_type = "HC0"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X", "Xdup"), ],
    rfo[rfo$term %in% c("Z", "X", "Xdup"), ]
  )

  ## HC1
  ro <- tidy(lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = dat, se_type = "HC1"))
  rfo <- tidy(lm_robust(Y ~ Z + X + Xdup, fixed_effects = ~ B + B2, data = dat, se_type = "HC1"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X", "Xdup"), ],
    rfo[rfo$term %in% c("Z", "X", "Xdup"), ]
  )

  ## HC2
  ro <- tidy(lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = dat, se_type = "HC2"))
  rfo <- tidy(lm_robust(Y ~ Z + X + Xdup, fixed_effects = ~ B + B2, data = dat, se_type = "HC2"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X", "Xdup"), ],
    rfo[rfo$term %in% c("Z", "X", "Xdup"), ]
  )

  ## HC3

  ro <- tidy(lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = dat, se_type = "HC3"))
  rfo <- tidy(lm_robust(Y ~ Z + X + Xdup, fixed_effects = ~ B + B2, data = dat, se_type = "HC3"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X", "Xdup"), ],
    rfo[rfo$term %in% c("Z", "X", "Xdup"), ]
  )

  ## CR0
  ro <- tidy(lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR0"))
  rfo <- tidy(lm_robust(Y ~ Z + X + Xdup, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR0"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X", "Xdup"), ],
    rfo[rfo$term %in% c("Z", "X", "Xdup"), ]
  )

  ## CR stata
  ro <- tidy(lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "stata"))
  rfo <- tidy(lm_robust(Y ~ Z + X + Xdup, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "stata"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X", "Xdup"), ],
    rfo[rfo$term %in% c("Z", "X", "Xdup"), ]
  )

  ## CR2
  ro <- tidy(lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR2"))
  rfo <- tidy(lm_robust(Y ~ Z + X + Xdup, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR2"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X", "Xdup"), ],
    rfo[rfo$term %in% c("Z", "X", "Xdup"), ]
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
