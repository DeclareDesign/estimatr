context("Estimator - iv_robust, fixed effects")

set.seed(43)
N <- 20
dat <- data.frame(
  Y = rnorm(N),
  X1 = rnorm(N),
  X2 = rnorm(N),
  Z = rbinom(N, 1, .5),
  B = factor(rep(1:2, times = c(8, 12))),
  B2 = factor(rep(1:4, times = c(3, 3, 4, 10))),
  cl = sample(1:4, size = N, replace = T),
  w = runif(N)
)
dat$Xdup <- dat$X
dat$Bdup <- dat$B

test_that("FE matches with multiple FEs and covars", {
  ## Classical
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "classical"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "classical"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC0
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "HC0"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "HC0"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC1
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "HC1"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "HC1"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC2
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "HC2"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "HC2"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC3
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "HC3"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "HC3"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR0
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR0"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR0"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR stata
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "stata"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "stata"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR2
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR2"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR2"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )
})

test_that("FE matches with weights", {
  ## Classical
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "classical"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "classical"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC0
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC0"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC0"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC1
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC1"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC1"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC2
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC2"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC2"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC3
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC3"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC3"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR0
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "CR0"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "CR0"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR stata
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "stata"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "stata"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR2
  ro <- tidy(iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "CR2"))
  rfo <- tidy(iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "CR2"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )
})
