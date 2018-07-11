context("Estimator - lm_robust, fixed effects")

set.seed(43)
N <- 20
dat <- data.frame(
  Y = rnorm(N),
  Y2 = rnorm(N),
  Z = rbinom(N, 1, .5),
  X = rnorm(N),
  B = factor(rep(1:2, times = c(8, 12))),
  B2 = factor(rep(1:4, times = c(3, 3, 4, 10))),
  cl = sample(1:4, size = N, replace = T),
  w = runif(N)
)
dat$Xdup <- dat$X
dat$Bdup <- dat$B

test_that("FE matches lm_robust with dummies", {
  ## One FE, one covar

  ## Classical
  ro <- tidy(lm_robust(Y ~ Z + factor(B), data = dat, se_type = "classical"))
  rfo <- tidy(lm_robust(Y ~ Z, fixed_effects = ~ B, data = dat, se_type = "classical"))

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
})

test_that("FE matches with multiple FEs and covars", {
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

  # Uncomment for perfect fits which reveal problems for our estimators
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
})

test_that("FE matches with weights", { ## Classical
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = dat, weights = w, se_type = "classical"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "classical"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC0
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC0"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC0"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC1
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC1"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC1"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC2
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC2"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC2"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC3
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC3"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC3"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR0
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "CR0"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "CR0"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR stata
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "stata"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "stata"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR2
  expect_error(
    rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "CR2")),
    "Cannot use `fixed_effects` with weighted"
  )

  # ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "CR2"))
  #
  # expect_equivalent(
  #   ro[ro$term %in% c("Z", "X"), ],
  #   rfo[rfo$term %in% c("Z", "X"), ]
  # )
})

test_that("FEs work with multiple outcomes", {
  ## Multiple Outcomes

  ## Classical
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "classical")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "classical")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC0
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "HC0")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC1
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "HC1")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC1")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC2
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "HC2")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC2")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC3
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, se_type = "HC3")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC3")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## CR0
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR0")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## CR stata
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "stata")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "stata")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## CR2
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR2")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR2")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## multiple outcomes with weights
  ## Classical
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, weights = w, se_type = "classical")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "classical")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC0
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC0")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC1
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC1")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC1")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC2
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC2")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC2")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC3
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC3")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC3")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## CR0
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "CR0")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "CR0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## CR stata
  ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "stata")
  rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "stata")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## CR2
  expect_error(
    rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "CR2"),
    "Cannot use `fixed_effects` with weighted"
  )
  # ro <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "CR2")
  #
  # expect_equivalent(
  #   tidy(ro)[ro$term %in% c("Z", "X"), ],
  #   tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  # )
  #
  # expect_true(
  #   max(
  #     unlist(ro[c("r.squared", "adj.r.squared")]) -
  #       unlist(rfo[c("r.squared", "adj.r.squared")])
  #   ) < 1e-7
  # )

})


test_that("FEs work with missingness", {

  # In outcome
  datmiss <- dat
  datmiss$Y[5] <- NA
  datmiss$B[1] <- NA

  ## Classical
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, se_type = "classical")
  expect_warning(
    rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = datmiss, se_type = "classical"),
    "Some observations have missingness in the fixed effects but not in the outcome or covariates."
  )

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC0
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, se_type = "HC0")
  expect_warning(
    rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = datmiss, se_type = "HC0"),
    "Some observations have missingness in the fixed effects but not in the outcome or covariates."
)

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC1
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, se_type = "HC1")
  expect_warning(
    rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = datmiss, se_type = "HC1"),
    "Some observations have missingness in the fixed effects but not in the outcome or covariates."
  )

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC2
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, se_type = "HC2")
  expect_warning(
    rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = datmiss, se_type = "HC2"),
    "Some observations have missingness in the fixed effects but not in the outcome or covariates."
  )

  # Check to make sure with only one FE when missingness works
  expect_warning(
    lm_robust(Y ~ Z + X, fixed_effects = ~ B, data = datmiss, se_type = "HC2"),
    "Some observations have missingness in the fixed effects but not in the outcome or covariates."
  )

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC3
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, se_type = "HC3")
  expect_warning(
    rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = datmiss, se_type = "HC3"),
    "Some observations have missingness in the fixed effects but not in the outcome or covariates."
  )
  lfo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss)

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    tidy(rfo)[rfo$term %in% c("Z", "X"), c("std.error")],
    sqrt(diag(sandwich::vcovHC(lfo, type = "HC3"))[2:3])
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## CR0
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = datmiss, se_type = "CR0")
  expect_warning(
    rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = datmiss, se_type = "CR0"),
    "Some observations have missingness in the fixed effects but not in the outcome or covariates."
  )

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## CR stata
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = datmiss, se_type = "stata")
  expect_warning(
    rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = datmiss, se_type = "stata"),
    "Some observations have missingness in the fixed effects but not in the outcome or covariates."
  )

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## CR2
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = datmiss, se_type = "CR2")
  expect_warning(
    rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = datmiss, se_type = "CR2"),
    "Some observations have missingness in the fixed effects but not in the outcome or covariates."
  )

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )
})

test_that("FEs handle collinear FEs", {
  ## Collinear factors
  ## Classical
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), data = dat, se_type = "classical"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, data = dat, se_type = "classical"))
  # lo <- lfe::felm(Y ~ Z + X|B + B2 + Bdup, data = dat)
  #
  # mtcars$cyl2 <- mtcars$cyl
  # # Doesn't properly deal if collinearity not in first two
  # lfe:::summary.felm(lfe::felm(mpg ~ hp | cyl + am + cyl2, data = mtcars))$coefficients
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
  # lfe:::summary.felm(lfe::felm(mpg ~ hp | cyl + cyl2 + am, data = mtcars))$coefficients
  # tidy(lm_robust(mpg ~ hp + factor(cyl) + factor(cyl3) + factor(am), data = mtcars, se_type = "classical"))[2,]

  expect_equivalent(
    ro$estimate[ro$term %in% c("Z", "X")],
    rfo$estimate[rfo$term %in% c("Z", "X")]
  )

  # SE is wrong because we count the FEs incorrectly for the finite sample correction
  expect_false(
    all(
      ro$std.error[ro$term %in% c("Z", "X")] ==
        rfo$std.error[rfo$term %in% c("Z", "X")]
    )
  )

  ## HC0
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), data = dat, se_type = "HC0"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, data = dat, se_type = "HC0"))

  # SE is right for HC0 (as there is no finite sample correction), but DoF is still wrong
  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), c("estimate", "std.error")],
    rfo[rfo$term %in% c("Z", "X"), c("estimate", "std.error")]
  )

  # SE is right, but FE is wrong because we count the FEs incorrectly for the finite sample correction
  expect_false(
    all(
      ro$df[ro$term %in% c("Z", "X")] ==
        rfo$df[rfo$term %in% c("Z", "X")]
    )
  )

  ## HC1
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), data = dat, se_type = "HC1"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, data = dat, se_type = "HC1"))

  expect_equivalent(
    ro$estimate[ro$term %in% c("Z", "X")],
    rfo$estimate[rfo$term %in% c("Z", "X")]
  )

  # SE is wrong because we count the FEs incorrectly for the finite sample correction
  expect_false(
    all(
      ro$std.error[ro$term %in% c("Z", "X")] ==
        rfo$std.error[rfo$term %in% c("Z", "X")]
    )
  )

  ## HC2
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), data = dat, se_type = "HC2"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, data = dat, se_type = "HC2"))

  # HC2 or HC3 work because we can get the collinearity in the FEs for free as we have to invert
  # UtU anyways (where U is cbind(X, FE_dummy_mat))
  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## HC3
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), data = dat, se_type = "HC3"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, data = dat, se_type = "HC3"))

  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR0
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), clusters = cl, data = dat, se_type = "CR0"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, clusters = cl, data = dat, se_type = "CR0"))

  # DoF for CR0 works, unlike HC0, because our DoF for CR0 is N_clust - 1, not N - total_rank
  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )

  ## CR stata
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), clusters = cl, data = dat, se_type = "stata"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, clusters = cl, data = dat, se_type = "stata"))

  expect_equivalent(
    ro$estimate[ro$term %in% c("Z", "X")],
    rfo$estimate[rfo$term %in% c("Z", "X")]
  )

  # SE is wrong because we count the FEs incorrectly for the finite sample correction
  expect_false(
    all(
      ro$std.error[ro$term %in% c("Z", "X")] ==
        rfo$std.error[rfo$term %in% c("Z", "X")]
    )
  )

  ## CR2
  ro <- tidy(lm_robust(Y ~ Z + X + factor(B) + factor(Bdup) + factor(B2), clusters = cl, data = dat, se_type = "CR2"))
  rfo <- tidy(lm_robust(Y ~ Z + X, fixed_effects = ~ B + Bdup + B2, clusters = cl, data = dat, se_type = "CR2"))

  # As with HC2, we get the count of FEs here for free
  expect_equivalent(
    ro[ro$term %in% c("Z", "X"), ],
    rfo[rfo$term %in% c("Z", "X"), ]
  )
})

test_that("FEs work with collinear covariates", {
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

  # write.csv(mtcars,
  #           file = 'tests/testthat/mtcars.csv',
  #           row.names = F)

  stata_ests <- read.table(
    "stata-fe-ests.txt",
    col.names = c("model", "var", "fstat"),
    stringsAsFactors = FALSE
  )

  mtcars$w <- mtcars$drat / 5

  estimatr_mat <- matrix(NA, 6, 1)

  rfo <- lm_robust(mpg ~ hp, mtcars, fixed_effects = ~ carb, se_type = "classical") # areg mpg hp, absorb(carb)
  estimatr_mat[1, ] <- c(rfo$std.error ^ 2)

  rfo <- tidy(lm_robust(mpg ~ hp, mtcars, fixed_effects = ~ carb, se_type = "HC1")) # areg mpg hp, absorb(carb) rob
  estimatr_mat[2, ] <- c(rfo$std.error ^ 2)

  rfo <- tidy(lm_robust(mpg ~ hp, mtcars, fixed_effects = ~ carb, clusters = cyl, se_type = "stata")) # areg mpg hp, absorb(carb) cl(cyl)
  estimatr_mat[3, ] <- c(rfo$std.error ^ 2)


  rfo <- tidy(lm_robust(mpg ~ hp, mtcars, fixed_effects = ~ carb, weights = w, se_type = "classical")) # areg mpg hp [aweight=w], absorb(carb)
  estimatr_mat[4, ] <- c(rfo$std.error ^ 2)

  rfo <- tidy(lm_robust(mpg ~ hp, mtcars, fixed_effects = ~ carb, weights = w, se_type = "HC1")) # areg mpg hp [aweight=w], absorb(carb) rob
  estimatr_mat[5, ] <- c(rfo$std.error ^ 2)

  rfo <- tidy(lm_robust(mpg ~ hp, mtcars, fixed_effects = ~ carb, weights = w, clusters = cyl, se_type = "stata")) # areg mpg hp [aweight=w], absorb(carb) cl(cyl)
  estimatr_mat[6, ] <- c(rfo$std.error ^ 2)

  expect_equal(
    estimatr_mat[, 1],
    stata_ests[, 2]
  )

})


test_that("FEs give correct projected F-stats", {

  feo <- lfe::felm(Y ~ Z + X | B + B2, data = dat)
  sfeo <- lfe:::summary.felm(feo)
  sfeor <- lfe:::summary.felm(feo, robust = TRUE)

  cfeo <- lfe::felm(Y ~ Z + X | B + B2 | 0 | cl, data = dat)
  sfeoc <- lfe:::summary.felm(cfeo, robust = TRUE)

  # classical
  rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "classical")

  expect_equivalent(
    tidy(rfo)[rfo$term %in% c("Z", "X"), c("estimate", "std.error", "statistic", "p.value")],
    as.data.frame(sfeo$coefficients[, 1:4])
  )

  expect_equivalent(
    rfo[c("r.squared", "adj.r.squared", "proj_r.squared", "proj_adj.r.squared")],
    sfeo[c("r.squared", "adj.r.squared", "P.r.squared", "P.adj.r.squared")]
  )

  expect_equivalent(
    rfo[["proj_fstatistic"]],
    sfeo[["P.fstat"]][c("F", "df1", "df2")]
  )

  ## HC1
  rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, data = dat, se_type = "HC1")

  expect_equivalent(
    tidy(rfo)[rfo$term %in% c("Z", "X"), c("estimate", "std.error", "statistic", "p.value")],
    as.data.frame(sfeor$coefficients[, 1:4])
  )

  expect_equivalent(
    rfo[c("r.squared", "adj.r.squared", "proj_r.squared", "proj_adj.r.squared")],
    sfeor[c("r.squared", "adj.r.squared", "P.r.squared", "P.adj.r.squared")]
  )

  expect_equivalent(
    rfo[["proj_fstatistic"]],
    sfeor[["P.fstat"]][c("F", "df1", "df2")]
  )

  ## CR stata
  rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "stata")

  # Different pval because lfe doesn't use J-1 as it's DoF
  expect_equivalent(
    tidy(rfo)[rfo$term %in% c("Z", "X"), c("estimate", "std.error")],
    as.data.frame(sfeoc$coefficients[, c(1, 2)])
  )

  expect_equivalent(
    rfo[c("r.squared", "adj.r.squared", "proj_r.squared", "proj_adj.r.squared")],
    sfeoc[c("r.squared", "adj.r.squared", "P.r.squared", "P.adj.r.squared")]
  )

  expect_equivalent(
    rfo[["proj_fstatistic"]],
    sfeoc[["P.fstat"]][c("F", "df1", "df2")]
  )

})

test_that("FE matches lm_robust with one block", {
  # In outcome
  datmiss <- dat
  datmiss$Y[5] <- NA
  datmiss$oneB <- "A"

  ## Classical
  ro <- lm_robust(Y ~ Z + X, data = datmiss, se_type = "classical")
  rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ oneB, data = datmiss, se_type = "classical")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC0
  ro <- lm_robust(Y ~ Z + X, data = datmiss, se_type = "HC0")
  rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ oneB, data = datmiss, se_type = "HC0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC1
  ro <- lm_robust(Y ~ Z + X, data = datmiss, se_type = "HC1")
  rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ oneB, data = datmiss, se_type = "HC1")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC2
  ro <- lm_robust(Y ~ Z + X, data = datmiss, se_type = "HC2")
  rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ oneB, data = datmiss, se_type = "HC2")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC3
  ro <- lm_robust(Y ~ Z + X, data = datmiss, se_type = "HC3")
  rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ oneB, data = datmiss, se_type = "HC3")
  lfo <- lm(Y ~ Z + X, data = datmiss)

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equivalent(
    tidy(rfo)[rfo$term %in% c("Z", "X"), c("std.error")],
    sqrt(diag(sandwich::vcovHC(lfo, type = "HC3"))[2:3])
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## CR0
  ro <- lm_robust(Y ~ Z + X, clusters = cl, data = datmiss, se_type = "CR0")
  rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ oneB, clusters = cl, data = datmiss, se_type = "CR0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## CR stata
  ro <- lm_robust(Y ~ Z + X, clusters = cl, data = datmiss, se_type = "stata")
  rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ oneB, clusters = cl, data = datmiss, se_type = "stata")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## CR2
  ro <- lm_robust(Y ~ Z + X, clusters = cl, data = datmiss, se_type = "CR2")
  rfo <- lm_robust(Y ~ Z + X, fixed_effects = ~ oneB, clusters = cl, data = datmiss, se_type = "CR2")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("Z", "X"), ],
    tidy(rfo)[rfo$term %in% c("Z", "X"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## Error when combined with other blocks
  expect_error(
    lm_robust(Y ~ Z + X, fixed_effects = ~ oneB + B, data = datmiss),
    "Can't have a fixed effect with only one group AND multiple fixed effect variables"
  )
})

test_that("Handle perfect fits appropriately", {
  dat$Bsingle <- c(1, 2, rep(3:4, each = 9))
  rfo <- lm_robust(Y ~ X, fixed_effects = ~ Bsingle, data = dat)
  ro <- lm_robust(Y ~ X + factor(Bsingle), data = dat)

  expect_equivalent(
    tidy(rfo)[rfo$term == "X", ],
    tidy(ro)[ro$term == "X", ]
  )
})
