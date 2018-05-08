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
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "classical")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "classical")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  library(lfe)
  lfo <- felm(Y ~ X2 | B + B2 | (X1 ~ Z), data = dat)
  debugonce(lfe:::summary.felm)
  summary(lfo)
  ro[c("r.squared", "adj.r.squared", "fstatistic")]
  rfo[c("r.squared", "adj.r.squared", "fstatistic", "proj_r.squared", "proj_adj.r.squared", "proj_fstatistic")]
  ## HC0
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "HC0")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "HC0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  ## HC1
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "HC1")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "HC1")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  ## HC2
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "HC2")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "HC2")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  ## HC3
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "HC3")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "HC3")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  ## CR0
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR0")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  ## CR stata
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "stata")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "stata")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  ## CR2
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR2")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR2")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

})

test_that("FE matches with weights", {
  ## Classical
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "classical")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "classical")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  ## HC0
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC0")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  ## HC1
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC1")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC1")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  ## HC2
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC2")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC2")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  ## HC3
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC3")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC3")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  ## CR0
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "CR0")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "CR0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  ## CR stata
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "stata")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "stata")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

  ## CR2
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "CR2")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "CR2")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared", "fstatistic")],
    rfo[c("r.squared", "adj.r.squared", "fstatistic")]
  )

})
