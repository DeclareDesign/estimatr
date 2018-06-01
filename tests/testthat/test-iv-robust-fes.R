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

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC0
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "HC0")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "HC0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
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
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "HC1")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "HC1")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
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
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "HC2")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "HC2")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
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
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = "HC3")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "HC3")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
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
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR0")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
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
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "stata")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "stata")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
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
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, se_type = "CR2")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "CR2")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  )

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
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

  expect_equivalent(
    ro$fitted.values,
    rfo$fitted.values
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    rfo[c("r.squared", "adj.r.squared")]
  )

  ## HC0
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC0")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
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
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC1")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC1")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
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
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC2")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC2")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
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
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = "HC3")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC3")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
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
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "CR0")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "CR0")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
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
  ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "stata")
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "stata")

  expect_equivalent(
    tidy(ro)[ro$term %in% c("X1", "X2"), ],
    tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
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
    rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "CR2"),
    "Cannot use `fixed_effects` with weighted"
  )


  # ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = "CR2")
  #
  # expect_equivalent(
  #   tidy(ro)[ro$term %in% c("X1", "X2"), ],
  #   tidy(rfo)[rfo$term %in% c("X1", "X2"), ]
  # )
  #
  # expect_equal(
  #   ro[c("r.squared", "adj.r.squared")],
  #   rfo[c("r.squared", "adj.r.squared")]
  # )

})

test_that("IV FE matches lfe including proj r2", {
  ## unweighted

  ## Classical
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "classical")
  feo <- lfe::felm(Y ~ X2 | B + B2 | (X1 ~ Z), data = dat)
  sfeo <- lfe:::summary.felm(feo)

  expect_equivalent(
    tidy(rfo)[rfo$term %in% c("X1", "X2"), c("estimate", "std.error")],
    data.frame(feo$coefficients, feo$se)[c(2, 1), ]
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
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = "HC1")

  feo <- lfe::felm(Y ~ X2 | B + B2 | (X1 ~ Z), data = dat)
  sfeo <- lfe:::summary.felm(feo, robust = T)

  expect_equivalent(
    tidy(rfo)[rfo$term %in% c("X1", "X2"), c("estimate", "std.error")],
    data.frame(feo$coefficients, sqrt(diag(feo$robustvcv)))[c(2, 1), ]
  )

  expect_equivalent(
    rfo[c("r.squared", "adj.r.squared", "proj_r.squared", "proj_adj.r.squared")],
    sfeo[c("r.squared", "adj.r.squared", "P.r.squared", "P.adj.r.squared")]
  )

  expect_equivalent(
    rfo[["proj_fstatistic"]],
    sfeo[["P.fstat"]][c("F", "df1", "df2")]
  )

  ## CR1stata
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = "stata")
  feo <- lfe::felm(Y ~ X2 | B + B2 | (X1 ~ Z) | cl, data = dat)

  sfeo <- lfe:::summary.felm(feo)

  expect_equivalent(
    tidy(rfo)[rfo$term %in% c("X1", "X2"), c("estimate", "std.error")],
    data.frame(feo$coefficients, sqrt(diag(feo$clustervcv)))[c(2, 1), ]
  )

  expect_equivalent(
    rfo[c("r.squared", "adj.r.squared", "proj_r.squared", "proj_adj.r.squared")],
    sfeo[c("r.squared", "adj.r.squared", "P.r.squared", "P.adj.r.squared")]
  )

  expect_true(
    max(abs(
      rfo[["proj_fstatistic"]] - sfeo[["P.fstat"]][c("F", "df1", "df2")]
    )) < 1e-7
  )

  ## Weighted
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "classical")
  feo <- lfe::felm(Y ~ X2 | B + B2 | (X1 ~ Z), data = dat, weights = dat$w)
  sfeo <- lfe:::summary.felm(feo)

  expect_equivalent(
    tidy(rfo)[rfo$term %in% c("X1", "X2"), c("estimate", "std.error")],
    data.frame(feo$coefficients, feo$se)[c(2, 1), ]
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
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = "HC1")

  feo <- lfe::felm(Y ~ X2 | B + B2 | (X1 ~ Z), data = dat,  weights = dat$w)
  sfeo <- lfe:::summary.felm(feo, robust = T)

  expect_equivalent(
    tidy(rfo)[rfo$term %in% c("X1", "X2"), c("estimate", "std.error")],
    data.frame(feo$coefficients, sqrt(diag(feo$robustvcv)))[c(2, 1), ]
  )

  expect_equivalent(
    rfo[c("r.squared", "adj.r.squared", "proj_r.squared", "proj_adj.r.squared")],
    sfeo[c("r.squared", "adj.r.squared", "P.r.squared", "P.adj.r.squared")]
  )

  expect_equivalent(
    rfo[["proj_fstatistic"]],
    sfeo[["P.fstat"]][c("F", "df1", "df2")]
  )

  ## CR1stata
  rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = "stata")
  feo <- lfe::felm(Y ~ X2 | B + B2 | (X1 ~ Z) | cl, data = dat, weights = dat$w)

  sfeo <- lfe:::summary.felm(feo)

  expect_equivalent(
    tidy(rfo)[rfo$term %in% c("X1", "X2"), c("estimate", "std.error")],
    data.frame(feo$coefficients, sqrt(diag(feo$clustervcv)))[c(2, 1), ]
  )

  expect_equivalent(
    rfo[c("r.squared", "adj.r.squared", "proj_r.squared", "proj_adj.r.squared")],
    sfeo[c("r.squared", "adj.r.squared", "P.r.squared", "P.adj.r.squared")]
  )

  expect_true(
    max(abs(
      rfo[["proj_fstatistic"]] - sfeo[["P.fstat"]][c("F", "df1", "df2")]
    )) < 1e-6
  )

})
