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

  for (se_type in se_types) {
    ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, se_type = se_type)
    rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, se_type = se_type)


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

    # weights
    ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), data = dat, weights = w, se_type = se_type)
    rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, data = dat, weights = w, se_type = se_type)

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
  }

  for (se_type in cr_se_types) {
    ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, se_type = se_type)
    rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, se_type = se_type)

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

    # weights
    if (se_type %in% c("CR2", "CR3")) {
      expect_error(
        rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = se_type),
        "Cannot use `fixed_effects` with weighted CR2"
      )
    } else {
      ro <- iv_robust(Y ~ X1 + X2 + factor(B) + factor(B2) | Z + X2 + factor(B) + factor(B2), clusters = cl, data = dat, weights = w, se_type = se_type)
      rfo <- iv_robust(Y ~ X1 + X2 | Z + X2, fixed_effects = ~ B + B2, clusters = cl, data = dat, weights = w, se_type = se_type)

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
    }
  }
})

test_that("IV FE warns about diagnostics", {

  expect_warning(
    iv_robust(mpg ~ hp | wt, data = mtcars, fixed_effects = cyl, diagnostics = TRUE),
    "Will not return `diagnostics` if `fixed_effects` are used."
  )

})
