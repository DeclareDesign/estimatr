context("Helper - commarobust + starprep")

test_that("starprep works", {

  fit_1 <- lm(mpg ~ hp, data = mtcars)
  fit_2 <- lm(mpg ~ hp, data = mtcars)

  fit_1_r <- lm_robust(mpg ~ hp, data = mtcars)
  fit_2_r <- lm_robust(mpg ~ hp, data = mtcars)

  library(stargazer)
  stargazer(fit_1, fit_2,
            se = starprep(fit_1_r, fit_2_r),
            p = starprep(fit_1_r, fit_2_r, stat = "p"))

})

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
# In outcome
datmiss <- dat
datmiss$Y[5] <- NA
datmiss$B[1] <- NA

test_that("commarobust works with regular lm", {

  # expect cluster length error
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss)
  expect_error(
    clo <- commarobust(lo, clusters = datmiss$cl, se_type = "CR0"),
    "`clusters` must be the same length as the model data."
  )

  ## Classical
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, se_type = "classical")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss)
  clo <- commarobust(lo, se_type = "classical")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## HC0
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, se_type = "HC0")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss)
  clo <- commarobust(lo, se_type = "HC0")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## HC1
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, se_type = "HC1")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss)
  clo <- commarobust(lo, se_type = "HC1")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## HC2
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, se_type = "HC2")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss)
  clo <- commarobust(lo, se_type = "HC2")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## HC3
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, se_type = "HC3")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss)
  clo <- commarobust(lo, se_type = "HC3")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## CR0
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = datmiss, se_type = "CR0")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss)
  clo <- commarobust(lo, clusters = datmiss$cl[complete.cases(datmiss)], se_type = "CR0")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## CR stata
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = datmiss, se_type = "stata")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss)
  clo <- commarobust(lo, clusters = datmiss$cl[complete.cases(datmiss)], se_type = "stata")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## CR2
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = datmiss, se_type = "CR2")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss)
  clo <- commarobust(lo, clusters = datmiss$cl[complete.cases(datmiss)], se_type = "CR2")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

})

test_that("commarobust works with weighted lm", {
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, weights = w, se_type = "classical")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, se_type = "classical")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## HC0
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, weights = w, se_type = "HC0")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, se_type = "HC0")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## HC1
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, weights = w, se_type = "HC1")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, se_type = "HC1")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## HC2
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, weights = w, se_type = "HC2")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, se_type = "HC2")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## HC3
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, weights = w, se_type = "HC3")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, se_type = "HC3")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## CR0
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = datmiss, weights = w, se_type = "CR0")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, clusters = datmiss$cl[complete.cases(datmiss)], se_type = "CR0")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## CR stata
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = datmiss, weights = w, se_type = "stata")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, clusters = datmiss$cl[complete.cases(datmiss)], se_type = "stata")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## CR2
  ro <- lm_robust(Y ~ Z + X + factor(B) + factor(B2), clusters = cl, data = datmiss, weights = w, se_type = "CR2")
  lo <- lm(Y ~ Z + X + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, clusters = datmiss$cl[complete.cases(datmiss)], se_type = "CR2")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )
})

test_that("commarobust works with dependency, weighted lm", {
  ro <- lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = datmiss, weights = w, se_type = "classical")
  lo <- lm(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, se_type = "classical")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## HC0
  ro <- lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = datmiss, weights = w, se_type = "HC0")
  lo <- lm(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, se_type = "HC0")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## HC1
  ro <- lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = datmiss, weights = w, se_type = "HC1")
  lo <- lm(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, se_type = "HC1")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## HC2
  ro <- lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = datmiss, weights = w, se_type = "HC2")
  lo <- lm(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, se_type = "HC2")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## HC3
  ro <- lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = datmiss, weights = w, se_type = "HC3")
  lo <- lm(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, se_type = "HC3")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## CR0
  ro <- lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), clusters = cl, data = datmiss, weights = w, se_type = "CR0")
  lo <- lm(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, clusters = datmiss$cl[complete.cases(datmiss)], se_type = "CR0")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## CR stata
  ro <- lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), clusters = cl, data = datmiss, weights = w, se_type = "stata")
  lo <- lm(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, clusters = datmiss$cl[complete.cases(datmiss)], se_type = "stata")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )

  ## CR2
  ro <- lm_robust(Y ~ Z + X + Xdup + factor(B) + factor(B2), clusters = cl, data = datmiss, weights = w, se_type = "CR2")
  lo <- lm(Y ~ Z + X + Xdup + factor(B) + factor(B2), data = datmiss, weights = w)
  clo <- commarobust(lo, clusters = datmiss$cl[complete.cases(datmiss)], se_type = "CR2")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  expect_equal(
    ro$fstatistic,
    clo$fstatistic
  )

  expect_equal(
    ro[c("r.squared", "adj.r.squared")],
    clo[c("r.squared", "adj.r.squared")]
  )
})

test_that("Only works with lm, not mlm or glm", {
  expect_error(
    commarobust(glm(vs ~ hp, mtcars, family = binomial), "HC2"),
    "`model` must be an lm object"
  )

  expect_error(
    commarobust(lm(cbind(vs, mpg) ~ hp, data = mtcars), "HC2"),
    "`model` must be an lm object"
  )
})

