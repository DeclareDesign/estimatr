context("Helper - commarobust + starprep")

test_that("starprep works", {

  fit_1 <- lm(mpg ~ hp, data = mtcars)
  fit_2 <- lm(mpg ~ hp, data = mtcars)

  fit_1_r <- lm_robust(mpg ~ hp, data = mtcars)
  fit_2_r <- lm_robust(mpg ~ hp, data = mtcars)

  expect_output(
    stargazer::stargazer(fit_1, fit_2,
              type = "text",
              se = starprep(fit_1_r, fit_2),
              p = starprep(fit_1_r, fit_2, stat = "p.value")),
    "\\(0\\.015\\)\\s+\\(0\\.015\\)"
  )

  expect_output(
    stargazer::stargazer(fit_1, fit_2,
              type = "text",
              ci.custom = starprep(fit_1_r, fit_2, stat = "ci")),
    "\\(25\\.620\\, 34\\.578\\)\\s+\\(25\\.620\\, 34\\.578\\)"
  )

})

set.seed(43)
N <- 480
dat <- data.frame(
  Z = rbinom(N, 1, .5),
  X = rnorm(N),
  B = factor(rep(1:2, times = c(8, 12))),
  cl = sample(1:(N/4), size = N, replace = T),
  w = runif(N)
)
dat$Y <- dat$Z + dat$X + rnorm(N)
dat$Y2 = dat$Z + rnorm(N)
dat$Xdup <- dat$X
dat$Bdup <- dat$B
# In outcome
datmiss <- dat
datmiss$Y[5] <- NA
datmiss$B[1] <- NA

test_that("commarobust works with regular lm", {

  # expect cluster length error
  lo <- lm(Y ~ Z + X + factor(B), data = datmiss)
  expect_error(
    clo <- commarobust(lo, clusters = datmiss$cl, se_type = "CR0"),
    "`clusters` must be the same length as the model data."
  )

  ## Test unclustered SEs
  for (se_type in se_types) {
    ro <- lm_robust(Y ~ Z + X + factor(B), data = datmiss, se_type = se_type)
    lo <- lm(Y ~ Z + X + factor(B), data = datmiss)
    clo <- commarobust(lo, se_type = se_type)

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
  }

  ## Test clustered SEs
  for (se_type in cr_se_types) {

    ro <- lm_robust(Y ~ Z + X + factor(B), clusters = cl, data = datmiss, se_type = se_type)
    lo <- lm(Y ~ Z + X + factor(B), data = datmiss)
    clo <- commarobust(lo, clusters = datmiss$cl[complete.cases(datmiss)], se_type = se_type)

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
  }

  # Works with character, factor, and numeric clusters
  datmiss$cl_char <- sample(letters, size = nrow(datmiss), replace = TRUE)
  datmiss$cl_num <- sample(rnorm(3), size = nrow(datmiss), replace = TRUE)
  datmiss$cl_fac <- as.factor(datmiss$cl_char)

  ro <- lm_robust(Y ~ Z + X + factor(B), clusters = cl_char, data = datmiss, se_type = "CR2")
  lo <- lm(Y ~ Z + X + factor(B), data = datmiss)
  clo <- commarobust(lo, clusters = datmiss$cl_char[complete.cases(datmiss)], se_type = "CR2")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  ro <- lm_robust(Y ~ Z + X + factor(B), clusters = cl_num, data = datmiss, se_type = "CR2")
  lo <- lm(Y ~ Z + X + factor(B), data = datmiss)
  clo <- commarobust(lo, clusters = datmiss$cl_num[complete.cases(datmiss)], se_type = "CR2")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

  ro <- lm_robust(Y ~ Z + X + factor(B), clusters = cl_fac, data = datmiss, se_type = "CR2")
  lo <- lm(Y ~ Z + X + factor(B), data = datmiss)
  clo <- commarobust(lo, clusters = datmiss$cl_fac[complete.cases(datmiss)], se_type = "CR2")

  expect_equal(
    tidy(ro),
    tidy(clo)
  )

})

test_that("commarobust works with weighted lm", {
  # Test unclustered SEs
  for (se_type in se_types) {
    ro <- lm_robust(Y ~ Z + X + factor(B), data = datmiss, weights = w, se_type = se_type)
    lo <- lm(Y ~ Z + X + factor(B), data = datmiss, weights = w)
    clo <- commarobust(lo, se_type = se_type)

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
  }

  ## Test clustered SEs
  for (se_type in cr_se_types) {

    ro <- lm_robust(Y ~ Z + X + factor(B), clusters = cl, data = datmiss, weights = w, se_type = se_type)
    lo <- lm(Y ~ Z + X + factor(B), data = datmiss, weights = w)
    clo <- commarobust(lo, clusters = datmiss$cl[complete.cases(datmiss)], se_type = se_type)


    max(abs(clo$vcov - ro$vcov))
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
  }

})

test_that("commarobust works with dependency, weighted lm", {
  check_obj <- function(ro, clo, x) {
    if (x != "call") {
      print(x)
      expect_equal(ro[[x]], clo[[x]])
    }
  }

  for (se_type in se_types) {
    ro <- lm_robust(Y ~ Z + X + Xdup + factor(B), data = datmiss, weights = w, se_type = se_type)
    lo <- lm(Y ~ Z + X + Xdup + factor(B), data = datmiss, weights = w)
    clo <- commarobust(lo, se_type = se_type)

    capture_output(sapply(names(ro), check_obj, ro = ro, clo = clo))

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
  }

  for (se_type in cr_se_types) {
    ro <- lm_robust(Y ~ Z + X + Xdup + factor(B), clusters = cl, data = datmiss, weights = w, se_type = se_type)
    lo <- lm(Y ~ Z + X + Xdup + factor(B), data = datmiss, weights = w)
    clo <- commarobust(lo, clusters = datmiss$cl[complete.cases(datmiss)], se_type = se_type)

    capture_output(sapply(names(ro), check_obj, ro = ro, clo = clo))

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
  }
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

