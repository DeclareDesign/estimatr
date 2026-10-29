library(estimatr)

# df_type = "satterthwaite" with HC2: the Bell and McCaffrey (2002) degrees of
# freedom, taken from CR2 with every observation its own cluster. The variance
# is the default HC2 variance, so the only thing these tests can find moved is
# the degrees of freedom and what is computed from them.

d <- ext_data_ols()

# ---- the variance is untouched ----

test_that("df_type = 'satterthwaite' leaves the HC2 estimates and variance identical", {
  default <- lm_robust(y ~ x + z, data = d)
  bm <- lm_robust(y ~ x + z, data = d, df_type = "satterthwaite")
  expect_identical(bm$coefficients, default$coefficients)
  expect_identical(bm$vcov, default$vcov)
  expect_identical(bm$std.error, default$std.error)
  expect_false(isTRUE(all.equal(bm$df, default$df)))
})

test_that("df_type = NULL and 'residual' are the default HC2 fit", {
  default <- lm_robust(y ~ x + z, data = d)
  expect_identical(lm_robust(y ~ x + z, data = d, df_type = NULL)$df, default$df)
  expect_identical(lm_robust(y ~ x + z, data = d, df_type = "residual")$df, default$df)
})

test_that("df_type = 'satterthwaite' moves the p-values and intervals with the degrees of freedom", {
  bm <- lm_robust(y ~ x + z, data = d, df_type = "satterthwaite")
  tstat <- bm$coefficients / bm$std.error
  expect_equal(bm$p.value, 2 * pt(abs(tstat), bm$df, lower.tail = FALSE))
  expect_equal(bm$conf.high,
               bm$coefficients + qt(0.975, bm$df) * bm$std.error)
})

# ---- closed form ----

# For a difference in means the HC2 variance is s1^2/n1 + s0^2/n0, and Bell and
# McCaffrey's degrees of freedom are Satterthwaite's under a homoskedastic
# working model: Welch's formula with the two sample variances set equal.
test_that("the two-group degrees of freedom are the closed form", {
  set.seed(343)
  n1 <- 18
  n0 <- 12
  g <- data.frame(y = rnorm(n1 + n0), z = rep(0:1, c(n0, n1)))
  fit <- lm_robust(y ~ z, data = g, df_type = "satterthwaite")
  closed <- (1 / n1 + 1 / n0)^2 /
    ((1 / n1)^2 / (n1 - 1) + (1 / n0)^2 / (n0 - 1))
  expect_equal(unname(fit$df["z"]), closed)
})

# ---- #289 ----

# A contrast that only one group's observations identify gets the degrees of
# freedom of the fit on that group, whether it is a coefficient of the
# interacted model or a linear hypothesis on it.
test_that("an interacted contrast gets the subset fit's degrees of freedom (#289)", {
  skip_if_not_installed("car")
  set.seed(343)
  dat <- data.frame(Y = rnorm(40), X = rbinom(40, 1, 0.5), Z = rbinom(40, 1, 0.5))

  full <- lh_robust(Y ~ X * Z, data = dat, linear_hypothesis = "Z + X:Z = 0",
                    df_type = "satterthwaite")
  sub1 <- lm_robust(Y ~ Z, data = subset(dat, X == 1), df_type = "satterthwaite")
  sub0 <- lm_robust(Y ~ Z, data = subset(dat, X == 0), df_type = "satterthwaite")

  expect_equal(unname(full$lh$coefficients), unname(sub1$coefficients["Z"]))
  expect_equal(full$lh$df, unname(sub1$df["Z"]))
  expect_equal(unname(full$lm_robust$df["Z"]), unname(sub0$df["Z"]))
})

# ---- clubSandwich ----

test_that("HC2 Satterthwaite degrees of freedom match clubSandwich on lm", {
  skip_if_not_installed("clubSandwich")
  m <- lm(y ~ x + z, data = d)
  fit <- lm_robust(y ~ x + z, data = d, df_type = "satterthwaite")
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = seq_len(nrow(d)),
                                test = "Satterthwaite")
  expect_equal(unname(fit$df), ct$df_Satt, tolerance = LIVE_TOL)
  expect_equal(unname(fit$std.error), ct$SE, tolerance = LIVE_TOL)
})

test_that("HC2 Satterthwaite degrees of freedom match clubSandwich on ivreg", {
  skip_if_not_installed("clubSandwich")
  skip_if_not_installed("ivreg")
  di <- ext_data_fe_iv()
  fit <- iv_robust(y ~ en + x | inst + x, data = di, df_type = "satterthwaite")
  mi <- ivreg::ivreg(y ~ en + x | inst + x, data = di)
  ct <- clubSandwich::coef_test(mi, vcov = "CR2", cluster = seq_len(nrow(di)),
                                test = "Satterthwaite")
  expect_equal(unname(fit$df), ct$df_Satt, tolerance = LIVE_TOL)
  expect_equal(unname(fit$std.error), ct$SE, tolerance = LIVE_TOL)
})

# Weighted HC2 is clubSandwich's CR2 with singleton clusters at
# `inverse_var = TRUE`, and so are its degrees of freedom. At the default
# `inverse_var = FALSE` the working model is the identity, and both the
# variance and the degrees of freedom differ.
test_that("weighted HC2 Satterthwaite degrees of freedom match clubSandwich at inverse_var = TRUE", {
  skip_if_not_installed("clubSandwich")
  m <- lm(y ~ x + z, data = d, weights = w)
  fit <- lm_robust(y ~ x + z, data = d, weights = w, df_type = "satterthwaite")
  V <- clubSandwich::vcovCR(m, cluster = seq_len(nrow(d)), type = "CR2",
                            inverse_var = TRUE)
  ct <- clubSandwich::coef_test(m, vcov = V, test = "Satterthwaite")
  expect_equal(unname(fit$df), ct$df_Satt, tolerance = LIVE_TOL)
  expect_equal(unname(fit$std.error), ct$SE, tolerance = LIVE_TOL)
})

# clubSandwich has no `inverse_var` for ivreg. The degrees of freedom depend on
# the design alone, which in the second stage is the first-stage fitted values,
# so they are those of the weighted regression on that design.
test_that("weighted iv_robust degrees of freedom are those of the second-stage design", {
  di <- ext_data_fe_iv()
  fit <- iv_robust(y ~ en + x | inst + x, data = di, weights = w,
                   df_type = "satterthwaite")
  di$en_hat <- fitted(lm(en ~ inst + x, data = di, weights = w))
  second <- lm_robust(y ~ en_hat + x, data = di, weights = w,
                      df_type = "satterthwaite")
  expect_equal(unname(fit$df), unname(second$df), tolerance = LIVE_TOL)
})

test_that("lh_robust HC2 Satterthwaite degrees of freedom match clubSandwich's contrast", {
  skip_if_not_installed("car")
  skip_if_not_installed("clubSandwich")
  m <- lm(y ~ x + z, data = d, weights = w)
  fit <- lh_robust(y ~ x + z, data = d, weights = w, df_type = "satterthwaite",
                   linear_hypothesis = "x - z = 0")
  V <- clubSandwich::vcovCR(m, cluster = seq_len(nrow(d)), type = "CR2",
                            inverse_var = TRUE)
  ct <- clubSandwich::linear_contrast(m, vcov = V, contrasts = matrix(c(0, 1, -1), 1),
                                      test = "Satterthwaite")
  expect_equal(fit$lh$df, ct$df, tolerance = LIVE_TOL)
})

# ---- fixed effects ----

test_that("absorbed fixed effects give the degrees of freedom of explicit dummies", {
  di <- ext_data_fe_iv()
  absorbed <- lm_robust(y ~ x, data = di, fixed_effects = ~g,
                        df_type = "satterthwaite")
  dummies <- lm_robust(y ~ x + g, data = di, df_type = "satterthwaite")
  expect_equal(absorbed$df[["x"]], dummies$df[["x"]], tolerance = LIVE_TOL)
  expect_equal(absorbed$std.error[["x"]], dummies$std.error[["x"]],
               tolerance = LIVE_TOL)

  two_way <- lm_robust(y ~ x, data = di, fixed_effects = ~ g + cl,
                       df_type = "satterthwaite")
  two_way_dummies <- lm_robust(y ~ x + g + factor(cl), data = di,
                               df_type = "satterthwaite")
  expect_equal(two_way$df[["x"]], two_way_dummies$df[["x"]], tolerance = LIVE_TOL)

  absorbed_iv <- iv_robust(y ~ en | inst, data = di, fixed_effects = ~g,
                           df_type = "satterthwaite")
  dummies_iv <- iv_robust(y ~ en + g | inst + g, data = di,
                          df_type = "satterthwaite")
  expect_equal(absorbed_iv$df[["en"]], dummies_iv$df[["en"]], tolerance = LIVE_TOL)
})

# ---- pass-through and ci = FALSE ----

test_that("lm_lin passes df_type through", {
  dl <- d
  dl$x_c <- dl$x - mean(dl$x)
  lin <- lm_lin(y ~ z, covariates = ~x, data = d, df_type = "satterthwaite")
  direct <- lm_robust(y ~ z * x_c, data = dl, df_type = "satterthwaite")
  expect_equal(unname(lin$df), unname(direct$df), tolerance = LIVE_TOL)
})

test_that("df_type = 'satterthwaite' with ci = FALSE returns no inference", {
  fit <- lm_robust(y ~ x + z, data = d, df_type = "satterthwaite", ci = FALSE)
  default <- lm_robust(y ~ x + z, data = d, ci = FALSE)
  expect_identical(fit$std.error, default$std.error)
  expect_true(all(is.na(fit$p.value)))
  expect_identical(fit$df, default$df)
})

# ---- refusals ----

test_that("df_type refuses every pairing it does not implement", {
  expect_error(lm_robust(y ~ x, data = d, se_type = "HC1", df_type = "satterthwaite"),
               "not available with `se_type = \"HC1\"`")
  expect_error(lm_robust(y ~ x, data = d, se_type = "HC3", df_type = "satterthwaite"),
               "not available with `se_type = \"HC3\"`")
  expect_error(lm_robust(y ~ x, data = d, clusters = cl, se_type = "CR0",
                         df_type = "satterthwaite"),
               "not available with `se_type = \"CR0\"`")
  expect_error(lm_robust(y ~ x, data = d, clusters = cl, se_type = "CR2",
                         df_type = "residual"),
               "not available with `se_type = \"CR2\"`")
  expect_error(lm_robust(y ~ x, data = d, df_type = "welch"),
               "must be NULL, 'residual', or 'satterthwaite'")
  expect_error(lm_robust(y ~ x, data = d, df_type = c("residual", "satterthwaite")),
               "must be NULL, 'residual', or 'satterthwaite'")

  di <- ext_data_fe_iv()
  expect_error(lm_robust(y ~ x, data = di, weights = w, fixed_effects = ~g,
                         df_type = "satterthwaite"),
               "cannot be combined with `weights` and `fixed_effects`")
})

test_that("CR2 with df_type = 'satterthwaite' is the default CR2 fit", {
  default <- lm_robust(y ~ x, data = d, clusters = cl)
  named <- lm_robust(y ~ x, data = d, clusters = cl, df_type = "satterthwaite")
  expect_identical(named$df, default$df)
})
