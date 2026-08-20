library(estimatr)

# HC2 and HC3 with absorbed fixed effects.
#
# These used to be refused, on the grounds that they need hat values from the
# full [dummies | X] design and absorption only leaves the demeaned ones. That
# is true for two or more FE factors and false for one, where the hat value
# splits exactly:
#
#   h_ii = h_ii(demeaned X) + w_i / sum(w over i's group)
#
# so the answer is exact and costs one vector. The tests below pin the identity
# against the two references that matter: writing the dummies out by hand, and
# estimatr itself.

configs <- ref_data_leverage()

# ---- against explicit dummies ----

for (nm in names(configs)) {
  for (se in c("HC2", "HC3")) {
    test_that(paste0("absorbed FE equals explicit dummies: ", nm, ", ", se), {
      d <- configs[[nm]]
      fe  <- lm_robust(y ~ x + z, fixed_effects = ~ g, data = d, se_type = se)
      dum <- lm_robust(y ~ x + z + factor(g), data = d, se_type = se)
      expect_equal(unname(fe$std.error),
                   unname(dum$std.error[c("x", "z")]))
      expect_equal(unname(fe$coefficients),
                   unname(dum$coefficients[c("x", "z")]))
    })
  }
}

test_that("the identity survives weights", {
  for (nm in names(configs)) {
    d <- configs[[nm]]
    for (se in c("HC2", "HC3")) {
      fe  <- lm_robust(y ~ x + z, fixed_effects = ~ g, data = d,
                       se_type = se, weights = wts)
      dum <- lm_robust(y ~ x + z + factor(g), data = d,
                       se_type = se, weights = wts)
      expect_equal(unname(fe$std.error), unname(dum$std.error[c("x", "z")]),
                   info = paste(nm, se))
    }
  }
})

# ---- against estimatr ----

test_that("absorbed FE agrees with estimatr for HC2 and HC3", {
  for (nm in names(configs)) {
    d <- configs[[nm]]
    for (se in c("HC2", "HC3")) {
      a <- lm_robust(y ~ x + z, fixed_effects = ~ g, data = d, se_type = se)
      b <- ref(paste0("lev_", nm, "_", se))
      expect_equal(unname(a$std.error), unname(b$std.error),
                   tolerance = 1e-10, info = paste(nm, se))
    }
  }
})

test_that("absorbed FE agrees with estimatr under weights", {
  d <- configs$unbalanced
  for (se in c("HC2", "HC3")) {
    a <- lm_robust(y ~ x + z, fixed_effects = ~ g, data = d,
                   se_type = se, weights = wts)
    b <- ref(paste0("levw_unbalanced_", se))
    expect_equal(unname(a$std.error), unname(b$std.error), tolerance = 1e-10)
  }
})

# ---- defaults ----

test_that("one-way FE keeps the package default of HC2", {
  # It used to fall back to "stata" because HC2 was unaffordable. It is not
  # unaffordable any more, and HC2 is also what estimatr returns for this call,
  # so falling back would be a silent disagreement with the released package.
  d <- configs$balanced
  expect_equal(lm_robust(y ~ x + z, fixed_effects = ~ g, data = d)$se_type, "HC2")
})

test_that("two-way FE still falls back, since the identity does not hold there", {
  set.seed(3)
  d <- data.frame(y = rnorm(300), x = rnorm(300),
                  a = factor(rep(1:15, each = 20)),
                  b = factor(rep(1:20, times = 15)))
  # Reported as "HC1" rather than "stata": the fit resolves the alias before
  # storing it, which is also what `se_type = "stata"` reports with no FE.
  expect_equal(lm_robust(y ~ x, fixed_effects = ~ a + b, data = d)$se_type, "HC1")
})

# ---- what is still refused, and why ----

test_that("two-way FE refuses HC2 and HC3 rather than approximating", {
  set.seed(3)
  d <- data.frame(y = rnorm(300), x = rnorm(300),
                  a = factor(rep(1:15, each = 20)),
                  b = factor(rep(1:20, times = 15)))
  expect_error(lm_robust(y ~ x, fixed_effects = ~ a + b, data = d, se_type = "HC2"),
               "cannot be used with `fixed_effects`")
  expect_error(lm_robust(y ~ x, fixed_effects = ~ a + b, data = d, se_type = "HC3"),
               "cannot be used with `fixed_effects`")
})

test_that("CR2 is still refused even with one FE factor", {
  # CR2's adjustment comes from cluster-level blocks of the hat matrix, not
  # from h_ii, so the one-way identity says nothing about it.
  d <- configs$balanced
  d$cl <- factor(rep(1:10, length.out = nrow(d)))
  expect_error(
    lm_robust(y ~ x + z, fixed_effects = ~ g, clusters = cl, data = d, se_type = "CR2"),
    "cannot be used with `fixed_effects`"
  )
})

test_that("the two-way error still names the dummy escape hatch", {
  set.seed(3)
  d <- data.frame(y = rnorm(300), x = rnorm(300),
                  a = factor(rep(1:15, each = 20)),
                  b = factor(rep(1:20, times = 15)))
  expect_error(lm_robust(y ~ x, fixed_effects = ~ a + b, data = d, se_type = "HC2"),
               "factor\\(fe_var\\)")
})

# ---- the leverage identity itself ----

test_that("the FE leverage vector is the group weight share", {
  d <- configs$unbalanced
  md <- list(
    fixed_effects = d["g"],
    outcome = as.matrix(d$y),
    design_matrix = cbind("(Intercept)" = 1, x = d$x),
    weights = d$wts,
    terms = terms(y ~ x, data = d)
  )
  out <- estimatr:::demean_fes(md)
  expect_equal(out$fe_leverage, d$wts / ave(d$wts, d$g, FUN = sum))
})

test_that("no leverage vector is produced for two-way FE", {
  set.seed(3)
  d <- data.frame(y = rnorm(300), x = rnorm(300),
                  a = factor(rep(1:15, each = 20)),
                  b = factor(rep(1:20, times = 15)))
  md <- list(
    fixed_effects = d[c("a", "b")],
    outcome = as.matrix(d$y),
    design_matrix = cbind("(Intercept)" = 1, x = d$x),
    weights = NULL,
    terms = terms(y ~ x, data = d)
  )
  expect_null(estimatr:::demean_fes(md)$fe_leverage)
})
