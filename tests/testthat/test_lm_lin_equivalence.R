library(estimatr)

# lm_lin against an explicit lm_robust construction.
#
# Lin's estimator is a regression of the outcome on the treatment, the
# covariates centred at their mean, and the interaction of the two. That is a
# specification a user could write out by hand, so the whole of lm_lin can be
# checked against it, and it is worth checking that way because lm_lin builds
# the design internally: a mistake in the centring or in the interaction
# expansion would show up nowhere else.
#
# The weighted cells are the reason this file exists. With weights the centring
# point becomes the weighted mean, and the reference fixture in
# test_vs_estimatr.R compares 2.0 against 1.0.6 rather than against the
# specification, so an error inherited from 1.0.6 would pass there unnoticed.
# The grid below crosses the treatment kinds lm_lin accepts with the covariate
# kinds, weighted and unweighted, clustered and not.

lin_test_data <- function() {
  set.seed(1)
  n <- 300
  d <- data.frame(
    y  = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n),
    w  = runif(n, 0.5, 2),
    cl = rep(1:30, 10)
  )
  # Drawn rather than cycled: a treatment and a factor covariate that both
  # alternate with period two are collinear, and the fit is then rank deficient
  # for reasons that have nothing to do with lm_lin.
  d$zb <- rbinom(n, 1, 0.5)
  d$zf <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  d$zn <- sample(c(0, 2, 5), n, replace = TRUE)
  d$xf <- factor(sample(c("p", "q"), n, replace = TRUE))
  # lm_lin treats a treatment with more than two values as categorical, so the
  # hand-built reference has to do the same.
  d$znf <- factor(d$zn)
  d
}

# The Lin specification, written out.
lin_by_hand <- function(d, z, covariates, weighted, clustered, se_type) {
  cov_mat <- model.matrix(reformulate(covariates), d)[, -1, drop = FALSE]
  centre <- if (weighted) {
    apply(cov_mat, 2, weighted.mean, d$w)
  } else {
    colMeans(cov_mat)
  }
  centred <- sweep(cov_mat, 2, centre)
  colnames(centred) <- paste0(colnames(cov_mat), "_c")
  dd <- cbind(d, as.data.frame(centred))

  form <- reformulate(
    c(z, colnames(centred), paste0(z, ":", colnames(centred))),
    response = "y"
  )
  args <- list(formula = form, data = dd, se_type = se_type)
  if (weighted) args$weights <- dd$w
  if (clustered) args$clusters <- dd$cl
  do.call(lm_robust, args)
}

lin_grid <- expand.grid(
  z         = c("zb", "zf", "zn"),
  covariate = c("one", "two", "factor"),
  weighted  = c(FALSE, TRUE),
  clustered = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

test_that("lm_lin equals the hand-built Lin specification across the grid", {
  d <- lin_test_data()

  for (i in seq_len(nrow(lin_grid))) {
    g <- lin_grid[i, ]
    covariates <- switch(g$covariate,
      one    = "x1",
      two    = c("x1", "x2"),
      factor = c("x1", "xf")
    )
    se_type <- if (g$clustered) "CR2" else "HC2"
    label <- paste(g$z, g$covariate, if (g$weighted) "weighted" else "unweighted",
                   if (g$clustered) "clustered" else "iid", sep = "/")

    args <- list(
      formula = reformulate(g$z, response = "y"),
      covariates = reformulate(covariates),
      data = d, se_type = se_type
    )
    if (g$weighted) args$weights <- d$w
    if (g$clustered) args$clusters <- d$cl
    fit <- do.call(lm_lin, args)

    ref <- lin_by_hand(d, if (g$z == "zn") "znf" else g$z,
                       covariates, g$weighted, g$clustered, se_type)
    # lm_lin names the levels of a numeric treatment after the original
    # variable; the hand-built version goes through a factor of the same name.
    ref_names <- sub("^znf", "zn", names(coef(ref)))

    nms <- names(coef(fit))
    expect_setequal(nms, ref_names)

    # The whole coefficient vector and the whole variance matrix, not just the
    # treatment coefficient. An error in the interaction block would leave the
    # treatment coefficient untouched in the balanced cells.
    expect_equal(unname(coef(fit)), unname(coef(ref)[match(nms, ref_names)]),
                 tolerance = 1e-10, label = paste0(label, ": coefficients"))
    ref_vcov <- ref$vcov[match(nms, ref_names), match(nms, ref_names)]
    expect_equal(unname(fit$vcov), unname(ref_vcov),
                 tolerance = 1e-10, label = paste0(label, ": vcov"))
  }
})

test_that("every se_type agrees with the hand-built specification when weighted", {
  d <- lin_test_data()
  for (se in c("classical", "HC0", "HC1", "HC2", "HC3")) {
    fit <- lm_lin(y ~ zf, covariates = ~ x1 + x2, data = d, weights = w, se_type = se)
    ref <- lin_by_hand(d, "zf", c("x1", "x2"), weighted = TRUE,
                       clustered = FALSE, se_type = se)
    nms <- names(coef(fit))
    expect_equal(unname(coef(fit)), unname(coef(ref)[nms]),
                 tolerance = 1e-10, label = paste0(se, ": coefficients"))
    expect_equal(unname(fit$vcov), unname(ref$vcov[nms, nms]),
                 tolerance = 1e-10, label = paste0(se, ": vcov"))
  }
})

# ---- what the weights are allowed to do to the centring point ----

test_that("scaling every weight by a constant changes nothing", {
  d <- lin_test_data()
  a <- lm_lin(y ~ zb, covariates = ~ x1 + x2, data = d, weights = w)
  d$w <- d$w * 10
  b <- lm_lin(y ~ zb, covariates = ~ x1 + x2, data = d, weights = w)
  expect_equal(coef(a), coef(b), tolerance = 1e-12)
  expect_equal(a$std.error, b$std.error, tolerance = 1e-12)
  expect_equal(a$scaled_center, b$scaled_center, tolerance = 1e-12)
})

test_that("a zero weight is the same as removing the row", {
  d <- lin_test_data()
  set.seed(3)
  dropped <- sample(nrow(d), 40)
  d$w[dropped] <- 0
  a <- suppressWarnings(lm_lin(y ~ zb, covariates = ~ x1 + x2, data = d, weights = w))
  b <- lm_lin(y ~ zb, covariates = ~ x1 + x2, data = d[-dropped, ], weights = w)
  expect_equal(coef(a), coef(b), tolerance = 1e-10)
  expect_equal(a$std.error, b$std.error, tolerance = 1e-10)
  expect_equal(a$scaled_center, b$scaled_center, tolerance = 1e-10)
})

test_that("the centring point is the weighted mean over complete cases only", {
  # The case that would be wrong if the covariate mean were taken before the
  # rows with a missing outcome were dropped. The two candidate centres differ
  # by about 2e-2 here, so the comparison discriminates between them.
  d <- lin_test_data()
  set.seed(4)
  missing_y <- sample(nrow(d), 30)
  d$y[missing_y] <- NA
  fit <- lm_lin(y ~ zb, covariates = ~ x1 + x2, data = d, weights = w)

  complete <- d[!is.na(d$y), ]
  expect_equal(
    unname(fit$scaled_center),
    c(weighted.mean(complete$x1, complete$w), weighted.mean(complete$x2, complete$w)),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(
    unname(fit$scaled_center),
    c(weighted.mean(d$x1, d$w), weighted.mean(d$x2, d$w)),
    tolerance = 1e-4
  )))

  refit <- lm_lin(y ~ zb, covariates = ~ x1 + x2, data = complete, weights = w)
  expect_equal(coef(fit), coef(refit), tolerance = 1e-12)
  expect_equal(fit$std.error, refit$std.error, tolerance = 1e-12)
})

test_that("a missing covariate drops the row before centring", {
  d <- lin_test_data()
  set.seed(5)
  d$x1[sample(nrow(d), 25)] <- NA
  fit <- lm_lin(y ~ zb, covariates = ~ x1 + x2, data = d, weights = w)
  refit <- lm_lin(y ~ zb, covariates = ~ x1 + x2, data = d[!is.na(d$x1), ], weights = w)
  expect_equal(coef(fit), coef(refit), tolerance = 1e-12)
  expect_equal(fit$std.error, refit$std.error, tolerance = 1e-12)
  expect_equal(fit$scaled_center, refit$scaled_center, tolerance = 1e-12)
})

test_that("weighted and unweighted centring differ when the weights are informative", {
  # Guards the grid above: if lm_lin ignored the weights when centring, every
  # weighted cell would still agree with a hand-built reference that also
  # ignored them. It does not, and this says so.
  d <- lin_test_data()
  d$w <- exp(d$x1)
  fit_w <- lm_lin(y ~ zb, covariates = ~ x1 + x2, data = d, weights = w)
  fit_u <- lm_lin(y ~ zb, covariates = ~ x1 + x2, data = d)
  expect_gt(abs(fit_w$scaled_center[["x1"]] - fit_u$scaled_center[["x1"]]), 0.1)
})

# The factorial recipe in ?lm_lin's details (#131): the fully interacted
# regression on centered terms and lm_lin() on the four cells are the same
# model, and their interaction is the difference-in-differences of each
# cell's own covariate regression evaluated at the overall means.
test_that("a two-by-two factorial's interaction is the DiD of Lin-adjusted cell means", {
  set.seed(343)
  n <- 200
  d <- data.frame(z1 = rbinom(n, 1, 0.4), z2 = rbinom(n, 1, 0.6),
                  x1 = rnorm(n), x2 = runif(n))
  d$y <- with(d, z1 + 0.5 * z2 + 0.8 * z1 * z2 +
                x1 * (1 + z1 - z2 + 2 * z1 * z2) + x2 * (1 - z1) + rnorm(n))
  d <- transform(d, z1_c = z1 - mean(z1), z2_c = z2 - mean(z2),
                 x1_c = x1 - mean(x1), x2_c = x2 - mean(x2),
                 cell = factor(paste0(z1, z2)))

  at_means <- data.frame(x1 = mean(d$x1), x2 = mean(d$x2))
  m <- vapply(split(d, d$cell), function(s) {
    unname(predict(lm(y ~ x1 + x2, data = s), newdata = at_means))
  }, numeric(1))
  did <- m[["11"]] - m[["10"]] - m[["01"]] + m[["00"]]
  p2 <- mean(d$z2)
  z1_effect <- p2 * (m[["11"]] - m[["01"]]) + (1 - p2) * (m[["10"]] - m[["00"]])

  full <- lm_robust(y ~ z1_c * z2_c * (x1_c + x2_c), data = d)
  expect_equal(unname(coef(full)["z1_c:z2_c"]), did, tolerance = 1e-10)
  expect_equal(unname(coef(full)["z1_c"]), z1_effect, tolerance = 1e-10)

  cells <- lm_lin(y ~ cell, covariates = ~ x1 + x2, data = d)
  L <- c(cell11 = 1, cell10 = -1, cell01 = -1)
  V <- vcov(cells)[names(L), names(L)]
  expect_equal(sum(coef(cells)[names(L)] * L), did, tolerance = 1e-10)
  expect_equal(sqrt(drop(L %*% V %*% L)), unname(full$std.error["z1_c:z2_c"]),
               tolerance = 1e-10)

  # The shortcut documented as wrong with raw covariates, right with
  # pre-centered ones.
  raw <- lm_lin(y ~ z1, covariates = ~ z2 * (x1 + x2), data = d)
  pre <- lm_lin(y ~ z1, covariates = ~ z2 * (x1_c + x2_c), data = d)
  expect_gt(abs(coef(raw)[["z1:z2_c"]] - did), 0.05)
  expect_equal(coef(pre)[["z1:z2_c"]], did, tolerance = 1e-10)
})
