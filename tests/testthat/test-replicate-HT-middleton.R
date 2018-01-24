context("Verification - HT matches Joel Middleton code")

test_that("We match Joel's estimator", {

  # Code from Joel Middleton
  n <- 400

  # simple random assignment
  d.00 <- diag(rep(1, n))
  dmat <- rbind(
    cbind(d.00, -d.00)
    , cbind(-d.00, d.00)
  )
  d.tilde <- diag(rep(2, 2 * n))
  d.00 <- matrix(rep(-0.001251564, n ^ 2), ncol = n) + diag(rep(1 + 0.001251564, n))
  pmat <- diag(rep(.5, 2 * n)) %*% (dmat + 1) %*% diag(rep(.5, 2 * n))
  pmat.true <- pmat
  pmat[pmat == 0] <- 1
  d.tilde.wt <- d.tilde / pmat

  # complete random assignment
  dmat.CR <- round(rbind(
    cbind(d.00, -d.00)
    , cbind(-d.00, d.00)
  ), 10)
  pmat.CR <- diag(rep(.5, 2 * n)) %*% (dmat.CR + 1) %*% diag(rep(.5, 2 * n))
  pmat.CR.true <- pmat.CR
  pmat.CR[pmat.CR == 0] <- 1
  d.tilde.CR <- dmat.CR + (dmat.CR == -1) + diag(rep(1, 2 * n))
  d.tilde.wt.CR <- d.tilde.CR / pmat.CR

  # ourpmat <- declaration_to_condition_pr_mat(randomizr::declare_ra(N = 400, prob = 0.5, simple = F))
  # pmat.CR.true[1:5, 1:5]
  # ourpmat[1:5, 1:5]
  # pmat.CR.true[401:410, 401:410]
  # ourpmat[401:410, 401:410]
  # pmat.CR.true[401:410, 1:10]
  # ourpmat[401:410,1:10]

  ## DGP with truly random 0.5 chance for each unit being treated
  dat <-
    data.frame(
      p = 0.5,
      z = rbinom(n, 1, 0.5),
      y0 = rnorm(n, sd = 3)
    )

  # Constant treatment effects, SRS
  dat$y1 <- dat$y0 + 3
  Y <- c(-dat$y0, dat$y1)
  R <- c(1 - dat$z, dat$z)
  pi.inv <- (c(rep(1 / (1 - dat$p)), rep(1 / dat$p)))

  ht_est <- sum(R * pi.inv * Y) / n
  y1.hat <- dat$y0 + ht_est
  y0.hat <- dat$y1 - ht_est
  # true_ses_ht <- sqrt(t(Y)%*%dmat%*%Y)/n
  Y.hat <- R * Y + (1 - R) * c(-y0.hat, y1.hat)
  se_ht <- sqrt(t(Y * R) %*% d.tilde.wt %*% (Y * R)) / n
  se_constant_ht <- sqrt(t(Y.hat) %*% dmat %*% Y.hat) / n

  dat$y <- ifelse(dat$z == 1, dat$y1, dat$y0)

  # Simple random assignment
  ht_decl_o <- horvitz_thompson(
    y ~ z,
    data = dat,
    declaration = randomizr::declare_ra(
      N = nrow(dat),
      prob = dat$p[1],
      simple = TRUE
    )
  )

  # Second way to do same estimator, since it's SRS
  ht_prob_o <- horvitz_thompson(
    y ~ z,
    data = dat,
    condition_prs = p
  )

  expect_equal(
    tidy(ht_decl_o)[, c("coefficients", "se")],
    tidy(ht_prob_o)[, c("coefficients", "se")]
  )
  expect_equivalent(
    tidy(ht_decl_o)[, c("coefficients", "se")],
    c(ht_est, se_ht)
  )

  # Now with constant effects assumption
  ht_const_o <- horvitz_thompson(
    y ~ z,
    data = dat,
    declaration = randomizr::declare_ra(
      N = nrow(dat),
      prob = dat$p[1],
      simple = TRUE
    ),
    se_type = "constant"
  )

  expect_equivalent(
    tidy(ht_const_o)[, c("coefficients", "se")],
    c(ht_est, se_constant_ht)
  )

  ## Constant treatment effects, CRS
  dat$z <- sample(rep(0:1, each = n / 2))
  dat$y <- ifelse(dat$z == 1, dat$y1, dat$y0)

  R <- c(1 - dat$z, dat$z)
  pi.inv <- (c(rep(1 / (1 - dat$p)), rep(1 / dat$p)))

  ht_comp_est <- sum(R * pi.inv * Y) / n
  Y.hat <- R * Y + (1 - R) * c(-y0.hat, y1.hat)
  se_comp_ht <- sqrt(t(Y * R) %*% d.tilde.wt.CR %*% (Y * R)) / n
  se_comp_constant_ht <- sqrt(t(Y.hat) %*% dmat.CR %*% Y.hat) / n


  # complete random assignment
  ht_comp_decl_o <- horvitz_thompson(
    y ~ z,
    data = dat,
    declaration = randomizr::declare_ra(
      N = nrow(dat),
      prob = dat$p[1],
      simple = FALSE
    ),
    return_condition_pr_mat = T
  )
  # ht_comp_decl_o$condition_pr_mat[1:5, 1:5]
  # pmat.CR.true[1:5, 1:5]

  # Don't match right now because pmats are diff
  # expect_equal(
  #   tidy(ht_comp_decl_o)[, c("coefficients", "se")],
  #   c(ht_comp_est, se_comp_ht)
  # )

  # Does match if I use JM's pmat
  ht_comp_decl_o <- horvitz_thompson(
    y ~ z,
    data = dat,
    condition_pr_mat = pmat.CR.true
  )
  expect_equivalent(
    tidy(ht_comp_decl_o)[, c("coefficients", "se")],
    c(ht_comp_est, se_comp_ht)
  )

  # Now with constant effects assumption
  # ht_comp_const_o  <- horvitz_thompson(
  #   y ~ z,
  #   data = dat,
  #   declaration = randomizr::declare_ra(
  #     N = nrow(dat),
  #     prob = dat$p[1],
  #     simple = FALSE
  #   ),
  #   se_type = "constant"
  # )
  #
  # expect_equivalent(
  #   tidy(ht_comp_const_o)[, c("coefficients", "se")],
  #   c(ht_comp_est, se_comp_constant_ht)
  # )


  # ht_comp_const_o  <- horvitz_thompson(
  #   y ~ z,
  #   data = dat,
  #   condition_pr_mat = pmat.CR.true,
  #   se_type = "constant"
  # )
  # expect_equivalent(
  #   tidy(ht_comp_const_o)[, c("coefficients", "se")],
  #   c(ht_comp_est, se_comp_constant_ht)
  # )

  # Not matching so we error
  expect_error(
    horvitz_thompson(
      y ~ z,
      data = dat,
      condition_pr_mat = pmat.CR.true,
      se_type = "constant"
    ),
    "`se_type` = 'constant' only supported for simple"
  )
})
