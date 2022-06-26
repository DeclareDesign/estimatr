context("Estimator - iv_robust")

N <- 20
dat <- data.frame(
  y = rnorm(N),
  x = rnorm(N),
  x2 = rnorm(N),
  z2 = rnorm(N),
  w = runif(N),
  clust = sample(letters[1:3], size = N, replace = TRUE)
)
dat$z <- dat$x * 0.5 + rnorm(N)
dat$x1_c <- dat$x

test_that("iv_robust warnings and errors are correct", {
  expect_warning(
    ivro <- iv_robust(mpg ~ hp + cyl | am, data = mtcars, se_type = "HC0"),
    "More regressors than instruments"
  )

  expect_error(
    iv_robust(mpg ~ hp + cyl, data = mtcars),
    "Must specify a `formula` with both regressors and instruments."
  )
})

test_that("iv_robust matches AER + ivpack", {

  skip_if_not_installed("AER")
  skip_if_not_installed("ivpack")
  ivco <- iv_robust(y ~ x | z, data = dat, se_type = "classical")
  ivfit <- AER::ivreg(y ~ x | z, data = dat)
  ivo <- summary(ivfit)

  expect_equivalent(
    as.matrix(tidy(ivco)[, c("estimate", "std.error", "statistic", "p.value")]),
    coef(ivo)[, 1:4]
  )
  # Same as stata if you specify `small` as a stata option
  # which applies the N / N-k finite sample correction

  expect_equivalent(
    ivfit$fitted.values,
    ivco$fitted.values
  )

  # Stata defaults to HC0 as well, but does HC1 with `small`
  ivro <- iv_robust(y ~ x | z, data = dat, se_type = "HC0")
  capture_output(ivpackrob <- ivpack::robust.se(ivfit))

  expect_equivalent(
    as.matrix(tidy(ivro)[, c("estimate", "std.error", "statistic", "p.value")]),
    ivpackrob[, 1:4]
  )

  expect_equivalent(
    ivfit$fitted.values,
    ivro$fitted.values
  )

  # "Stata" clustered SEs are CR0, but they are the same as below with `small`
  ivclusto <- iv_robust(y ~ x | z, data = dat, se_type = "stata", clusters = clust)
  capture_output(ivpackclust <- ivpack::cluster.robust.se(ivfit, dat$clust))

  # Our p-values are bigger (ivpack is be using less conservative DF, we use J - 1 which
  # is what stata uses for clusters w/ `small` and in OLS)
  expect_equivalent(
    as.matrix(tidy(ivclusto)[, c("estimate", "std.error")]),
    ivpackclust[, c(1, 2)]
  )

  expect_equivalent(
    ivfit$fitted.values,
    ivclusto$fitted.values
  )

  # Weighting classical
  ivw <- AER::ivreg(y ~ x | z, data = dat, weights = w)
  ivcw <- iv_robust(y ~ x | z, data = dat, weights = w, se_type = "classical")
  ivregsum <- summary(ivcw)

  expect_equivalent(
    as.matrix(tidy(ivcw)[, c("estimate", "std.error", "statistic", "p.value")]),
    coef(ivregsum)[, 1:4]
  )

  expect_equivalent(
    ivw$fitted.values,
    ivcw$fitted.values
  )

  # HC0 weighted
  ivrw <- iv_robust(y ~ x | z, data = dat, weights = w, se_type = "HC0")
  capture_output(ivpackrobw <- ivpack::robust.se(ivw))

  expect_equivalent(
    as.matrix(tidy(ivrw)[, c("estimate", "std.error", "statistic", "p.value")]),
    ivpackrobw[, 1:4]
  )

  expect_equivalent(
    ivrw$fitted.values,
    ivcw$fitted.values
  )

  # Stata weighted
  ivclrw <- iv_robust(y ~ x | z, data = dat, clusters = clust, weights = w, se_type = "stata")
  ivclw <- AER::ivreg(y ~ x | z, data = dat, weights = w)
  capture_output(ivclwse <- ivpack::cluster.robust.se(ivclw, clusterid = dat$clust))

  expect_equivalent(
    as.matrix(tidy(ivclrw)[1:2, c("estimate", "std.error")]),
    ivclwse[, c(1, 2)]
  )

  expect_equivalent(
    ivclrw$fitted.values,
    ivclw$fitted.values
  )

  # Rank-deficiency
  # HC0
  ivdefr <- iv_robust(y ~ x + x1_c| z + z2, data = dat, se_type = "HC0")
  ivdef <- AER::ivreg(y ~ x + x1_c| z + z2, data = dat)
  capture_output(ivdefse <- ivpack::robust.se(ivdef))

  expect_equal(
    coef(ivdefr),
    coef(ivdef)
  )

  expect_equivalent(
    as.matrix(tidy(ivdefr)[1:2, c("estimate", "std.error", "statistic", "p.value")]),
    ivdefse[, 1:4]
  )

  expect_equivalent(
    ivdefr$fitted.values,
    ivdef$fitted.values
  )

  # # Does not work if instrument is collinear with other instrument
  # ivdefri <- iv_robust(y ~ z + z2| x + x1_c, data = dat, se_type = "HC0")
  # ivdefi <- AER::ivreg(y ~ z + z2| x + x1_c, data = dat)
  # ivdefsei <- ivpack::robust.se(ivdefi)
  #
  # # No longer equal!
  # expect_equal(
  #   coef(ivdefri),
  #   coef(ivdefi)
  # )

  # expect_equivalent(
  #   as.matrix(tidy(ivdefri)[1:2, c("estimate", "std.error", "statistic", "p.value")]),
  #   ivdefsei[, 1:4]
  # )

  # Stata
  ivdefclr <- iv_robust(y ~ x + x1_c | z + z2, data = dat, clusters = clust, se_type = "stata")
  ivdefcl <- AER::ivreg(y ~ x + x1_c | z + z2, data = dat)
  capture_output(ivdefclse <- ivpack::cluster.robust.se(ivdefcl, clusterid = dat$clust))

  expect_equal(
    coef(ivdefclr),
    coef(ivdefcl)
  )

  expect_equivalent(
    as.matrix(tidy(ivdefclr)[1:2, c("estimate", "std.error")]),
    ivdefclse[, c(1, 2)]
  )

  expect_equivalent(
    ivdefclr$fitted.values,
    ivdefcl$fitted.values
  )

  # HC0 Weighted
  ivdefrw <- iv_robust(y ~ x + x1_c| z + z2, weights = w, data = dat, se_type = "HC0")
  ivdefw <- AER::ivreg(y ~ x + x1_c| z + z2, weights = w, data = dat)
  capture_output(ivdefsew <- ivpack::robust.se(ivdefw))

  expect_equal(
    coef(ivdefrw),
    coef(ivdefw)
  )

  expect_equivalent(
    as.matrix(tidy(ivdefrw)[1:2, c("estimate", "std.error", "statistic", "p.value")]),
    ivdefsew[, 1:4]
  )

  expect_equivalent(
    ivdefrw$fitted.values,
    ivdefw$fitted.values
  )

  # Stata weighted
  ivdefclr <- iv_robust(y ~ x + x1_c | z + z2, data = dat, weights = w, clusters = clust, se_type = "stata")
  ivdefcl <- AER::ivreg(y ~ x + x1_c | z + z2, data = dat, weights = w)
  capture_output(ivdefclse <- ivpack::cluster.robust.se(ivdefcl, clusterid = dat$clust))

  expect_equal(
    coef(ivdefclr),
    coef(ivdefcl)
  )

  expect_equivalent(
    as.matrix(tidy(ivdefclr)[1:2, c("estimate", "std.error")]),
    ivdefclse[, c(1, 2)]
  )

  expect_equivalent(
    ivdefclr$fitted.values,
    ivdefcl$fitted.values
  )

  # F-stat fails properly with blocks of size 1
  set.seed(42)
  N <- 20
  dat <- data.frame(y = rnorm(N), x = rnorm(N), z = rnorm(N), bl = sample(letters, size = N, replace = T))
  ivr <- iv_robust(y ~ bl + x | bl + z, data = dat, se_type = "stata")
  expect_equivalent(
    ivr$fstatistic[1],
    NA_integer_
  )

})

test_that("iv_robust matches AER + clubSandwich", {
  skip_if_not_installed("AER")
  skip_if_not_installed("clubSandwich")
  skip_on_cran()

  # ClubSandwich IV tests
  for (se_type in cr_se_types) {

    ivfit <- AER::ivreg(y ~ x | z, data = dat)
    ivfitw <- AER::ivreg(y ~ x | z, weights = w, data = dat)

    # Standard IV models
    ivcr <- iv_robust(y ~ x | z, data = dat, clusters = clust, se_type = se_type)
    clubsand <- clubSandwich::coef_test(ivfit,
                                        vcov = ifelse(se_type == "stata", "CR1S", se_type),
                                        cluster = dat$clust,
                                        test = ifelse(se_type == "CR2", "Satterthwaite", "naive-t"))

    clubsand <- as.data.frame(clubsand)

    cols <- c("estimate", "std.error", "statistic", "df", "p.value")

    expect_equivalent(
      as.matrix(tidy(ivcr)[, cols]),
      as.matrix(clubsand[,-1])
    )

    expect_equivalent(
      ivfit$fitted.values,
      ivcr$fitted.values
    )

    # Weighted IV models
    ivcrw <- iv_robust(y ~ x | z, data = dat, clusters = clust, weights = w, se_type = se_type)
    clubsandw <- clubSandwich::coef_test(ivfitw,
                                         vcov = ifelse(se_type == "stata", "CR1S", se_type),
                                         cluster = dat$clust,
                                         test = ifelse(se_type == "CR2", "Satterthwaite", "naive-t"))

    clubsandw <- as.data.frame(clubsandw)

    expect_equivalent(
      as.matrix(tidy(ivcrw)[, cols]),
      as.matrix(clubsandw[,-1])
    )

    expect_equivalent(
      ivfitw$fitted.values,
      ivcrw$fitted.values
    )

    # Rank-deficiency
    ivfit_rd <- AER::ivreg(y ~ x + x1_c | z + z2, data = dat)
    ivcr_rd <- iv_robust(y ~ x + x1_c | z + z2, data = dat, clusters = clust, se_type = se_type)
    clubsand_rd <- clubSandwich::coef_test(ivfit_rd,
                                           vcov = ifelse(se_type == "stata", "CR1S", se_type),
                                           cluster = dat$clust,
                                           test = ifelse(se_type == "CR2", "Satterthwaite", "naive-t"))

    clubsand_rd <- as.data.frame(clubsand_rd)

    expect_equivalent(
      na.omit(as.matrix(tidy(ivcr_rd)[, cols])),
      na.omit(as.matrix(clubsand_rd[,-1]))
    )

    expect_equivalent(
      ivfit_rd$fitted.values,
      ivcr_rd$fitted.values
    )

    # Rank-deficient, weighted
    ivfitw_rd <- AER::ivreg(y ~ x + x1_c | z + z2, weights = w, data = dat)
    ivcrw_rd <- iv_robust(y ~ x + x1_c | z + z2, data = dat, weights = w, clusters = clust, se_type = se_type)
    clubsandw_rd <- clubSandwich::coef_test(ivfitw_rd,
                                            vcov = ifelse(se_type == "stata", "CR1S", se_type),
                                            cluster = dat$clust,
                                            test = ifelse(se_type == "CR2", "Satterthwaite", "naive-t"))

    clubsandw_rd <- as.data.frame(clubsandw_rd)

    expect_equivalent(
      na.omit(as.matrix(tidy(ivcrw_rd)[, cols])),
      na.omit(as.matrix(clubsandw_rd[,-1]))
    )

    expect_equivalent(
      ivfitw_rd$fitted.values,
      ivcrw_rd$fitted.values
    )
  }

})

test_that("iv_robust different specifications work", {
  skip_if_not_installed("AER")

  # More instruments than endog. regressors
  ivro <- iv_robust(mpg ~ wt | hp + cyl, data = mtcars, se_type = "HC0")
  ivo <- AER::ivreg(mpg ~ wt | hp + cyl, data = mtcars)
  capture_output(ivpo <- ivpack::robust.se(ivo))
  expect_equivalent(
    as.matrix(tidy(ivro)[, c("estimate", "std.error", "statistic", "p.value")]),
    ivpo[, 1:4]
  )

  # . notation for multiple exog, doesnt work!
  # ivro <- iv_robust(mpg ~ wt + hp + vs | . - vs + cyl, data = mtcars, se_type = "HC0")
  # ivo <- AER::ivreg(mpg ~ wt + hp + vs | . - vs + cyl, data = mtcars)
  # ivpo <- ivpack::robust.se(ivo)
  # expect_equivalent(
  #   as.matrix(tidy(ivro)[, c("estimate", "std.error", "statistic", "p.value")]),
  #   ivpo[, 1:4]
  # )

  # . notation in general
  ivro <- iv_robust(mpg ~ .| ., data = mtcars, se_type = "HC0")
  ivo <- AER::ivreg(mpg ~ . | ., data = mtcars)
  capture_output(ivpo <- ivpack::robust.se(ivo))

  expect_equivalent(
    as.matrix(tidy(ivro)[, c("estimate", "std.error", "statistic", "p.value")]),
    ivpo[, 1:4]
  )

})

test_that("S3 methods", {
  skip_if_not_installed("AER")


  ivo <- AER::ivreg(mpg ~ hp + cyl | wt + gear, data = mtcars)
  ivro <- iv_robust(mpg ~ hp + cyl | wt + gear, data = mtcars, se_type = "classical")

  expect_equal(
    vcov(ivro),
    vcov(ivo)
  )

  expect_is(
    tidy(ivro),
    "data.frame"
  )

  expect_equal(
    nrow(tidy(ivro)),
    3
  )

  summary(ivro)

  siv <- capture_output(
    summary(ivro),
    print = TRUE
  )

  expect_true(
    grepl(
      "iv\\_robust\\(formula = mpg \\~ hp \\+ cyl \\| wt \\+ gear, data = mtcars,",
      siv
    )
  )

  expect_true(
    grepl(
      "F\\-statistic\\: 33\\.73 on 2 and 29 DF,  p\\-value\\: 2\\.706e\\-08",
      siv
    )
  )

  capture_output(
    expect_equivalent(
      coef(summary(ivro)),
      print(ivro)
    )
  )

  expect_equivalent(
    ivro$fstatistic,
    summary(ivo)$waldtest[-2]
  )

  expect_equal(
    predict(ivro, newdata = mtcars),
    predict(ivo)
  )

  glance_ivro <- glance(ivro)

  expect_equal(nrow(glance_ivro), 1)

  expect_equal(
    colnames(glance(ivro)),
    c("r.squared", "adj.r.squared", "df.residual", "nobs", "se_type",
      "statistic", "p.value", "statistic.weakinst", "p.value.weakinst",
      "statistic.endogeneity", "p.value.endogeneity", "statistic.overid",
      "p.value.overid")
  )

  # no intercept
  ivo <- AER::ivreg(mpg ~ hp + cyl + 0 | wt + gear, data = mtcars)
  ivro <- iv_robust(mpg ~ hp + cyl + 0 | wt + gear, data = mtcars, se_type = "classical")

  expect_equivalent(
    ivro$fstatistic,
    summary(ivo)$waldtest[-2]
  )

})

test_that("IV diagnostics", {

  # Load stata diagnostics
  stata_diags <- read.table(
    "stata-iv-diagnostics.txt",
    col.names = c("formula", "weights", "options", "diag", "df1", "df2", "val", "pval"),
    sep = ";",
    na = ".",
    stringsAsFactors = FALSE
  )

  formulae <- c(
    "(hp = wt)" = mpg ~ hp | wt,
    # mpg ~ 0 + hp | wt,
    "(hp = wt)0" = mpg ~ 0 + hp | 0 + wt,
    "(hp am = wt gear)" = mpg ~ hp + am | wt + gear,
    # mpg ~ 0 + hp + am | wt + gear,
    "(hp am = wt gear)0" = mpg ~ 0 + hp + am | 0 + wt + gear,
    "gear (hp = wt)" = mpg ~ hp + gear | wt + gear,
    # mpg ~ 0 + hp + gear | wt + gear,
    "gear (hp = wt)0" = mpg ~ 0 + hp + gear | 0 + wt + gear,
    "gear (hp = wt am)" = mpg ~ hp + gear | wt + gear + am,
    # mpg ~ 0 + hp + gear | wt + gear + am,
    # mpg ~ hp + gear | 0 + wt + gear + am,
    "gear (hp = wt am)0" = mpg ~ 0 + hp + gear | 0 + wt + gear + am
  )

  for (f_n in names(formulae)) {
    f <- formulae[[f_n]]
    ivro <- iv_robust(f, data = mtcars, se_type = "classical", diagnostics = TRUE)
    aer_ivro <- summary(AER::ivreg(f, data = mtcars), diagnostics = TRUE)

    # Sargan stat seems to be wrong for AER for this model (-ve critical value)
    if (f_n == "gear (hp = wt am)0") {
      expect_equivalent(
        build_ivreg_diagnostics_mat(ivro)[1:2, ],
        aer_ivro[["diagnostics"]][1:2, ]
      )
    } else {
      expect_equivalent(
        build_ivreg_diagnostics_mat(ivro),
        aer_ivro[["diagnostics"]]
      )
    }

    stata_diag <- subset(
      stata_diags,
      formula == gsub("0", "", f_n) & (grepl("noconstant", options) == grepl("0", f_n))
    )

    expect_equivalent(
      build_ivreg_diagnostics_mat(ivro, stata = TRUE),
      as.matrix(stata_diag[
        grepl("small", stata_diag$options) & nchar(stata_diag$weights) == 0,
        c("df1", "df2", "val", "pval")
      ]),
      tolerance = 1e-6
    )

    # With weights, don't match `overid` test, as we don't report it
    ivrow <- iv_robust(f, data = mtcars, se_type = "classical", weights = drat, diagnostics = TRUE)
    ivrow_diag_mat <- build_ivreg_diagnostics_mat(ivrow, stata = TRUE)
    expect_equivalent(ivrow_diag_mat[nrow(ivrow_diag_mat), ], rep(NA_real_, 4))

    expect_equivalent(
      ivrow_diag_mat[-nrow(ivrow_diag_mat), ],
      as.matrix(stata_diag[
        grepl("small", stata_diag$options) & nchar(stata_diag$weights) > 0 & stata_diag$diag != "overid",
        c("df1", "df2", "val", "pval")
      ]),
      tolerance = 1e-6
    )

    ivro_hc1 <- iv_robust(f, data = mtcars, se_type = "HC1", diagnostics = TRUE)
    ivrow_hc1 <- iv_robust(f, data = mtcars, se_type = "HC1", weights = drat, diagnostics = TRUE)

    expect_equivalent(
      build_ivreg_diagnostics_mat(ivro_hc1, stata = TRUE),
      as.matrix(stata_diag[
        grepl("rob", stata_diag$options) & nchar(stata_diag$weights) == 0,
        c("df1", "df2", "val", "pval")
      ]),
      tolerance = 1e-6
    )

    # Again, no overid test reported with weights
    ivrow_hc1_diag_mat <- build_ivreg_diagnostics_mat(ivrow_hc1, stata = TRUE)
    expect_equivalent(ivrow_hc1_diag_mat[nrow(ivrow_hc1_diag_mat), ], rep(NA_real_, 4))

    expect_equivalent(
      ivrow_hc1_diag_mat[-nrow(ivrow_hc1_diag_mat), ],
      as.matrix(stata_diag[
        grepl("rob", stata_diag$options) & nchar(stata_diag$weights) > 0 & stata_diag$diag != "overid",
        c("df1", "df2", "val", "pval")
      ]),
      tolerance = 1e-6
    )

    # tolerance higher here due to larger values in general
    ivro_crs <- iv_robust(f, data = mtcars, se_type = "stata", clusters = cyl, diagnostics = TRUE)
    ivro_crs_diag_mat <- build_ivreg_diagnostics_mat(ivro_crs, stata = TRUE)

    ivrow_crs <- iv_robust(f, data = mtcars, se_type = "stata", clusters = cyl, weights = drat, diagnostics = TRUE)
    ivrow_crs_diag_mat <- build_ivreg_diagnostics_mat(ivrow_crs, stata = TRUE)
    expect_equivalent(ivrow_crs_diag_mat[nrow(ivrow_crs_diag_mat), ], rep(NA_real_, 4))

    expect_equivalent(
      ivro_crs_diag_mat[-nrow(ivro_crs_diag_mat), ],
      as.matrix(stata_diag[
        grepl("cluster", stata_diag$options) & nchar(stata_diag$weights) == 0 & stata_diag$diag != "overid",
        c("df1", "df2", "val", "pval")
      ]),
      tolerance = 1e-3
    )

    # Stata doesn't report overid test with clusters
    expect_equivalent(
      ivrow_crs_diag_mat[-nrow(ivrow_crs_diag_mat), ],
      as.matrix(stata_diag[
        grepl("cluster", stata_diag$options) & nchar(stata_diag$weights) > 0 & stata_diag$diag != "overid",
        c("df1", "df2", "val", "pval")
        ]),
      tolerance = 1e-3
    )

    # Sanity check unmatched diagnostics
    for (se_type in se_types) {
      ivro <- iv_robust(f, data = mtcars, se_type = se_type, diagnostics = TRUE)
      diagnostic_mat <- build_ivreg_diagnostics_mat(ivro)
      expect_true(
        all(diagnostic_mat[1:2, ] > 0) & all(diagnostic_mat[3, ] >= 0 | is.na(diagnostic_mat[3, ]))
      )
    }

    # Test default se_type
    ivro <- iv_robust(f, data = mtcars, diagnostics = TRUE)
    diagnostic_mat <- build_ivreg_diagnostics_mat(ivro)
    expect_true(
      all(diagnostic_mat[1:2, ] > 0) & all(diagnostic_mat[3, ] >= 0 | is.na(diagnostic_mat[3, ]))
    )

    for (cr_se_type in cr_se_types) {
      ivro <- iv_robust(f, data = mtcars, se_type = cr_se_type, clusters = cyl, diagnostics = TRUE)
      diagnostic_mat <- build_ivreg_diagnostics_mat(ivro)
      expect_true(
        all(diagnostic_mat[1:2, ] > 0) & all(diagnostic_mat[3, ] >= 0 | is.na(diagnostic_mat[3, ]))
      )
    }

  }

})
