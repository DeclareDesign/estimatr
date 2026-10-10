# Records fixest and plm answers for tests/testthat/test_vs_fixest_plm.R.
#
# Run by hand, not part of R CMD check. The values are frozen rather than
# computed live because both packages change their small-sample defaults
# between releases: a live comparison would fail on someone else's release
# note, which is noise about fixest rather than signal about estimatr.
# sandwich and clubSandwich are compared live instead, because they are stable,
# single-purpose, and cheap enough to carry in Suggests.
#
#   Rscript data-raw/make_external_reference.R

library(fixest)
library(plm)
library(blkvar)
# blkvar calls dplyr::n() without importing it.
suppressMessages(library(dplyr))

source("tests/testthat/helper-external.R")

values <- list()

# ---- fixest ----------------------------------------------------------------
#
# `ssc(adj = TRUE, fixef.K = "full")` is the configuration that corresponds to
# estimatr's small-sample correction: absorbed parameters are counted in K,
# including when they are nested inside a cluster. `fixef.K = "nested"`, which
# is fixest's default, does not count nested absorbed parameters and differs
# from estimatr by about 2e-3 on the nested data below.
ssc_full <- ssc(adj = TRUE, fixef.K = "full", cluster.adj = TRUE)

d_fe <- ext_data_fe()

fe_fits <- list(
  fe1_hetero   = feols(y ~ x + z | g,      data = d_fe, vcov = "hetero", ssc = ssc_full),
  fe2_hetero   = feols(y ~ x + z | g + g2, data = d_fe, vcov = "hetero", ssc = ssc_full),
  fe1_cluster  = feols(y ~ x + z | g,      data = d_fe, cluster = ~cl,   ssc = ssc_full),
  fe2_cluster  = feols(y ~ x + z | g + g2, data = d_fe, cluster = ~cl,   ssc = ssc_full),
  fe1_iid      = feols(y ~ x + z | g,      data = d_fe, vcov = "iid",    ssc = ssc_full),
  fe1_w_hetero = feols(y ~ x + z | g,      data = d_fe, weights = ~w,
                       vcov = "hetero", ssc = ssc_full)
)

d_nest <- ext_data_fe_nested()
fe_fits$nested_cluster <- feols(y ~ x + z | g, data = d_nest,
                                cluster = ~cl, ssc = ssc_full)
fe_fits$nested_cluster_K_nested <-
  feols(y ~ x + z | g, data = d_nest, cluster = ~cl,
        ssc = ssc(adj = TRUE, fixef.K = "nested", cluster.adj = TRUE))

for (nm in names(fe_fits)) {
  f <- fe_fits[[nm]]
  values[[paste0("fixest_", nm)]] <- list(
    coefficients = coef(f)[c("x", "z")],
    std.error    = se(f)[c("x", "z")]
  )
}

# ---- fixest, the 2SLS path -------------------------------------------------
#
# The fixed effects are absorbed out of both stages. fixest does that by
# alternating projections and estimatr by its own solver, so what the
# comparison corroborates is the absorption on the 2SLS path, which no other
# external reference in the suite reaches (review C4).
d_iv <- ext_data_fe_iv()

iv_fits <- list(
  iv_fe1_iid = feols(y ~ x | g | en ~ inst, data = d_iv, vcov = "iid",
                          ssc = ssc_full),
  iv_fe1_hetero = feols(y ~ x | g | en ~ inst, data = d_iv, vcov = "hetero",
                          ssc = ssc_full),
  iv_fe1_cluster = feols(y ~ x | g | en ~ inst, data = d_iv, cluster = ~cl,
                          ssc = ssc_full),
  iv_fe2_hetero = feols(y ~ x | g + cl | en ~ inst, data = d_iv,
                          vcov = "hetero", ssc = ssc_full),
  iv_fe1_w_hetero = feols(y ~ x | g | en ~ inst, data = d_iv, weights = ~w,
                          vcov = "hetero", ssc = ssc_full)
)

for (nm in names(iv_fits)) {
  f <- iv_fits[[nm]]
  values[[paste0("fixest_", nm)]] <- list(
    coefficients = coef(f)[c("fit_en", "x")],
    std.error = se(f)[c("fit_en", "x")]
  )
}

# ---- plm -------------------------------------------------------------------
#
# The within estimator with Arellano's cluster-robust variance and no
# small-sample adjustment is CR0 on the absorbed design. plm reaches it through
# panel machinery that shares no code with estimatr.
pd <- pdata.frame(d_nest, index = c("g", "tm"))
pm <- plm(y ~ x + z, data = pd, model = "within")

values$plm_within_arellano_hc0 <- list(
  coefficients = coef(pm)[c("x", "z")],
  std.error    = sqrt(diag(plm::vcovHC(pm, method = "arellano",
                                       type = "HC0", cluster = "group")))[c("x", "z")]
)

# ---- blkvar ----------------------------------------------------------------
#
# blkvar is Pashley and Miratrix's own implementation of the hybrid variance,
# and it is on GitHub only. A live comparison therefore skips wherever blkvar
# is absent, which is every CI platform, so the one check of that variance
# against its authors never ran there (review C1).
#
# The record carries the data as well as the answer. The designs are assigned
# by randomizr::block_ra(), so a record keyed by the seed alone would come to
# describe a different assignment on a randomizr release while still looking
# current, and the comparison would then be against blkvar's answer to another
# question.
for (i in seq_along(ext_blocked_designs())) {
  des <- ext_blocked_designs()[[i]]
  d_bl <- ext_data_blocked(des$block_sizes, des$m_each, seed = i)
  theirs <- blkvar::block_estimator(Yobs = d_bl$y, Z = d_bl$z,
                                    B = factor(d_bl$bl), method = "hybrid_p",
                                    throw.warnings = FALSE)
  values[[paste0("blkvar_hybrid_p_", i)]] <- list(
    data = d_bl,
    coefficients = theirs$ATE_hat,
    std.error = theirs$se_est
  )
}

out <- list(
  values = values,
  versions = c(
    fixest = as.character(packageVersion("fixest")),
    plm = as.character(packageVersion("plm")),
    blkvar = as.character(packageVersion("blkvar")),
    # The assignments the blkvar designs were drawn with are randomizr's.
    randomizr = as.character(packageVersion("randomizr")),
    R = paste0(R.version$major, ".", R.version$minor)
  )
)

saveRDS(out, "tests/testthat/fixtures/external_reference.rds")
message("recorded ", length(values), " external reference values")
message(paste(names(out$versions), out$versions, sep = " ", collapse = " | "))
