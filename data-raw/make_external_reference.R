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

out <- list(
  values = values,
  versions = c(
    fixest = as.character(packageVersion("fixest")),
    plm    = as.character(packageVersion("plm")),
    R      = paste0(R.version$major, ".", R.version$minor)
  )
)

saveRDS(out, "tests/testthat/fixtures/external_reference.rds")
message("recorded ", length(values), " external reference values")
message(paste(names(out$versions), out$versions, sep = " ", collapse = " | "))
