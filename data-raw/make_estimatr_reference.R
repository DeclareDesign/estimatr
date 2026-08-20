# Record estimatr 1.0.6's answers so the comparison tests do not have to call it.
#
# Run by hand, from the package root, with estimatr installed:
#
#   Rscript data-raw/make_estimatr_reference.R
#
# Writes tests/testthat/fixtures/estimatr_reference.rds. Nothing in the check
# runs this, and estimatr does not need to be installed to run the tests.
#
# Rerun it only to change what is recorded. A rerun that changes an existing
# value means estimatr's answer moved, which is a finding rather than a
# refresh: read the diff before committing it.

REQUIRED_ESTIMATR <- "1.0.6"

if (!requireNamespace("estimatr", quietly = TRUE)) {
  stop("estimatr is not installed; nothing to record from.")
}
if (!requireNamespace("randomizr", quietly = TRUE)) {
  stop("randomizr is not installed; the Horvitz-Thompson designs need it.")
}
have <- as.character(utils::packageVersion("estimatr"))
if (have != REQUIRED_ESTIMATR) {
  stop("estimatr ", have, " is installed but the reference is defined against ",
       REQUIRED_ESTIMATR, ". Install that version or change REQUIRED_ESTIMATR ",
       "deliberately.")
}

library(estimatr)
source("tests/testthat/helper-data.R")

# Keep the parts of a fit that a test can compare: atomic vectors and matrices,
# and nested fits such as lh_robust's two components. Everything else is either
# unserialisable in a useful way (terms, call) or an input rather than an
# answer (the model frame).
DROP <- c("terms", "call", "model", "outcome", "terms_regressors", "weights")

record <- function(x) {
  if (is.atomic(x)) return(x)
  if (is.data.frame(x)) return(as.list(x))
  if (is.list(x)) {
    x <- x[setdiff(names(x), DROP)]
    keep <- vapply(x, function(el) is.atomic(el) || is.list(el), logical(1))
    return(lapply(x[keep], record))
  }
  NULL
}

values <- list()
put <- function(key, value) {
  if (key %in% names(values)) stop("duplicate reference key: ", key)
  values[[key]] <<- value
  invisible(NULL)
}

# ---- test_vs_estimatr.R ----

dat <- ref_data_vs()

put("lmr_HC0",       record(estimatr::lm_robust(y ~ x + z, data = dat, se_type = "HC0")))
put("lmr_HC1",       record(estimatr::lm_robust(y ~ x + z, data = dat, se_type = "HC1")))
put("lmr_HC2",       record(estimatr::lm_robust(y ~ x + z, data = dat)))
put("lmr_HC3",       record(estimatr::lm_robust(y ~ x + z, data = dat, se_type = "HC3")))
put("lmr_classical", record(estimatr::lm_robust(y ~ x + z, data = dat, se_type = "classical")))
put("lmr_stata",     record(estimatr::lm_robust(y ~ x + z, data = dat, se_type = "stata")))

put("lmr_CR2",      record(estimatr::lm_robust(y ~ x + z, data = dat, clusters = cl)))
put("lmr_CR0",      record(estimatr::lm_robust(y ~ x + z, data = dat, clusters = cl, se_type = "CR0")))
put("lmr_cl_stata", record(estimatr::lm_robust(y ~ x + z, data = dat, clusters = cl, se_type = "stata")))

put("lmr_w",    record(estimatr::lm_robust(y ~ x + z, data = dat, weights = w)))
put("lmr_w_cl", record(estimatr::lm_robust(y ~ x + z, data = dat, weights = w, clusters = cl)))

put("lmr_mv",    record(estimatr::lm_robust(cbind(y, y2) ~ x + z, data = dat)))
put("lmr_noint", record(estimatr::lm_robust(y ~ 0 + x + z, data = dat)))

put("lin",    record(estimatr::lm_lin(y ~ z, covariates = ~ x, data = dat)))
put("lin_cl", record(estimatr::lm_lin(y ~ z, covariates = ~ x, data = dat, clusters = cl)))
put("lin_w",  record(estimatr::lm_lin(y ~ z, covariates = ~ x, data = dat, weights = w)))
put("lin_mv", record(estimatr::lm_lin(y ~ z, covariates = ~ x + y2, data = dat)))

put("iv_HC2",       record(estimatr::iv_robust(y ~ z | iv, data = dat)))
put("iv_classical", record(estimatr::iv_robust(y ~ z | iv, data = dat, se_type = "classical")))
put("iv_CR2",       record(estimatr::iv_robust(y ~ z | iv, data = dat, clusters = cl)))
put("iv_diag",      record(estimatr::iv_robust(y ~ z | iv, data = dat, diagnostics = TRUE)))

put("dim_std", record(estimatr::difference_in_means(y ~ z, data = dat)))
put("dim_cl",  record(estimatr::difference_in_means(y ~ cl_z, clusters = cl, data = dat)))
put("dim_bl",  record(estimatr::difference_in_means(y ~ z_block, blocks = bl, data = dat)))
put("dim_w",   record(estimatr::difference_in_means(y ~ z, weights = w, data = dat)))

dp <- ref_data_pairs()
put("dim_mp", record(estimatr::difference_in_means(y ~ z, blocks = pr, data = dp)))

put("lh", record(estimatr::lh_robust(y ~ x + z, data = dat,
                                     linear_hypothesis = "z + 2*x = 0")))

put("lmr_subset", record(estimatr::lm_robust(y ~ x, data = dat, subset = z == 1)))

dat_f <- dat
dat_f$bl_f <- factor(dat_f$bl)
put("lmr_factor", record(estimatr::lm_robust(y ~ z + bl_f, data = dat_f)))

put("lmr_none", record(estimatr::lm_robust(y ~ x + z, data = dat, se_type = "none")))

# ---- test_fixed_effects.R ----

fe <- ref_data_fe()

put("fe_HC1", record(estimatr::lm_robust(y ~ z + x, data = fe, fixed_effects = ~bl, se_type = "HC1")))
put("fe_HC0", record(estimatr::lm_robust(y ~ z + x, data = fe, fixed_effects = ~bl, se_type = "HC0")))
put("fe_cl_stata", record(estimatr::lm_robust(y ~ z + x, data = fe, fixed_effects = ~bl,
                                              clusters = cl, se_type = "stata")))

# HC2, HC3 and CR2 under `fixed_effects`. 1.0.6 answered all of these by
# expanding the dummies, and so does this version, so they are recorded rather
# than left to a same-session dummy comparison alone.
put("fe_2way_HC2", record(estimatr::lm_robust(y ~ z + x, data = fe,
                                              fixed_effects = ~ bl + cl, se_type = "HC2")))
put("fe_2way_HC3", record(estimatr::lm_robust(y ~ z + x, data = fe,
                                              fixed_effects = ~ bl + cl, se_type = "HC3")))
put("fe_cl_CR2", record(estimatr::lm_robust(y ~ z + x, data = fe, fixed_effects = ~ bl,
                                            clusters = cl, se_type = "CR2")))
put("fe_iv_cl_CR2", record(estimatr::iv_robust(y ~ z | iv, data = fe, fixed_effects = ~ bl,
                                               clusters = cl, se_type = "CR2")))

fe_w <- ref_data_fe_weighted()
put("fe_w_HC1", record(estimatr::lm_robust(y ~ z + x, data = fe_w, fixed_effects = ~bl,
                                           weights = w, se_type = "HC1")))

# ---- test_fe_leverage.R ----

configs <- ref_data_leverage()
for (nm in names(configs)) {
  for (se in c("HC2", "HC3")) {
    fit <- estimatr::lm_robust(y ~ x + z, fixed_effects = ~ g,
                               data = configs[[nm]], se_type = se)
    put(paste0("lev_", nm, "_", se), record(fit))
  }
}
for (se in c("HC2", "HC3")) {
  fit <- estimatr::lm_robust(y ~ x + z, fixed_effects = ~ g,
                             data = configs$unbalanced, se_type = se,
                             weights = wts)
  put(paste0("levw_unbalanced_", se), record(fit))
}

# ---- test_horvitz_thompson.R ----

# The assignment vector is drawn here and recorded alongside the answer, so the
# test reads both rather than redrawing. Redrawing would make the comparison
# depend on how much RNG the rest of the file had consumed first, which is a
# property of test execution order and not of either package.
ht_dat <- ref_data_ht()
N <- nrow(ht_dat)

ht_record <- function(key, decl, dat, seed) {
  set.seed(seed)
  dat$Z <- randomizr::conduct_ra(decl)
  fit <- estimatr::horvitz_thompson(y ~ Z, data = dat, ra_declaration = decl)
  put(key, list(Z = dat$Z,
                coefficients = fit$coefficients,
                std.error = fit$std.error))
}

ht_record("ht_simple", randomizr::declare_ra(N = N, prob = 0.4, simple = TRUE),
          ht_dat, seed = 101)
ht_record("ht_complete", randomizr::declare_ra(N = N, m = 16), ht_dat, seed = 102)

d_bl <- ht_dat
d_bl$bl <- rep(1:4, each = 10)
ht_record("ht_blocked", randomizr::declare_ra(N = N, blocks = d_bl$bl),
          d_bl, seed = 103)

d_bl2 <- ht_dat
d_bl2$bl2 <- rep(1:2, each = 20)
d_bl2$cl2 <- rep(1:10, each = 4)
ht_record("ht_blocked_noninteger",
          randomizr::declare_ra(N = N, blocks = d_bl2$bl2, clusters = d_bl2$cl2),
          d_bl2, seed = 104)

d_cl <- ht_dat
d_cl$cl <- rep(1:8, each = 5)
ht_record("ht_cl_complete", randomizr::declare_ra(N = N, clusters = d_cl$cl),
          d_cl, seed = 105)
ht_record("ht_cl_simple",
          randomizr::declare_ra(N = N, clusters = d_cl$cl, simple = TRUE),
          d_cl, seed = 106)

set.seed(107)
perm <- replicate(20, {
  z <- rep(0L, N)
  z[sample(N, 16)] <- 1L
  z
})
decl_perm <- randomizr::declare_ra(permutation_matrix = perm)
d_perm <- ht_dat
d_perm$Z <- perm[, 1L]
fit_perm <- estimatr::horvitz_thompson(y ~ Z, data = d_perm, ra_declaration = decl_perm)
put("ht_perm", list(permutation_matrix = perm,
                    Z = d_perm$Z,
                    coefficients = fit_perm$coefficients,
                    std.error = fit_perm$std.error))

set.seed(108)
d_named <- data.frame(y = rnorm(20), Z = c(rep(1L, 10L), rep(0L, 10L)))
fit_named <- estimatr::horvitz_thompson(
  y ~ Z, data = d_named,
  condition_prs = ifelse(d_named$Z == 1, 0.5, 0.5)
)
put("ht_named", list(data = d_named,
                     coefficients = fit_named$coefficients,
                     std.error = fit_named$std.error))

# ---- test_postestimation.R ----

post <- ref_data_post()
ht_post <- estimatr::horvitz_thompson(y ~ z, data = post, condition_prs = 0.5)
put("post_ht", list(
  vcov = as.numeric(vcov(ht_post)),
  nobs = nobs(ht_post),
  confint = unname(confint(ht_post)),
  confint_90 = unname(confint(ht_post, level = 0.90)),
  glance_names = names(generics::glance(ht_post))
))

lin_dat <- ref_data_lin()
for (case in ref_lin_cases()) {
  f <- reformulate(case$z, "y")
  fit <- estimatr::lm_lin(f, covariates = case$cov, data = lin_dat)
  put(paste0("post_predict_", case$lbl),
      unname(estimatr:::predict.lm_robust(fit, newdata = lin_dat)))
}

fit_zf <- estimatr::lm_lin(y ~ zf, covariates = ~ x, data = lin_dat)
sub <- lin_dat[lin_dat$zf != "c", ]
put("post_predict_subset",
    unname(estimatr:::predict.lm_robust(fit_zf, newdata = sub)))

nd <- lin_dat
nd$zf <- factor(rep("a", nrow(nd)), levels = "a")
put("post_predict_one_level",
    unname(estimatr:::predict.lm_robust(fit_zf, newdata = nd)))

# ---- the return-object surface ----

# Six fields went missing from 2.0's fitted objects without anything reporting
# it, and one that survived (`proj_fstatistic`) came back wrong. None of it was
# visible to a NAMESPACE diff, to the estimator comparisons, or to a grep of
# the reverse dependencies; it surfaced only because a revdep read one of the
# missing fields. Recording the whole surface is what turns that class of
# regression into a test failure.

surface_d <- ref_surface_data()
for (nm in names(ref_surface_fits(surface_d))) {
  fit <- ref_surface_fits(surface_d)[[nm]]
  put(paste0("surface_", nm), record(fit))
}

# ---- write ----

fixture <- list(
  estimatr_version = have,
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  values = values
)

dir.create("tests/testthat/fixtures", showWarnings = FALSE, recursive = TRUE)
saveRDS(fixture, "tests/testthat/fixtures/estimatr_reference.rds", version = 2)

message("recorded ", length(values), " reference values from estimatr ", have)
