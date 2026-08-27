# The source of the three speed tables in vignettes/estimatr2.0.Rmd, and of the
# two-way figures in NEWS.md. Rbuildignored.
#
# Run: Rscript data-raw/benchmark_vignette.R
#
# Two things about this script are load-bearing and were learned the hard way.
#
# First, every row runs in its own R process. The 1.0.6 cells that build a large
# dummy or joint-probability matrix are bound by garbage collection rather than
# by arithmetic, so a cell's time depends on how much the heap had already grown
# when it ran: Horvitz-Thompson at N = 1000 reads 9.6 ms measured after the
# fixed-effects rows and 53.1 ms in a fresh process. The second is what a user
# meets. A benchmark whose value depends on which benchmark ran before it is not
# a measurement.
#
# Second, both versions read one saved data set rather than each generating its
# own from a seed, so a change in either package's RNG use cannot move the data.

# Each worker is a separate process with its own tempdir(), so the shared
# directory travels to them in the environment rather than being recomputed.
S <- Sys.getenv("ESTIMATR_BENCH_DIR")
if (!nzchar(S)) {
  S <- file.path(tempdir(), "estimatr_bench")
  Sys.setenv(ESTIMATR_BENCH_DIR = S)
}
# The 1.0.6 build is expensive and identical every run, so it is cached outside
# tempdir() and reused.
LIB106 <- file.path(tools::R_user_dir("estimatr", "cache"), "lib106")
DATA   <- file.path(S, "bench_data.rds")
dir.create(S, showWarnings = FALSE, recursive = TRUE)

# Row definitions ----
# label, a thunk, and replications. Replications fall where one call is slow.
build_rows <- function(d, v2) {
  ht <- function(nm) {
    dd <- d$ht[[nm]]
    decl <- randomizr::declare_ra(N = nrow(dd))
    # 1.0.6 takes the declaration as `ra_declaration`; 2.0 folded all five
    # design arguments into `condition_prs`.
    arg <- if (v2) list(condition_prs = decl) else list(ra_declaration = decl)
    function() do.call(horvitz_thompson, c(list(formula = y ~ z, data = dd), arg))
  }
  list(
    hc2_1000   = list("lm_robust, HC2, n = 1000",       function() lm_robust(y ~ z + x, data = d$d1000), 200),
    cls_1000   = list("lm_robust, classical, n = 1000", function() lm_robust(y ~ z + x, data = d$d1000, se_type = "classical"), 200),
    cr2_100    = list("lm_robust, CR2, 100 clusters",   function() lm_robust(y ~ z + x, data = d$dcl, clusters = cl), 200),
    cr2_w      = list("lm_robust, CR2, weighted",       function() lm_robust(y ~ z + x, data = d$dcl, clusters = cl, weights = w), 200),
    cr2_5000   = list("lm_robust, CR2, n = 5000, 500 clusters", function() lm_robust(y ~ z + x, data = d$dcl5, clusters = cl), 100),
    hc2_100    = list("lm_robust, HC2, n = 100",        function() lm_robust(y ~ z + x, data = d$d100), 200),
    lin        = list("lm_lin",                          function() lm_lin(y ~ z, covariates = ~ x, data = d$d1000), 200),
    iv         = list("iv_robust",                       function() iv_robust(y ~ z + x2 | z + inst, data = d$iv), 200),
    dim        = list("difference_in_means",             function() difference_in_means(y ~ z, data = d$d1000), 200),
    fe50_hc2   = list("50 blocks, n = 1000, HC2",        function() lm_robust(y ~ z + x, data = d$b50,   fixed_effects = ~ bl), 200),
    fe500_hc2  = list("500 blocks, n = 10,000, HC2",     function() lm_robust(y ~ z + x, data = d$b500,  fixed_effects = ~ bl), if (v2) 100 else 30),
    fe2000_hc2 = list("2000 blocks, n = 40,000, HC2",    function() lm_robust(y ~ z + x, data = d$b2000, fixed_effects = ~ bl), if (v2) 100 else 5),
    fe50_hc1   = list("50 blocks, HC1",                  function() lm_robust(y ~ z + x, data = d$b50,   fixed_effects = ~ bl, se_type = "HC1"), 200),
    fe500_hc1  = list("500 blocks, HC1",                 function() lm_robust(y ~ z + x, data = d$b500,  fixed_effects = ~ bl, se_type = "HC1"), if (v2) 100 else 60),
    fe2000_hc1 = list("2000 blocks, HC1",                function() lm_robust(y ~ z + x, data = d$b2000, fixed_effects = ~ bl, se_type = "HC1"), if (v2) 100 else 30),
    ht200      = list("complete, N = 200",               ht("200"),  if (v2) 200 else 60),
    ht1000     = list("complete, N = 1000",              ht("1000"), if (v2) 200 else 60),
    ht3000     = list("complete, N = 3000",              ht("3000"), if (v2) 200 else 20)
  )
}
ROW_KEYS <- c("hc2_1000", "cls_1000", "cr2_100", "cr2_w", "cr2_5000", "hc2_100",
              "lin", "iv", "dim", "fe50_hc2", "fe500_hc2", "fe2000_hc2",
              "fe50_hc1", "fe500_hc1", "fe2000_hc1", "ht200", "ht1000", "ht3000")

# Worker: one row, one process ----
args <- commandArgs(TRUE)
if (length(args) == 2L) {
  libp <- args[1]; key <- args[2]
  if (nzchar(libp)) .libPaths(c(libp, .libPaths()))
  suppressMessages({library(estimatr); library(microbenchmark); library(randomizr)})
  ver <- as.character(packageVersion("estimatr"))
  row <- build_rows(readRDS(DATA), ver >= "2")[[key]]
  # microbenchmark resolves its expression in its own calling frame, so the
  # thunk has to be a local of that frame rather than a global.
  time_it <- function(f, reps) microbenchmark(f(), times = reps)
  fit <- row[[2]]()
  t <- time_it(row[[2]], row[[3]])$time / 1e6
  saveRDS(list(label = row[[1]], med = median(t),
               est = tryCatch(c(coef(fit)), error = function(e) NA_real_),
               se  = tryCatch(c(fit$std.error), error = function(e) NA_real_)),
          file.path(S, sprintf("row_%s_%s.rds", if (ver >= "2") "200" else "106", key)))
  quit(save = "no")
}

# Orchestrator ----
set.seed(20260827)
mk <- function(n, ncl = NULL, nbl = NULL) {
  d <- data.frame(z = rbinom(n, 1, 0.5), x = rnorm(n))
  d$y <- 0.3 * d$z + 0.5 * d$x + rnorm(n)
  if (!is.null(ncl)) d$cl <- rep(seq_len(ncl), length.out = n)
  if (!is.null(nbl)) d$bl <- rep(seq_len(nbl), length.out = n)
  d
}
iv <- mk(1000)
iv$inst <- rnorm(1000)
iv$x2   <- 0.7 * iv$inst + rnorm(1000)
iv$y    <- 0.3 * iv$z + 0.4 * iv$x2 + rnorm(1000)
dat <- list(d1000 = mk(1000), d100 = mk(100),
            dcl   = within(mk(1000, ncl = 100), w <- runif(1000, 0.5, 2)),
            dcl5  = mk(5000, ncl = 500), b50 = mk(1000, nbl = 50),
            b500  = mk(10000, nbl = 500), b2000 = mk(40000, nbl = 2000), iv = iv)
dat$ht <- lapply(stats::setNames(c(200, 1000, 3000), c("200", "1000", "3000")),
                 function(n) data.frame(y = rnorm(n), z = rep(c(1, 0), each = n / 2)))
saveRDS(dat, DATA)

if (!dir.exists(file.path(LIB106, "estimatr"))) {
  dir.create(LIB106, showWarnings = FALSE, recursive = TRUE)
  message("installing estimatr 1.0.6 into ", LIB106, " (CRAN source)")
  install.packages("estimatr", lib = LIB106, repos = "https://cloud.r-project.org",
                   type = "source", INSTALL_opts = "--no-docs")
}

me <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
for (lib in c("", LIB106)) {
  for (k in ROW_KEYS) {
    message(if (nzchar(lib)) "1.0.6 " else "2.0.0 ", k)
    system2(file.path(R.home("bin"), "Rscript"), c(shQuote(me), shQuote(lib), shQuote(k)),
            stdout = NULL, stderr = NULL)
  }
}

# Tables ----
get <- function(v, k) readRDS(file.path(S, sprintf("row_%s_%s.rds", v, k)))
fmt <- function(x) if (x >= 1000) formatC(round(x), big.mark = ",", format = "d") else sprintf("%.2f", x)
spd <- function(a, b) {
  r <- a / b
  if (r >= 100) paste0(formatC(round(r), big.mark = ",", format = "d"), "x") else sprintf("%.1fx", r)
}
main <- ROW_KEYS[1:9]
cat("\n| task | 1.0.6 | 2.0.0 | speedup |\n|---|---:|---:|---:|\n")
for (k in main) {
  a <- get("106", k); b <- get("200", k)
  cat(sprintf("| `%s` | %s ms | %s ms | %s |\n", a$label, fmt(a$med), fmt(b$med), spd(a$med, b$med)))
}
cat("\n| design | 1.0.6, default (HC2) | 2.0.0, default (HC2) | speedup | 1.0.6, HC1 | 2.0.0, HC1 | speedup |\n|---|---:|---:|---:|---:|---:|---:|\n")
for (r in list(c("fe50_hc2", "fe50_hc1", "50 blocks, n = 1000"),
               c("fe500_hc2", "fe500_hc1", "500 blocks, n = 10,000"),
               c("fe2000_hc2", "fe2000_hc1", "2000 blocks, n = 40,000"))) {
  a2 <- get("106", r[1]); b2 <- get("200", r[1])
  a1 <- get("106", r[2]); b1 <- get("200", r[2])
  cat(sprintf("| %s | %s ms | %s ms | %s | %s ms | %s ms | %s |\n", r[3],
              fmt(a2$med), fmt(b2$med), spd(a2$med, b2$med),
              fmt(a1$med), fmt(b1$med), spd(a1$med, b1$med)))
}
cat("\n| design | 1.0.6 | 2.0.0 | speedup |\n|---|---:|---:|---:|\n")
for (k in c("ht200", "ht1000", "ht3000")) {
  a <- get("106", k); b <- get("200", k)
  cat(sprintf("| %s | %s ms | %s ms | %s |\n", a$label, fmt(a$med), fmt(b$med), spd(a$med, b$med)))
}

# The vignette's agreement claim, which is a promise and not a description.
gaps  <- vapply(ROW_KEYS, function(k) {
  a <- get("106", k); b <- get("200", k)
  max(c(abs(a$est - b$est), abs(a$se - b$se)), na.rm = TRUE)
}, 0)
ident <- vapply(ROW_KEYS, function(k) {
  a <- get("106", k); b <- get("200", k)
  identical(a$est, b$est) && identical(a$se, b$se)
}, TRUE)
gains <- vapply(ROW_KEYS, function(k) get("106", k)$med / get("200", k)$med, 0)
cat(sprintf("\nbit-identical rows: %d of %d\nworst gap: %.1e (%s)\nsmallest gain: %.1fx (%s)\n",
            sum(ident), length(ROW_KEYS), max(gaps), ROW_KEYS[which.max(gaps)],
            min(gains), ROW_KEYS[which.min(gains)]))
stopifnot(max(gaps) < 1e-10, min(gains) > 1)
cat("\nthe vignette's agreement and no-row-slower claims hold\n")
