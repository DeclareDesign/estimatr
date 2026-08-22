# Read the log produced by run-stata-iv-hc2.do and answer the question it was
# sent to answer.
#
#   Rscript data-raw/stata/parse-stata-iv-hc2-log.R path/to/stata-iv-hc2.log
#
# Everything the do-file records is on lines of the form RESULT|key|value, so
# the rest of the log (Stata's own echo of each command) is ignored.

args <- commandArgs(trailingOnly = TRUE)
log_path <- if (length(args)) args[1] else "stata-iv-hc2.log"
stopifnot(file.exists(log_path))

lines <- readLines(log_path, warn = FALSE)
hits <- grep("^RESULT\\|", lines, value = TRUE)
parts <- strsplit(sub("^RESULT\\|", "", hits), "|", fixed = TRUE)
vals <- setNames(trimws(vapply(parts, function(p) paste(p[-1], collapse = "|"), "")),
                 vapply(parts, `[`, "", 1))

get <- function(key, numeric = TRUE) {
  if (!key %in% names(vals)) return(NULL)
  v <- vals[[key]]
  if (!numeric) return(v)
  suppressWarnings(as.numeric(v))
}
mat <- function(prefix, k) {
  out <- matrix(NA_real_, k, k)
  for (i in seq_len(k)) for (j in seq_len(k)) {
    v <- get(sprintf("%s.V_%d_%d", prefix, i, j))
    if (!is.null(v)) out[i, j] <- v
  }
  out
}
rel <- function(a, b) {
  if (is.null(a) || is.null(b) || anyNA(a) || anyNA(b)) return(NA_real_)
  max(abs(a - b) / pmax(abs(b), 1e-300))
}
say <- function(...) cat(sprintf(...), "\n", sep = "")
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

say("Stata %s (%s, %s), run %s",
    get("meta.stata_version", FALSE) %||% "?", get("meta.flavor", FALSE) %||% "?",
    get("meta.os", FALSE) %||% "?", get("meta.date", FALSE) %||% "?")

# ---- did the data match? ----
say("\n== data ==")
d <- read.csv(file.path(dirname(log_path), "mtcars.csv"))
if (!file.exists(file.path(dirname(log_path), "mtcars.csv"))) {
  d <- transform(mtcars, w = drat / 5)
}
ok <- TRUE
for (v in c("mpg", "hp", "am", "wt", "gear", "cyl", "drat", "carb")) {
  s <- get(paste0("data.sum.", v))
  if (is.null(s)) next
  if (!isTRUE(all.equal(s, sum(d[[v]]), tolerance = 1e-12))) {
    say("  MISMATCH on %s: Stata %.10g vs ours %.10g", v, s, sum(d[[v]])); ok <- FALSE
  }
}
say("  N = %s; column sums %s", get("data.N", FALSE) %||% "?",
    if (ok) "all match ours" else "DO NOT match")

# ---- the question ----
say("\n== does this Stata implement HC2/HC3 for ivregress? ==")
for (tag in c("iv_just_hc2", "iv_just_hc3", "iv_over_hc2", "iv_over_hc3",
              "iv_just_hc2_w", "weak_hc2", "weak_hc3")) {
  rc <- get(paste0(tag, ".rc"))
  if (is.null(rc)) { say("  %-14s (not in log)", tag); next }
  say("  %-14s %s", tag,
      if (rc == 0) "ACCEPTED" else sprintf("refused (return code %d)", rc))
}

# ---- if accepted, which convention? ----
accepted <- names(vals)[grepl("\\.rc$", names(vals))]
accepted <- sub("\\.rc$", "", accepted)
accepted <- accepted[vapply(accepted, function(t) isTRUE(get(paste0(t, ".rc")) == 0), TRUE)]
# Only the 2SLS fits: the two conventions coincide for OLS, so an OLS row
# says nothing about which one Stata chose.
hc_fits <- grep("^(iv|weak)_.*hc[23]", accepted, value = TRUE)

if (length(hc_fits)) {
  say("\n== which convention did Stata pick? ==")
  for (tag in hc_fits) {
    dataset <- if (grepl("^weak", tag)) "weak" else "mtcars"
    k <- get(paste0(tag, ".k"))
    if (is.null(k)) next
    which_hc <- if (grepl("hc3", tag)) "hc3" else "hc2"
    stata_v <- mat(tag, k)
    a <- mat(sprintf("%s.mata.%s_a", dataset, which_hc), k)
    b <- mat(sprintf("%s.mata.%s_b", dataset, which_hc), k)
    ra <- rel(stata_v, a); rb <- rel(stata_v, b)
    verdict <- if (is.na(ra) || is.na(rb)) "inconclusive (missing values)" else
      if (ra < 1e-8 && rb >= 1e-8) "matches estimatr's second-stage leverage" else
      if (rb < 1e-8 && ra >= 1e-8) "matches sandwich's influence leverage" else
      if (ra < rb) "closer to estimatr" else "closer to sandwich"
    say("  %-14s vs estimatr %9.2e | vs sandwich %9.2e -> %s", tag, ra, rb, verdict)
  }
} else {
  say("\n== which convention did Stata pick? ==")
  say("  Not applicable: this Stata refused every HC2/HC3 request for ivregress.")
}

# ---- anchors we already hold ----
say("\n== anchors (these should already agree with what we have) ==")
if (requireNamespace("estimatr", quietly = TRUE)) {
  dd <- transform(mtcars, w = drat / 5)
  ord <- c("hp", "am", "(Intercept)")   # Stata reports _cons last
  anchors <- list(
    iv_just_unadj = estimatr::iv_robust(mpg ~ hp + am | wt + gear, data = dd,
                                        se_type = "classical"),
    iv_just_rob   = estimatr::iv_robust(mpg ~ hp + am | wt + gear, data = dd,
                                        se_type = "HC1"),
    ols_hc2       = estimatr::lm_robust(mpg ~ hp, data = dd, se_type = "HC2"),
    ols_hc3       = estimatr::lm_robust(mpg ~ hp, data = dd, se_type = "HC3")
  )
  for (tag in names(anchors)) {
    k <- get(paste0(tag, ".k"))
    if (is.null(k)) next
    o <- if (grepl("^ols", tag)) c("hp", "(Intercept)") else ord
    r <- rel(mat(tag, k), anchors[[tag]]$vcov[o, o])
    say("  %-14s relative difference %9.2e %s", tag, r,
        if (!is.na(r) && r < 1e-4) "ok" else "CHECK THIS")
  }
}

# ---- the jackknife ----
say("\n== jackknife ==")
for (ds in c("mtcars", "weak")) {
  tag <- if (ds == "mtcars") "iv_just_jack" else "weak_jack"
  k <- get(paste0(tag, ".k"))
  if (is.null(k) || !isTRUE(get(paste0(tag, ".rc")) == 0)) next
  sj <- mat(tag, k); mj <- mat(paste0(ds, ".mata.jack"), k)
  say("  %-7s Stata vce(jackknife) vs our Mata jackknife: %9.2e", ds, rel(sj, mj))
  for (conv in c("a", "b")) {
    v3 <- mat(sprintf("%s.mata.hc3_%s", ds, conv), k)
    nm <- if (conv == "a") "estimatr HC3" else "sandwich HC3"
    if (!anyNA(v3) && !anyNA(sj)) {
      say("    %-13s vs jackknife SEs: %+6.1f%% (mean over coefficients)",
          nm, 100 * (mean(sqrt(diag(v3)) / sqrt(diag(sj))) - 1))
    } else {
      say("    %-13s not computable on this data", nm)
    }
  }
}

# ---- the discriminating fact ----
say("\n== leverage bounds ==")
for (ds in c("mtcars", "weak")) {
  ha <- get(paste0(ds, ".mata.max_ha")); hb <- get(paste0(ds, ".mata.max_hb"))
  if (is.null(ha)) next
  say("  %-7s max second-stage leverage %8.4f | max influence leverage %10.4f%s",
      ds, ha, hb, if (!is.null(hb) && hb > 1) "   <- exceeds 1, sandwich HC2 undefined" else "")
}
r2 <- get("weak.first_stage_r2")
if (!is.null(r2)) say("  weak-iv first-stage R-squared: %.4f", r2)
say("")
