# How many assertions compare estimatr against an implementation built
# independently of it?
#
# NEWS.md, README.md, cran-comments.md and vignette("estimatr2.0") all state
# that number and its six components. They were once wrong, because the count
# had been taken as the file total of `test_blocked_variance.R` and
# `test_lm_lin_equivalence.R`, and most of the assertions in those two files
# are internal-consistency checks rather than comparisons. This script is the
# measurement those documents quote, so the number can be re-derived rather
# than remembered.
#
# The rule: an assertion counts when one side of it is an estimatr result and
# the other is a value produced by something with no shared lineage with this
# package. Properties of estimatr checked against itself do not count, nor do
# algebraic identities that never touch an estimatr result, nor do fixture
# sanity checks.
#
# Run with:  Rscript data-raw/count_assertions.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))

# Per test_that block: "all" when every assertion in it is a comparison, or a
# count with the reason the rest are excluded. A block whose name is missing
# from this list, or a name here matching no block, is an error rather than a
# silent zero, so renaming or adding a test forces the count to be revisited.
components <- list(
  sandwich = list(
    file = "test_vs_sandwich.R",
    blocks = list(
      "lm_robust HC0-HC3 match sandwich::vcovHC" = "all",
      "lm_robust classical matches stats::vcov" = "all",
      "weighted lm_robust HC0-HC3 match sandwich::vcovHC" = "all",
      "weighted lm_robust classical matches stats::vcov" = "all",
      "CR0 matches sandwich::vcovCL with no cluster adjustment" = "all",
      "se_type 'stata' matches sandwich::vcovCL with the Stata corrections" = "all",
      "weighted cluster-robust matches sandwich::vcovCL" = "all",
      # Both sides are estimatr: CR0 rescaled by the Stata correction factor.
      "the CR0 and Stata cluster corrections differ as expected" = 0,
      "iv_robust HC0, HC1 and CR0 match sandwich on AER::ivreg" = "all",
      # The third assertion compares ivreg's hatvalues to a hand-built
      # projection; no estimatr value appears on either side of it.
      "iv_robust HC2 and HC3 match sandwich on an ivreg::ivreg fit" = 2,
      # Two of eight compare iv_robust's vcov against a hand-built sandwich
      # filling; the rest are properties of the two leverage conventions.
      "iv_robust uses second-stage leverage, not the influence diagonal" = 2,
      # Algebra on hand-built matrices; estimatr is not involved.
      "the 2SLS hat matrix is idempotent but not symmetric" = 0,
      "second-stage leverage stays in [0, 1] where influence leverage does not" = 0,
      "the two conventions differ by 8.6% at HC2 and 18.5% at HC3 on mtcars" = 0
    )
  ),
  clubSandwich = list(
    file = "test_vs_clubsandwich.R",
    blocks = list(
      "CR2 matches clubSandwich::vcovCR, balanced and unbalanced clusters" = "all",
      "weighted CR2 matches clubSandwich::vcovCR" = "all",
      "CR0 matches clubSandwich::vcovCR type CR0" = "all",
      "CR2 degrees of freedom match clubSandwich Satterthwaite" = "all",
      # Properties of estimatr's own degrees of freedom.
      "CR2 degrees of freedom differ across coefficients" = 0,
      "CR2 with absorbed fixed effects matches clubSandwich on the dummy expansion" = "all",
      "iv_robust CR2 matches clubSandwich on AER::ivreg" = "all"
    )
  ),
  Stata = list(
    file = "test_vs_stata.R",
    blocks = list(
      # One assertion checks that the fixture holds the models the file expects.
      "lm_robust reproduces Stata's reg" = 32,
      "weighted HC2 and HC3 differ from Stata by the documented amount" = "all",
      "the weighted-leverage divergence is confined to HC2 and HC3" = "all",
      # One fixture check, and one assertion that areg reports no F statistic.
      "lm_robust with fixed_effects reproduces Stata's areg" = 6,
      "iv_robust reproduces Stata's ivregress 2sls" = 39,
      "weighted 2SLS root MSE follows AER::ivreg rather than Stata" = "all"
    )
  ),
  fixest_plm = list(
    file = "test_vs_fixest_plm.R",
    blocks = list(
      "one-way absorption matches fixest" = "all",
      "two-way absorption matches fixest" = "all",
      "absorbed parameters are counted in the cluster correction" = "all",
      "the within estimator matches plm with Arellano's variance" = "all",
      # Metadata on the recorded fixture, not a comparison.
      "the recorded reference names the versions it came from" = 0
    )
  ),
  blkvar = list(
    file = "test_blocked_variance.R",
    blocks = list(
      "a block with a singleton arm is estimable" = 0,
      "mixed block sizes no longer collapse to matched pairs" = 0,
      "small blocks of varying size are estimable" = 0,
      "standard errors match blkvar, the authors' own package" = "all",
      # These two check Pashley and Miratrix's equations 4 and 5 written out by
      # hand. They are comparisons, but not against blkvar, so they are left
      # out rather than folded into a component named after that package.
      "all-big designs keep the plain blocked estimator" = 0,
      "matched pairs keep the matched-pairs estimator and df" = 0,
      "a single small block is refused rather than given zero variance" = 0,
      "a small block holding half the small-block sample is refused" = 0,
      "blocks with no variation in treatment still error" = 0,
      "hybrid variance is conservative and covers at the nominal rate" = 0,
      "a block with one treated cluster is refused, not silently wrong" = 0,
      "matched-pair clustered designs are exempt from that check" = 0,
      "clustered blocks with two or more clusters per arm still work" = 0,
      "within-block cluster variance is exact when both arms have 2+ clusters" = 0
    )
  ),
  lin_by_hand = list(
    file = "test_lm_lin_equivalence.R",
    blocks = list(
      "lm_lin equals the hand-built Lin specification across the grid" = "all",
      "every se_type agrees with the hand-built specification when weighted" = "all",
      # Internal consistency of lm_lin's own centring, not the Lin
      # specification: two lm_lin fits, or a weighted mean computed inline.
      "scaling every weight by a constant changes nothing" = 0,
      "a zero weight is the same as removing the row" = 0,
      "the centring point is the weighted mean over complete cases only" = 0,
      "a missing covariate drops the row before centring" = 0,
      "weighted and unweighted centring differ when the weights are informative" = 0
    )
  )
)

files <- vapply(components, function(x) x$file, character(1))
# One run of the whole suite. Running it twice in the same process returns
# 2,237 rather than 2,236, from one assertion in `test_lm_robust.R` that is
# stable when that file is run on its own, so the number depends on session
# history and only a fresh process reproduces what `R CMD check` reports.
suite <- as.data.frame(testthat::test_local(".", reporter = "silent", stop_on_failure = FALSE))
observed <- suite[c("file", "test", "nb")]

count_component <- function(spec) {
  in_file <- observed[observed$file == spec$file, ]
  declared <- names(spec$blocks)
  missing <- setdiff(in_file$test, declared)
  unknown <- setdiff(declared, in_file$test)
  if (length(missing) || length(unknown)) {
    stop(spec$file, " has changed since this count was written.\n",
         if (length(missing)) paste0("  not classified: ", paste(missing, collapse = "; "), "\n"),
         if (length(unknown)) paste0("  no longer present: ", paste(unknown, collapse = "; "), "\n"),
         call. = FALSE)
  }
  n <- vapply(declared, function(b) {
    rule <- spec$blocks[[b]]
    total <- in_file$nb[in_file$test == b]
    if (identical(rule, "all")) total else as.integer(rule)
  }, numeric(1))
  c(external = sum(n), file_total = sum(in_file$nb))
}

counts <- vapply(components, count_component, numeric(2))

cat("\nAssertions comparing estimatr against an independent implementation\n\n")
print(data.frame(
  component = colnames(counts),
  file = unname(files),
  external = counts["external", ],
  file_total = counts["file_total", ],
  row.names = NULL
))
cat("\nexternal total: ", sum(counts["external", ]),
    "   (file totals sum to ", sum(counts["file_total", ]), ")\n", sep = "")

cat("suite total:    ", sum(suite$nb), "\n", sep = "")
cat("recorded 1.0.6: ", sum(suite$nb[suite$file == "test_vs_estimatr.R"]), "\n", sep = "")
