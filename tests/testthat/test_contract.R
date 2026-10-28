# The post-estimation contract ----
#
# inst/contract/contract.csv records, for each probe of a stats or sandwich
# generic on each of ten fits, what estimatr promises: that the call agrees
# with the same call on lm() or ivreg() (or a stated reference), refuses with
# estimatr's own message, errors rather than return something, returns a
# value there is nothing to hold to, or is a known gap. `none` makes no
# promise. A `known_wrong` cell that starts agreeing fails here, so that a fix
# is entered in the contract rather than passing unnoticed.
#
# Each probe runs in its own fresh process holding estimatr and whatever
# namespace that probe reaches, and the references in one that never loads
# estimatr: clubSandwich replaces estimatr's bread() and model.frame() when it
# is loaded second, and performance registers a logLik() for iv_robust. The
# packages are checked with system.file() because skip_if_not_installed()
# loads the namespace into this session.

test_that("every post-estimation promise in the contract holds", {
  skip_on_cran()
  for (pkg in c("sandwich", "clubSandwich", "ivreg")) {
    skip_if(!nzchar(system.file(package = pkg)), paste(pkg, "is not installed"))
  }
  dir <- system.file("contract", package = "estimatr")
  # Under load_all() the children load the same source tree, which is the
  # one place inst/contract sits two levels below a DESCRIPTION
  dev <- file.exists(file.path(dirname(dirname(dir)), "DESCRIPTION"))
  env <- new.env()
  sys.source(file.path(dir, "fixtures.R"), envir = env)
  sys.source(file.path(dir, "observe.R"), envir = env)
  obs <- env$contract_observe(dir, source_dir = if (dev) dirname(dirname(dir)))

  for (i in seq_len(nrow(obs))) {
    row <- obs[i, ]
    lab <- paste0(row$probe, " on ", row$fit, ": promised ", row$promise,
                  ", observed ", row$observed,
                  if (!is.na(row$message)) paste0(" (", row$message, ")"))
    switch(row$promise,
      agrees = ,
      returns = expect_identical(row$observed, row$promise, info = lab),
      errors = expect_identical(row$observed, "errors", info = lab),
      refuses = expect_true(row$observed == "errors" &&
                              grepl(row$error_pattern, row$message, fixed = TRUE),
                            info = lab),
      known_wrong = expect_false(row$observed == "agrees", info = lab),
      none = NULL,
      stop("unknown promise ", row$promise, " in contract row ", i)
    )
  }
})
