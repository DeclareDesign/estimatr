# Runs every contract cell, each probe in its own fresh process, and returns
# the contract with what was observed beside each promise:
#
#   observed   "agrees" (value within tolerance of the reference), "differs",
#              "errors", or "returns" (a value, with no reference to hold it to)
#   rel_diff   the relative difference from the reference, when there is one
#   message    the error message, or the warnings, for reading
#
# Given a source directory, the estimatr processes load that tree with
# pkgload::load_all(); otherwise they attach the installed estimatr. Sourcing
# this file requires fixtures.R to have been sourced first.

contract_observe <- function(dir, source_dir = NULL) {
  contract <- read.csv(file.path(dir, "contract.csv"),
                       stringsAsFactors = FALSE, na.strings = "")
  if (!is.null(source_dir)) source_dir <- normalizePath(source_dir)
  rscript <- file.path(R.home("bin"), "Rscript")
  runner <- file.path(dir, "run_cells.R")
  # The parent's library paths, so a child under R CMD check finds the copy
  # of estimatr being checked rather than the user's. Set in this process
  # rather than through system2(env =), which Windows pastes after the
  # command as an argument
  old_libs <- Sys.getenv("R_LIBS", unset = NA)
  on.exit(if (is.na(old_libs)) Sys.unsetenv("R_LIBS") else
    Sys.setenv(R_LIBS = old_libs))
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))

  run <- function(args) {
    err <- tempfile(fileext = ".txt")
    status <- system2(rscript, c(shQuote(runner), args),
                      stdout = FALSE, stderr = err)
    if (status != 0) {
      stop("contract process failed: ", paste(args, collapse = " "), "\n",
           paste(utils::tail(readLines(err, warn = FALSE), 20), collapse = "\n"))
    }
  }

  ref_file <- tempfile(fileext = ".rds")
  run(c("reference", shQuote(ref_file)))
  refs <- readRDS(ref_file)

  est <- list()
  for (probe in unique(contract$probe)) {
    out <- tempfile(fileext = ".rds")
    run(c("estimatr", probe, shQuote(out),
          if (!is.null(source_dir)) shQuote(source_dir)))
    est[[probe]] <- readRDS(out)
  }

  contract$observed <- NA_character_
  contract$rel_diff <- NA_real_
  contract$message <- NA_character_
  for (i in seq_len(nrow(contract))) {
    e <- est[[contract$probe[i]]][[contract$fit[i]]]
    notes <- if (length(e$warnings)) paste("warning:", e$warnings, collapse = "; ")
    if (!e$ok) {
      contract$observed[i] <- "errors"
      contract$message[i] <- e$value
      next
    }
    contract$message[i] <- if (length(notes)) notes else NA_character_
    if (is.na(contract$reference[i])) {
      contract$observed[i] <- "returns"
      next
    }
    r <- refs[[paste(contract$probe[i], contract$fit[i], sep = "|")]]
    if (!r$ok) {
      contract$observed[i] <- "returns"
      contract$message[i] <- paste("reference errors:", r$value)
      next
    }
    rd <- contract_rel_diff(e$value, r$value)
    tol <- if (is.na(contract$tolerance[i])) 1e-6 else contract$tolerance[i]
    contract$rel_diff[i] <- rd
    contract$observed[i] <- if (!is.na(rd) && rd <= tol) "agrees" else "differs"
  }
  contract
}
