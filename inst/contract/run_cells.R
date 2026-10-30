# Evaluates contract cells in the process it is run in. Usage:
#
#   Rscript run_cells.R estimatr <probe> <out.rds> [<source dir>]
#   Rscript run_cells.R reference <out.rds>
#
# The estimatr mode evaluates one probe's `call` on each of the eleven fits, with
# nothing loaded but estimatr and whatever namespace that call reaches. Given a
# source directory it loads estimatr with pkgload::load_all(), so that a
# development session tests its own tree rather than the installed copy.
# The reference mode never loads estimatr.

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value = TRUE))))
source(file.path(here, "fixtures.R"))
contract <- read.csv(file.path(here, "contract.csv"), stringsAsFactors = FALSE,
                     na.strings = "")
d <- contract_data()

if (args[1] == "estimatr") {
  if (length(args) >= 4) {
    suppressMessages(pkgload::load_all(args[4], quiet = TRUE,
                                       export_all = FALSE))
  } else {
    suppressMessages(library(estimatr))
  }
  fits <- contract_fits(d)
  rows <- contract[contract$probe == args[2], ]
  out <- list()
  for (i in seq_len(nrow(rows))) {
    env <- new.env(parent = globalenv())
    env$fit <- fits[[rows$fit[i]]]
    out[[rows$fit[i]]] <- contract_eval(rows$call[i], env)
  }
  saveRDS(out, args[3])
} else {
  refs <- contract_reference_fits(d)
  V <- lapply(setNames(names(refs), names(refs)), function(id) {
    tryCatch(contract_reference_vcov(refs[[id]], id, d), error = function(e) NULL)
  })
  rows <- contract[!is.na(contract$reference), ]
  out <- list()
  for (i in seq_len(nrow(rows))) {
    id <- rows$fit[i]
    env <- new.env(parent = globalenv())
    env$fit <- env$ref <- refs[[id]]
    env$ref_positive <- contract_reference_positive(d)
    env$V <- V[[id]]
    env$fit_id <- id
    expr <- if (rows$reference[i] == "same") rows$call[i] else rows$reference[i]
    out[[paste(rows$probe[i], id, sep = "|")]] <- contract_eval(expr, env)
  }
  saveRDS(out, args[2])
}
