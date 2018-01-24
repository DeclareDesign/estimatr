#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/


#' Extra logging on na.omit handler
#'
#' @param object a data.frame
#' @param ... unused
#'
#' @return a normal \code{omit} object, with the extra attribute \code{why_omit},
#' which contains the leftmost column containing an NA for each row that was dropped, by
#' column name, if any were dropped.
#'
#' @export
#' @seealso \code{\link{na.omit}}
na.omit_detailed.data.frame <- function(object, ...) {
  n <- length(object)
  omit <- logical(nrow(object))
  vars <- colnames(object)
  row_names <- attr(object, "row.names")
  why_omit <- list()

  for (j in vars) {
    x <- object[[j]]
    if (!is.atomic(x)) {
      next
    }
    x <- is.na(x)
    d <- dim(x)

    # special case for nested df and matrices
    if (length(d) == 2L && d[2L] > 1L) {
      for (ii in d[2L]:2L) {
        x[, ii - 1L] <- x[, ii] | x[, ii - 1L]
      }
      x <- x[, 1L]
    }

    why_omit[[j]] <- which(x & !omit)
    # TODO omit/!omit can be factored out of loop, and why_omit can be rebuilt using set logic

    omit <- omit | x
  }
  object <- object[!omit, , drop = FALSE]

  if (any(omit > 0L)) {
    temp <- setNames(which(omit), row_names[omit])
    attr(temp, "class") <- c("omit", "detailed")
    attr(temp, "why_omit") <- Filter(length, why_omit)
    attr(object, "na.action") <- temp
  }
  object
}

# NJF 10/18
# Silly microbenchmark to make sure I didn't make it slower
# df <- expand.grid(x=c(1:100, NA), y=c(1:5, NA), z=c(1:8, NA), q=c(NA,2:5))

# microbenchmark(stock=na.omit(df), hack1=na.omit_detailed.data.frame(df))
## Unit: milliseconds
## expr      min       lq     mean   median       uq        max neval
## stock 6.114132 6.184318 7.744881 6.232744 6.961491 101.823530   100
## hack1 5.360638 5.480531 6.525075 5.694078 7.752104   9.323943   100
