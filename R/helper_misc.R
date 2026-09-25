add_cis_pvals <- function(return_frame, alpha, ci, ttest = TRUE) {
  if (ci) {
    if (alpha <= 0 || alpha >= 1) {
      stop("`alpha` must be numeric between 0 and 1")
    }

    return_frame$statistic <- with(return_frame, coefficients / std.error)

    if (ttest) {
      if (any(return_frame$df <= 0, na.rm = TRUE)) {
        warning(
          "Some degrees of freedom have been estimated as negative or zero\n",
          "p-values and confidence intervals may not be calculated"
        )

        return_frame$df <- ifelse(return_frame$df <= 0, NA, return_frame$df)
      }

      return_frame$p.value <- with(
        return_frame,
        2 * pt(abs(statistic), df = df, lower.tail = FALSE)
      )

      crit_se <- with(return_frame, qt(1 - alpha / 2, df = df) * std.error)
    } else {
      return_frame$p.value <- with(
        return_frame,
        2 * pnorm(abs(statistic), lower.tail = FALSE)
      )

      crit_se <- with(return_frame, qnorm(1 - alpha / 2) * std.error)

      return_frame$df <- NA
    }

    return_frame$conf.low <- with(return_frame, coefficients - crit_se)
    return_frame$conf.high <- with(return_frame, coefficients + crit_se)

    return(as.list(return_frame))
  } else {
    return_frame$p.value <- NA
    return_frame$statistic <- NA
    return_frame$conf.low <- NA
    return_frame$conf.high <- NA

    return(as.list(return_frame))
  }
}

lm_return <- function(return_list, model_data, formula,
                      lin_interactions = NULL, lin_treatment = NULL) {

  # A collinear column is dropped and comes back as an NA coefficient. Saying
  # so is the difference between a user reading the NA correctly and reading
  # it as a bug, since the remaining coefficients are then conditional on a
  # different set of regressors than they asked for (estimatr #411). It is a
  # message and not a warning because stats::lm() signals neither, and because
  # a warning raised inside a grouped dplyr verb nested in another one crashes
  # dplyr 1.2.1 while it builds the warning's group label.
  coefs <- return_list[["coefficients"]]
  na_coefs <- is.na(coefs)
  if (any(na_coefs)) {
    dropped <- if (is.matrix(coefs))
      rownames(coefs)[apply(na_coefs, 1, any)]
      else names(coefs)[na_coefs]
    msg <- paste0(
      "Some coefficients are collinear with other regressors and were ",
      "dropped, and are returned as NA: ",
      paste(dropped, collapse = ", "), "."
    )

    # A dropped treatment interaction is not the same event as a dropped
    # covariate, and naming the column alone does not say so. lm_lin()'s
    # treatment coefficient is the effect at the covariate means only while
    # every covariate is interacted; once an interaction goes, that covariate's
    # slope is constrained equal across arms, and where the dependency runs
    # through the treatment column itself the reported value is the one the
    # drop rule happened to leave behind. On one published replication the
    # intercept, the treatment indicator, a covariate and its interaction are
    # an exact four-column dependency in every subgroup cell where a category
    # is empty in one arm; the fitted values and r.squared agree to twelve
    # digits under either drop while the treatment coefficient moves from
    # -1.56 to -2.00, so no goodness-of-fit channel can see it. `lin_interactions` names each
    # interaction column, the treatment column it was built from, and the
    # covariate it centres, and is NULL everywhere but lm_lin().
    if (!is.null(lin_interactions)) {
      # Only where the covariate's own main effect SURVIVES. When a covariate
      # is a duplicate or is constant, its main effect and its interaction drop
      # together, the remaining fit is lm_lin() on the covariates that are
      # left, and the treatment coefficient is still the effect at their means.
      lost <- lin_interactions$interaction %in% dropped &
        !(lin_interactions$covariate %in% dropped)
      affected <- setdiff(unique(lin_interactions$treatment[lost]), dropped)
      if (length(affected) > 0) {
        lost_covs <- unique(lin_interactions$covariate[lost])
        # Named while the list is short enough to read. A saturated factor can
        # lose twenty interactions at once, and the first sentence has already
        # printed every one of them.
        covs_txt <- if (length(lost_covs) > 3) {
          paste(length(lost_covs), "covariates")
        } else {
          paste(lost_covs, collapse = ", ")
        }
        msg <- paste0(
          msg,
          " A dropped treatment interaction is not a dropped covariate: with ",
          covs_txt,
          " no longer interacted with the treatment, ",
          paste(affected, collapse = ", "),
          if (length(affected) > 1) {
            " are no longer effects at the covariate means, since "
          } else {
            " is no longer the effect at the covariate means, since "
          },
          if (length(lost_covs) > 1) "their slopes are " else "its slope is ",
          "now constrained equal across treatment arms. Where the dependency ",
          "runs through the treatment column itself, the reported value also ",
          "depends on which column was dropped, and the effect at the ",
          "covariate means is not identified."
        )
      }
    }

    message(msg)

    # A dropped treatment indicator is the one case here that is worth the
    # louder channel. The other drops leave a fit whose treatment coefficient
    # means something narrower than it did; this one leaves no treatment
    # coefficient at all, under a name that is still printed with an NA beside
    # it, and there is nothing a caller can read off the fit instead.
    dropped_treat <- intersect(dropped, lin_treatment)
    if (length(dropped_treat) > 0) {
      warning(
        "lm_lin() dropped the treatment indicator ",
        paste(dropped_treat, collapse = ", "),
        ", so the fit carries no treatment effect for ",
        if (length(dropped_treat) > 1) "those arms" else "that arm",
        ". A treatment indicator is collinear with the rest of the design ",
        "only where its arm is empty or where the covariates reproduce it ",
        "exactly, and neither is a fit to read an effect from.",
        call. = FALSE
      )
    }
  }

  if (!is.null(model_data)) {
    return_list[["contrasts"]] <- attr(model_data$design_matrix, "contrasts")
    return_list[["terms"]] <- model_data$terms
    return_list[["xlevels"]] <- model_data$xlevels
    return_list[["felevels"]] <- model_data$felevels
    return_list[["weights"]] <- model_data$weights
    if (is.matrix(model_data$outcome) &&
        is.character(colnames(model_data$outcome))) {
      return_list[["outcome"]] <- colnames(model_data$outcome)
    } else {
      return_list[["outcome"]] <- deparse(formula[[2]], nlines = 5)
    }
  }

  if (is.matrix(return_list[["std.error"]]) &&
      ncol(return_list[["std.error"]]) > 1) {
    dimnames(return_list[["std.error"]]) <- dimnames(return_list[["coefficients"]])
  } else {
    return_list[["coefficients"]] <- drop(return_list[["coefficients"]])
    nms <- c("std.error", "statistic", "p.value", "df", "conf.low", "conf.high")
    for (nm in nms) {
      if (length(return_list[[nm]]) > 1 || !is.na(return_list[[nm]])) {
        return_list[[nm]] <- setNames(
          drop(return_list[[nm]]),
          names(return_list[["coefficients"]])
        )
      }
    }
  }
  # The observation names were taken off the design matrix in
  # clean_model_data() so that nothing in the fit had to carry them, and this
  # is where they go back on. `fitted.values` is the only field that keeps
  # them, as in 1.0.6; `weights` is named to match, as it was before.
  obs_names <- if (!is.null(model_data)) model_data[["obs_names"]] else NULL
  if (return_list[["weighted"]] && !is.null(return_list[["weights"]])) {
    names(return_list[["weights"]]) <- obs_names
  }
  fitted <- drop(return_list[["fitted.values"]])
  if (!is.null(obs_names)) {
    if (is.matrix(fitted)) rownames(fitted) <- obs_names else names(fitted) <- obs_names
  }
  return_list[["fitted.values"]] <- fitted
  return_list[["ei.iv"]] <- drop(return_list[["ei.iv"]])
  # Residuals come back unnamed, as they did in 1.0.6. They no longer acquire
  # names to lose: the design matrix arrives here without row names at all.
  return_list[["residuals"]] <- drop(return_list[["residuals"]])
  return(return_list)
}

dim_like_return <- function(return_list, alpha, formula, conditions) {
  return_list[["alpha"]] <- alpha

  treat_condition <- conditions[[2]]

  add_label <- !(conditions[[2]] == 1 && conditions[[1]] == 0)

  fterms <- terms(formula)
  coef_name <- labels(fterms)

  if (add_label) {
    return_list[["term"]] <- paste0(
      coef_name,
      conditions[[2]]
    )
  } else {
    return_list[["term"]] <- coef_name
  }

  return_list[["outcome"]] <- deparse(formula[[2]], nlines = 5)

  names(return_list[["coefficients"]]) <-
    names(return_list[["std.error"]]) <-
    names(return_list[["p.value"]]) <-
    names(return_list[["df"]]) <- return_list[["term"]]

  return_list[["condition2"]] <- conditions[[2]]
  return_list[["condition1"]] <- conditions[[1]]

  return_list[["vcov"]] <- matrix(
    data = return_list[["std.error"]] ^ 2,
    dimnames = list(return_list[["term"]], return_list[["term"]])
  )

  return(return_list)
}

parse_conditions <- function(treatment, condition1, condition2, estimator) {
  if (is.factor(treatment)) {
    condition_names <- levels(droplevels(treatment))
  } else {
    condition_names <- sort(unique(treatment))
  }

  if (any(!(c(condition1, condition2) %in% condition_names))) {
    stop("`condition1` and `condition2` must be values found in the treatment")
  }

  n_conditions <- length(condition_names)

  conditions <- list(NULL, NULL)

  if (n_conditions > 2) {
    if (is.null(condition1) || is.null(condition2)) {
      stop(
        "Treatment has > 2 values; must specify both `condition1` and ",
        "`condition2` or use a treatment with only 2 values"
      )
    } else {
      conditions[1:2] <- c(condition1, condition2)
    }
  } else if (n_conditions == 2) {
    if (is.null(condition1) && is.null(condition2)) {
      conditions[1:2] <- condition_names
    } else if (!is.null(condition2) && is.null(condition1)) {
      conditions[1:2] <- c(setdiff(condition_names, condition2), condition2)
    } else if (!is.null(condition1) && is.null(condition2)) {
      conditions[1:2] <- c(condition1, setdiff(condition_names, condition1))
    } else {
      conditions[1:2] <- c(condition1, condition2)
    }
  } else if (n_conditions == 1) {
    stop(
      "Must have more than one value in treatment"
    )
  }

  return(conditions)
}

check_clusters_blocks <- function(data) {
  if (!is.null(data$cluster)) {
    one_block_per_clust <-
      tapply(data$block, data$cluster, function(x) all(x == x[1]))

    if (any(!one_block_per_clust)) {
      stop("All `clusters` must be contained within `blocks`")
    }

    clust_per_block <- tapply(
      data$cluster,
      data$block,
      function(x) length(unique(x))
    )
  } else {
    clust_per_block <- tabulate(as.factor(data$block))
  }

  return(clust_per_block)
}

#' The `_c` names lm_lin gives its centred covariates
#'
#' Shared by [lm_lin()], which builds the design, and `predict.lm_robust()`,
#' which rebuilds it for `newdata`. The two must agree exactly, since predict
#' lines the design up against the coefficients by name; keeping one copy of
#' the rule is what stops them drifting apart.
#'
#' A name containing an interaction or a function call is parenthesised first,
#' so `poly(x, 2)1` becomes `(poly(x, 2)1)_c` rather than `poly(x, 2)1_c`.
#'
#' @param x Character vector of covariate column names.
#' @return Character vector of the same length.
#' @keywords internal
#' @noRd
lin_covar_names <- function(x) {
  paste0(
    ifelse(grepl("\\:|(^.+\\()", x), paste0("(", x, ")"), x),
    "_c"
  )
}

#' Recover the absorbed group effects of a one-way fixed-effects fit
#'
#' `predict()` needs them to put the group contribution back for `newdata`, and
#' they are otherwise unrecoverable once the design matrix has been demeaned.
#' `fitted - X'b` is constant within group by construction, so one value per
#' group is the whole answer.
#'
#' Only defined for one-way FE: with several factors the sum is identified but
#' the parts are not, so there is nothing to map a `newdata` row onto. Returns
#' `NULL` in that case, and for multivariate outcomes.
#'
#' @param fitted Full-model fitted values.
#' @param coefficients The fitted coefficient vector.
#' @param model_data The object returned by `demean_fes()`, which holds
#'   `Xoriginal` and the undemeaned `fixed_effects` matrix.
#' @return A named numeric vector, one entry per group, or `NULL`.
#' @keywords internal
#' @noRd
absorbed_group_effects <- function(fitted, coefficients, model_data) {
  fe <- model_data[["fixed_effects"]]
  if (ncol(fe) != 1L || NCOL(fitted) != 1L) return(NULL)
  # The value is constant within group, so the group effect is the value at
  # the group's first row. tapply() would build the full split to read one
  # element out of each piece; match() finds those rows in a single pass, and
  # matching integer codes rather than the level strings is three times faster
  # again at n = 100,000.
  g <- fe[, 1L]
  lv <- model_data[["fe_level_names"]][[1L]]
  # 1.0.6 built these from a character matrix, so it returned them in string
  # order: a factor with levels 1..30 came back 1, 10, 11, ... The codes carry
  # the data's own order instead, so sort here rather than silently renaming
  # and reordering a returned field. Sorting `length(lv)` strings is nothing;
  # the string matching over all n observations, which this replaces, was not.
  ord <- order(as.character(lv))
  # Evaluate at the representative rows only. The group effect is constant
  # within group, so forming `fitted - X b` for all n observations to read one
  # value per group out of it does n-by-k work, and allocates two n-length
  # vectors, to produce `length(lv)` numbers.
  # `match(ord, g)` hashes all n codes to find one row per group. The codes
  # are contiguous 1..L, so a scatter assignment gets a representative row in
  # one pass and no hash table. It lands on the group's LAST row rather than
  # its first, which is the same number: the group effect is constant within
  # group, which is the whole reason this is evaluated at one row at all.
  rep_row <- integer(length(lv))
  rep_row[g] <- seq_along(g)
  # A level with no rows keeps its 0 from the initialisation, and 0 is a
  # DROP index in R, not a missing one: it would silently shorten the result
  # rather than leave a hole. `match()` returned NA here, and the effects come
  # back NA for that level with the vector still one entry per level, which is
  # what `dim<-` below and every consumer expect. Reachable whenever a factor
  # carries a declared level the data never uses.
  rep_row[rep_row == 0L] <- NA_integer_
  idx <- rep_row[ord]
  # `coefficients` carries an NA wherever a column was dropped for collinearity,
  # and NA times anything is NA, so multiplying the full `Xoriginal` by the full
  # vector made every group effect NA on any rank-deficient fit: the same defect
  # b96ab8b fixed in the F statistic, in a different field. A dropped column
  # contributes nothing to the fitted values, so both sides drop it.
  keep <- !is.na(coefficients)
  effects <- fitted[idx] -
    drop(model_data[["Xoriginal"]][idx, keep, drop = FALSE] %*% coefficients[keep])
  lv <- lv[ord]
  # tapply() returned a one-dimensional array, and that shape is part of the
  # return surface: `predict()` and the 1.0.6 comparison both see it.
  dim(effects) <- length(lv)
  dimnames(effects) <- list(paste0(colnames(fe), lv))
  effects
}

# Full-model R^2 for an absorbed fixed-effects fit, on the original outcome.
#
# The FE branches used `mean(yoriginal)` and raw residuals, so a weighted fit
# reported an unweighted R^2 (0.6564 against 0.6432 from the same model with
# dummies and from `lm()`), and a multivariate outcome was pooled into one
# number where the dummy fit gives one per column. `w` is the RAW weight
# vector, as `demean_fes()` leaves it, or NULL.
fe_r2 <- function(yoriginal, residuals_proj, w) {
  y <- as.matrix(yoriginal)
  e <- as.matrix(residuals_proj)
  n <- nrow(y)
  if (is.null(w)) {
    ymw <- y - rep(colMeans(y), each = n)
    list(tss = colSums(ymw * ymw), rss = colSums(e * e))
  } else {
    ymw <- y - rep(drop(crossprod(w, y) / sum(w)), each = n)
    list(tss = colSums(ymw * ymw * w), rss = colSums(e * e * w))
  }
}
