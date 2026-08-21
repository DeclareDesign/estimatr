#' @importFrom rlang f_rhs %||% quo_is_missing eval_tidy sym as_label
#' @importFrom stats terms model.matrix model.response model.extract
clean_model_data <- function(data, datargs, estimator = "") {

  data <- if (quo_is_missing(data)) NULL else eval_tidy(data)

  mfargs <- Filter(Negate(quo_is_missing), datargs)

  m_formula <- eval_tidy(mfargs[["formula"]])
  m_formula_env <- environment(m_formula)

  stopifnot("`formula` argument must be a formula" = inherits(m_formula, "formula"))

  mfargs <- as.list(mfargs)

  args_ignored <- c("se_type", "fixed_effects")
  to_process <- setdiff(
    names(mfargs),
    c(
      setdiff(names(formals(stats::model.frame.default)), "subset"),
      args_ignored
    )
  )

  for (da in to_process) {
    name <- sprintf(".__%s%%%d__", da, sample.int(.Machine$integer.max, 1))
    m_formula_env[[name]] <- eval_tidy(mfargs[[da]], data = data)
    mfargs[[da]] <- sym(name)
  }

  # fixed_effects: evaluate and store as factor matrix in formula env so
  # model.frame can attach it as an ancillary variable
  if ("fixed_effects" %in% names(mfargs)) {
    name <- sprintf(".__fixed_effects%%%d__", sample.int(.Machine$integer.max, 1))
    fe_formula <- eval_tidy(mfargs[["fixed_effects"]], data = data)
    # 1.x also accepted a bare column name or an already-evaluated grouping
    # vector. estimatr issue #304 asked for the formula to be enforced; a
    # warning that names the argument and the expected form does that without
    # breaking working code, so the vector is still accepted. `RCT` on CRAN
    # passes one. Anything that is neither a formula nor a grouping vector is
    # still an error.
    if (!inherits(fe_formula, "formula")) {
      if (!(is.atomic(fe_formula) || is.factor(fe_formula)) ||
          is.matrix(fe_formula) || is.null(fe_formula)) {
        stop(
          "`fixed_effects` must be a one-sided formula, such as `~ blockID` or ",
          "`~ block + year`. You passed an object of class ",
          paste(class(fe_formula), collapse = "/"), "."
        )
      }
      # The term name has to come from the unevaluated expression, since the
      # value is a bare vector with no name of its own. 1.0.6 named it by
      # deparsing the same expression, so `felevels` and the absorbed
      # `fixed_effects` keep the names they had.
      fe_name <- as_label(mfargs[["fixed_effects"]])
      warning(
        "Passing `fixed_effects` a bare grouping vector is deprecated and will ",
        "be removed in a future version. Use a one-sided formula instead:\n",
        "  fixed_effects = ~ ", fe_name,
        call. = FALSE
      )
      fe_mf <- stats::setNames(
        data.frame(as.factor(fe_formula), stringsAsFactors = TRUE),
        fe_name
      )
    } else {
      fe_mf <- stats::model.frame.default(fe_formula, data = data, na.action = NULL)
    }
    # The fixed effects travel as an INTEGER CODE matrix, with the level names
    # carried alongside. An earlier version put a character matrix here, which
    # every consumer then coerced back to factors: a round trip that cost
    # around 4 ms per fit at n = 100,000, doubled the memory the model frame
    # had to carry through NA handling, and re-derived the levels by string
    # sort, so a factor with levels 1..30 came back ordered 1, 10, 11. Codes
    # keep the order the data had, so nothing has to be corrected afterwards.
    fe_fac <- lapply(fe_mf, function(x) if (is.factor(x)) x else as.factor(x))
    fe_level_names <- lapply(fe_fac, levels)
    m_formula_env[[name]] <- sapply(fe_fac, as.integer)
    mfargs[["fixed_effects"]] <- sym(name)
  } else {
    fe_level_names <- NULL
  }

  # IV needs Formula (for the | separator); everything else uses plain formula
  # and avoids the ~80µs overhead of Formula::as.Formula + model.frame.Formula
  if (estimator == "iv") {
    mfargs[["formula"]] <- Formula::as.Formula(m_formula)
  } else {
    mfargs[["formula"]] <- m_formula
  }

  mf <- eval_tidy(rlang::quo((stats::model.frame)(
    !!!mfargs,
    data = data,
    na.action = na_omit_detailed,
    drop.unused.levels = TRUE
  )))

  local({
    na.action <- attr(mf, "na.action")
    why_omit <- attr(na.action, "why_omit")

    missing_warning <- c(
      "Some observations have missingness in the %s variable(s) but not in ",
      "the outcome or covariates. These observations have been dropped."
    )

    to_check_if_missing <- c("cluster", "block", "weights", "fixed_effects")

    for (x in to_check_if_missing) {
      if (!is.null(why_omit[[sprintf("(%s)", x)]])) {
        warning(sprintf(missing_warning, x))
      }
    }
  })

  if (!is.null(attr(terms(mf), "Formula_without_dot"))) {
    formula <- attr(terms(mf), "Formula_without_dot")
  } else {
    formula <- eval_tidy(mfargs[["formula"]])
  }

  # For IV, select only the regressor RHS (rhs=1) from the Formula object.
  # For non-IV, use terms(mf) which already has '.' expanded against the data.
  design_terms <- if (estimator == "iv") terms(formula, rhs = 1) else terms(mf)
  ret <- list(
    outcome = stats::model.response(mf, type = "numeric"),
    design_matrix = stats::model.matrix(design_terms, data = mf),
    formula = formula
  )

  if (estimator == "iv") {
    if (length(formula)[2] != 2) {
      stop(
        "Must specify a `formula` with both regressors and instruments. For ",
        "example, `formula = y ~ x1 + x2 | x1 + z2`.\n\nSee ?iv_robust."
      )
    }
    ret[["instrument_matrix"]] <- stats::model.matrix(terms(formula, rhs = 2), data = mf)
    ret[["terms_regressors"]] <- terms(formula, rhs = 1)
  } else if (estimator == "dim") {
    ret[["original_treatment"]] <- mf[, colnames(mf) == all.vars(terms(mf)[[3]])[1]]
  }

  ret[["weights"]] <- stats::model.extract(mf, "weights")
  if (any(ret[["weights"]] < 0)) {
    stop("`weights` must not be negative")
  }

  ret[["cluster"]] <- stats::model.extract(mf, "cluster")
  if (!is.null(ret[["cluster"]]) &&
      !is.factor(ret[["cluster"]]) &&
      !is.integer(ret[["cluster"]])) {
    ret[["cluster"]] <- as.factor(ret[["cluster"]])
  }

  ret[["block"]] <- stats::model.extract(mf, "block")

  ret[["fixed_effects"]] <- stats::model.extract(mf, "fixed_effects")
  if (!is.null(ret[["fixed_effects"]])) {
    ret[["fixed_effects"]] <- as.matrix(ret[["fixed_effects"]])
    # A level can disappear between capture and here, either because every one
    # of its rows was dropped for missingness or because the factor arrived
    # carrying levels it never used. Downstream code sizes a rank correction
    # off the number of levels, so an absent level must not be counted, and
    # fe_leverage() indexes by code and needs them contiguous. Renumber only
    # the columns that actually have a gap.
    for (k in seq_len(ncol(ret[["fixed_effects"]]))) {
      col <- ret[["fixed_effects"]][, k]
      u <- sort(unique(col))
      if (length(u) != u[length(u)]) {
        ret[["fixed_effects"]][, k] <- match(col, u)
        fe_level_names[[k]] <- fe_level_names[[k]][u]
      }
    }
  }
  ret[["fe_level_names"]] <- fe_level_names

  ret[["terms"]] <- attr(mf, "terms")
  dcs <- attr(ret[["terms"]], "dataClasses")
  drop_vars <- c("(block)", "(cluster)", "(fixed_effects)")
  attr(ret[["terms"]], "dataClasses") <- dcs[setdiff(names(dcs), drop_vars)]
  ret[["xlevels"]] <- .getXlevels(ret[["terms"]], mf)

  return(ret)
}

# The fixed effects as a data frame of factors, rebuilt from the integer codes
# and the level names. Only `model.matrix()` needs factors, so this is on the
# CR2 path and the instrument-demeaning path and nowhere else; everything on
# the ordinary path works from the codes directly.
#' @noRd
fe_factors <- function(model_data) {
  codes <- model_data[["fixed_effects"]]
  lv <- model_data[["fe_level_names"]]
  out <- lapply(seq_len(ncol(codes)), function(k) {
    factor(codes[, k], levels = seq_along(lv[[k]]), labels = lv[[k]])
  })
  stats::setNames(as.data.frame(out, stringsAsFactors = FALSE), colnames(codes))
}

# Demean outcome and design matrix by fixed effects (Frisch-Waugh-Lovell).
# FE demeaning via alternating projections in C++ (demean_cpp in lm_robust_helper.cpp).
# One-way FE converges in exactly 1 iteration; multi-way iterates to eps = 1e-8.
# Returns model_data with demeaned outcome and design matrix (intercept absorbed).
#' @importFrom stats ave
#' @noRd
demean_fes <- function(model_data) {
  codes     <- model_data[["fixed_effects"]]
  fe_codes  <- lapply(seq_len(ncol(codes)), function(k) codes[, k])
  fe_levels <- vapply(model_data[["fe_level_names"]], length, 0L)
  # At this point model_data[["weights"]] is the RAW weight vector (not yet
  # sqrt-transformed — that happens in lm_robust_fit / prep_lm_data).
  # WLS group means use raw weights as the metric.
  w       <- model_data[["weights"]] %||% numeric(0L)
  has_int <- attr(model_data$terms, "intercept")

  y_dm <- demean_cpp(as.matrix(model_data[["outcome"]]), fe_codes, w)
  colnames(y_dm) <- colnames(model_data[["outcome"]])
  model_data[["outcome"]] <- y_dm

  X         <- model_data[["design_matrix"]]
  keep_cols <- if (has_int) colnames(X) != "(Intercept)" else rep(TRUE, ncol(X))
  X_sub     <- X[, keep_cols, drop = FALSE]
  # Kept so lm_robust can recover the absorbed group effects, which predict()
  # needs and which are otherwise unrecoverable once X is demeaned.
  model_data[["Xoriginal"]] <- X_sub
  model_data[["design_matrix"]] <- if (ncol(X_sub) > 0L) {
    X_dm <- demean_cpp(X_sub, fe_codes, w)
    colnames(X_dm) <- colnames(X_sub)
    X_dm
  } else {
    X_sub
  }

  model_data[["fe_levels"]] <- fe_levels
  # `fe_level_names` needs no repair here any more. It is captured in the order
  # the data had and renumbered alongside the codes in clean_model_data(), so
  # by this point its order and its membership are both already right.
  # `fe_codes` and the raw weights are what fe_leverage() needs, and it is
  # called only when the requested se_type wants it.
  model_data[["fe_codes"]] <- fe_codes
  model_data
}

# Per-observation leverage contributed by the absorbed fixed effects, exactly:
# the diagonal of P_D, the projection onto the span of the FE dummies.
#
# The hat values of the full design decompose as
#
#   P_[X | D] = P_D + P_{M_D X}
#
# so h_ii = (this vector) + the hat value of the DEMEANED X, which the fitter
# already computes. That is what makes HC2 and HC3 available under
# `fixed_effects` without ever materialising the n-by-g dummy matrix.
#
# For ONE factor P_D is diagonal and this is just w_i / sum(w in group i).
#
# For TWO OR MORE it is not, but it is still cheap. Write D with the factor
# that has the MOST levels in full dummies and the rest contrast-coded; then
# A = D'WD has a diagonal leading block, and block inversion needs only the
# Schur complement S, which is `sum_{k>1}(g_k - 1)` square -- as small as the
# design allows. `col(D)`, and so P_D, does not depend on that choice.
#
# S is inverted through its eigendecomposition rather than a solve, so a
# DISCONNECTED design (where D is rank deficient) gives the pseudo-inverse and
# the right projection anyway.
#
# An earlier version of this file recorded that no such decomposition existed
# beyond one factor. What had been tested was the sum of each factor's own
# within-group share, which ignores the cross terms and is wrong in the third
# decimal place. The identity above is exact; see test_fe_leverage.R.
# Returns a list with `rank`, the exact rank of the FE design, and `leverage`,
# `diag(P_D)` -- or NULL for the leverage when `leverage = FALSE`, since only
# HC2 and HC3 need the vector while every se_type needs the rank, and the
# vector's per-observation gathers are the expensive part.
#' @noRd
fe_leverage <- function(fe_codes, w = NULL, leverage = TRUE) {
  n <- length(fe_codes[[1L]])
  if (is.null(w) || length(w) == 0L) w <- rep(1, n)
  g <- vapply(fe_codes, max, 0L)

  ord <- order(g, decreasing = TRUE)
  fe_codes <- fe_codes[ord]
  g <- g[ord]

  a <- fe_codes[[1L]]
  g1 <- unname(g[1L])
  A11 <- as.vector(rowsum(w, a, reorder = TRUE))

  K <- length(g)
  offs <- if (K > 1L) cumsum(c(0L, g[-1L] - 1L)) else 0L
  p <- offs[K]

  # `p` is the number of contrast columns the factors after the first
  # contribute. A factor with a single level contributes none: it is an
  # intercept, and the first factor's full dummies already span it. When every
  # later factor is like that -- and when there is only one factor at all --
  # the design reduces to the one-way case, where P_D is diagonal.
  if (p == 0L) {
    return(list(rank = g1,
                leverage = if (leverage) w / A11[a] else NULL))
  }
  ck <- lapply(2:K, function(k) fe_codes[[k]])
  keep <- lapply(ck, function(z) z > 1L)
  colof <- lapply(seq_len(K - 1L), function(k) offs[k] + ck[[k]] - 1L)

  # Weighted cross-tabulation of two integer codes, as (index, value) pairs.
  xtab <- function(i1, i2, n1, wts) {
    s <- rowsum(wts, i1 + (i2 - 1L) * n1, reorder = TRUE)
    list(idx = as.integer(rownames(s)), val = as.vector(s))
  }

  A12 <- matrix(0, g1, p)
  for (k in seq_len(K - 1L)) {
    kp <- keep[[k]]
    s <- xtab(a[kp], ck[[k]][kp] - 1L, g1, w[kp])
    A12[cbind((s$idx - 1L) %% g1 + 1L, offs[k] + (s$idx - 1L) %/% g1 + 1L)] <- s$val
  }

  A22 <- matrix(0, p, p)
  for (k in seq_len(K - 1L)) {
    for (l in k:(K - 1L)) {
      kp <- keep[[k]] & keep[[l]]
      if (!any(kp)) next
      nk <- g[k + 1L] - 1L
      s <- xtab(ck[[k]][kp] - 1L, ck[[l]][kp] - 1L, nk, w[kp])
      rr <- offs[k] + (s$idx - 1L) %% nk + 1L
      cc <- offs[l] + (s$idx - 1L) %/% nk + 1L
      A22[cbind(rr, cc)] <- A22[cbind(rr, cc)] + s$val
      if (l != k) A22[cbind(cc, rr)] <- A22[cbind(cc, rr)] + s$val
    }
  }

  C <- A12 / A11
  es <- eigen(A22 - crossprod(A12, C), symmetric = TRUE)
  tol <- max(es$values, 0) * 1e-12
  # rank(D) = g1 + rank(S): the leading block of A is diagonal and positive, so
  # it contributes g1, and the Schur complement contributes the rest. The
  # eigenvalues are already here, so the exact rank of a nested or disconnected
  # FE design costs nothing extra. The nominal `sum(levels) - K + 1` overstates
  # it whenever one factor is partly spanned by the others.
  rank_S <- sum(es$values > tol)
  if (!leverage) {
    return(list(rank = g1 + rank_S, leverage = NULL))
  }
  S_inv <- es$vectors %*% (ifelse(es$values > tol, 1 / es$values, 0) * t(es$vectors))
  Tm <- C %*% S_inv
  q <- rowSums(Tm * C)

  cross <- numeric(n)
  quad <- numeric(n)
  for (k in seq_len(K - 1L)) {
    kp <- keep[[k]]
    cross[kp] <- cross[kp] + Tm[cbind(a[kp], colof[[k]][kp])]
    for (l in seq_len(K - 1L)) {
      kl <- kp & keep[[l]]
      if (any(kl)) {
        quad[kl] <- quad[kl] + S_inv[cbind(colof[[k]][kl], colof[[l]][kl])]
      }
    }
  }

  list(rank = g1 + rank_S,
       leverage = w * (1 / A11[a] + q[a] - 2 * cross + quad))
}

# HC2, HC3 and CR2 are built from the hat values of the full [X | FE dummies]
# design. Under a single factor those decompose into `fe_leverage` (see
# demean_fes) and no dummy matrix is needed, but under two or more factors, and
# for CR2 under any number, there is no such shortcut: the dummies have to be
# materialised, exactly as estimatr 1.0.6 did. That costs roughly O(g^3) in the
# number of levels, so it is built only when the requested `se_type` actually
# needs it and never to satisfy a default.
#' @importFrom stats model.matrix
#' @noRd
fe_dummy_matrix <- function(model_data) {
  fe_df     <- fe_factors(model_data)
  fe_levels <- vapply(fe_df, nlevels, 0L)

  # A factor with one level contributes an intercept, not a set of contrasts,
  # and model.matrix() would return a zero-column matrix for it. Alongside
  # other factors its span is already inside theirs, so it is dropped; on its
  # own the intercept column is the whole design. fe_leverage() takes the same
  # view, so every se_type answers the same question.
  if (any(fe_levels == 1L)) {
    if (all(fe_levels == 1L)) {
      return(matrix(
        1,
        nrow  = nrow(fe_df),
        ncol  = 1L,
        dimnames = list(rownames(fe_df),
                        paste0(names(fe_df)[1L], levels(fe_df[[1L]])))
      ))
    }
    fe_df <- fe_df[, fe_levels > 1L, drop = FALSE]
  }

  model.matrix(~ 0 + ., data = fe_df)
}

# Whether a fit needs the dummy matrix above. `se_type` here is the value the
# user asked for, before check_se_type() fills in a default: a default must
# never trigger the expansion, which is what keeps absorption fast.
#' @noRd
needs_fe_dummies <- function(se_type, model_data) {
  # HC2 and HC3 come from fe_leverage() for any number of factors. CR2 needs
  # the whole per-cluster block of the hat matrix, not just its diagonal, and
  # has no such decomposition, so it is the only case left that expands.
  isTRUE(se_type == "CR2")
}

# Whether the fit wants the leverage vector. A NULL `se_type` is the default,
# which under `fixed_effects` is HC2, so it does.
#' @noRd
needs_fe_leverage <- function(se_type) {
  is.null(se_type) || se_type %in% c("HC2", "HC3")
}

# Demean an arbitrary matrix by the same FE structure as model_data.
# Used in iv_robust to demean the instrument matrix.
demean_matrix_by_fes <- function(mat, model_data) {
  codes    <- model_data[["fixed_effects"]]
  fe_codes <- lapply(seq_len(ncol(codes)), function(k) codes[, k])
  w        <- model_data[["weights"]] %||% numeric(0L)
  has_int  <- !is.null(colnames(mat)) && "(Intercept)" %in% colnames(mat)

  orig_colnames <- colnames(mat)
  mat <- demean_cpp(mat, fe_codes, w)
  colnames(mat) <- orig_colnames
  if (has_int) mat <- mat[, colnames(mat) != "(Intercept)", drop = FALSE]
  mat
}
