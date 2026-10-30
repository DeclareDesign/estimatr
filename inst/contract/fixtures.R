# The post-estimation contract's fixtures: one data set, the eleven estimatr fits
# each contract row names in its `fit` column, and the reference fits and
# variances its `reference` column is evaluated against. Sourced by
# run_cells.R, which is run once per probe in a fresh process.
#
# The reference fits come from lm(), ivreg, sandwich, and clubSandwich and are
# evaluated in a process that never loads estimatr: clubSandwich replaces
# estimatr's bread() and model.frame() methods when it is loaded second.

contract_data <- function() {
  set.seed(343)
  n <- 60
  d <- data.frame(
    x = rnorm(n), x2 = rnorm(n),
    g = rep(1:10, each = 6),
    f = factor(rep(1:4, length.out = n)),
    iv = rbinom(n, 1, 0.5)
  )
  d$z <- rbinom(n, 1, plogis(d$iv - 0.5 + 0.3 * d$x))
  d$y <- d$x + d$z + rnorm(n)
  d$y2 <- d$x2 + rnorm(n)
  d$w <- runif(n)
  d$w[c(1, 7, 20)] <- 0
  d$x_c <- d$x - mean(d$x)
  d$x2_c <- d$x2 - mean(d$x2)
  d
}

contract_fits <- function(d) list(
  lm_plain = lm_robust(y ~ x + z, d),
  lm_wt0 = lm_robust(y ~ x + z, d, weights = w),
  lm_cl = lm_robust(y ~ x + z, d, clusters = g),
  lm_fe = lm_robust(y ~ x + z, d, fixed_effects = ~f),
  lm_mv = lm_robust(cbind(y, y2) ~ x + z, d),
  lm_lin = lm_lin(y ~ z, ~x, d),
  iv_plain = iv_robust(y ~ z + x | iv + x, d),
  iv_cl = iv_robust(y ~ z + x | iv + x, d, clusters = g),
  dim = difference_in_means(y ~ z, d),
  ht = horvitz_thompson(y ~ z, d, condition_prs = c("0" = 0.5, "1" = 0.5)),
  lh = lh_robust(y ~ x + z, data = d, linear_hypothesis = "z = 0")
)

# dim's and ht's references are the regressions they equal; they are
# meaningful only where a row says so. lh's reference is the regression its
# lm_robust part fits; the hypothesis table has no lm() counterpart. lm_wt0's variances are referred to
# `ref_positive`, the same regression on the positive-weight rows, because
# sandwich counts zero-weight rows as observations.
contract_reference_fits <- function(d) list(
  lm_plain = lm(y ~ x + z, d),
  lm_wt0 = lm(y ~ x + z, d, weights = w),
  lm_cl = lm(y ~ x + z, d),
  lm_fe = lm(y ~ x + z + f, d),
  lm_mv = lm(cbind(y, y2) ~ x + z, d),
  lm_lin = lm(y ~ z + x_c + z:x_c, d),
  iv_plain = ivreg::ivreg(y ~ z + x | iv + x, data = d),
  iv_cl = ivreg::ivreg(y ~ z + x | iv + x, data = d),
  dim = lm(y ~ z, d),
  ht = lm(y ~ z, d),
  lh = lm(y ~ x + z, d)
)

contract_reference_positive <- function(d) {
  lm(y ~ x + z, d[d$w > 0, ], weights = w)
}

# The variance each estimatr fit reports by default: HC2, or CR2 when
# clustered. ht has no regression counterpart.
contract_reference_vcov <- function(ref, fit_id, d) {
  if (fit_id == "ht") return(NULL)
  if (fit_id == "lm_wt0") ref <- contract_reference_positive(d)
  if (fit_id %in% c("lm_cl", "iv_cl")) {
    return(as.matrix(clubSandwich::vcovCR(ref, cluster = d$g, type = "CR2")))
  }
  sandwich::vcovHC(ref, type = "HC2")
}

# Normal-theory t intervals on the reference variance, with the fit's
# residual degrees of freedom; CR2 fits use clubSandwich's Satterthwaite df.
contract_reference_confint <- function(ref, V, fit_id, d) {
  if (fit_id %in% c("lm_cl", "iv_cl")) {
    ci <- clubSandwich::conf_int(ref, vcov = "CR2", cluster = d$g)
    out <- as.matrix(ci[, c("CI_L", "CI_U")])
    rownames(out) <- names(coef(ref))
    return(out)
  }
  b <- contract_reference_coef(ref)
  se <- sqrt(diag(V))
  if (fit_id == "lm_wt0") ref <- contract_reference_positive(d)
  q <- qt(0.975, df.residual(ref))
  out <- cbind(b - q * se, b + q * se)
  rownames(out) <- names(b)
  out
}

# The coefficients as one vector; for several outcomes named "outcome:term",
# as vcov() names them on both lm() and estimatr fits.
contract_reference_coef <- function(ref) {
  b <- coef(ref)
  if (!is.matrix(b)) return(b)
  setNames(as.vector(b),
           paste0(rep(colnames(b), each = nrow(b)), ":", rownames(b)))
}

# The estimates and standard errors coef(summary()) reports; for several
# outcomes one table per outcome, named as summary() names them.
contract_reference_coef_table <- function(ref, V) {
  b <- coef(ref)
  if (!is.matrix(b)) return(cbind(Estimate = b, "Std. Error" = sqrt(diag(V))))
  se <- matrix(sqrt(diag(V)), nrow(b))
  tables <- lapply(seq_len(ncol(b)), function(j) {
    cbind(Estimate = b[, j], "Std. Error" = se[, j])
  })
  setNames(tables, paste("Response", colnames(b)))
}

# The robust Wald F that summary() reports for the slopes.
contract_reference_fstat <- function(ref, V, fit_id, d) {
  if (fit_id == "lm_wt0") ref <- contract_reference_positive(d)
  b <- coef(ref)[-1]
  q <- length(b)
  c(value = drop(t(b) %*% solve(V[-1, -1, drop = FALSE]) %*% b) / q,
    numdf = q, dendf = df.residual(ref))
}

# Compares an estimatr value with its reference: numeric values as matrices,
# aligned on row and column names where both have them, so that a fit
# reporting only some of lm()'s coefficients (dim, fixed effects) is compared
# on those it reports. Returns NA when the two cannot be put side by side.
contract_as_matrix <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x[vapply(x, is.numeric, TRUE)])
  # One table per outcome, stacked with the outcome in the row names, so that
  # the columns can still be matched by name
  if (is.list(x) && length(x) && all(vapply(x, is.matrix, TRUE))) {
    x <- do.call(rbind, Map(function(m, nm) {
      rownames(m) <- paste0(nm, ":", rownames(m))
      m
    }, x, names(x)))
  }
  if (is.list(x)) x <- unlist(x)
  if (!is.numeric(x) && !is.logical(x)) return(NULL)
  if (is.null(dim(x))) {
    x <- matrix(as.numeric(x), ncol = 1, dimnames = list(names(x), NULL))
  }
  storage.mode(x) <- "double"
  x
}

contract_rel_diff <- function(a, b) {
  if (is.null(a) || is.null(b)) return(if (is.null(a) && is.null(b)) 0 else NA)
  if (is.character(a) || is.character(b) || is.factor(a) || is.factor(b)) {
    return(if (identical(as.character(a), as.character(b))) 0 else NA)
  }
  A <- contract_as_matrix(a)
  B <- contract_as_matrix(b)
  if (is.null(A) || is.null(B)) return(NA)
  if (!is.null(rownames(A)) && !is.null(rownames(B))) {
    rows <- intersect(rownames(A), rownames(B))
    if (length(rows) == 0) return(NA)
    A <- A[rows, , drop = FALSE]
    B <- B[rows, , drop = FALSE]
  }
  if (!is.null(colnames(A)) && !is.null(colnames(B)) &&
      length(intersect(colnames(A), colnames(B))) > 0) {
    cols <- intersect(colnames(A), colnames(B))
    A <- A[, cols, drop = FALSE]
    B <- B[, cols, drop = FALSE]
  }
  if (!identical(dim(A), dim(B)) || length(A) == 0) return(NA)
  if (any(is.na(A) != is.na(B))) return(NA)
  fin <- is.finite(A) & is.finite(B)
  if (!any(fin)) return(0)
  max(abs(A[fin] - B[fin])) / max(abs(B[fin]), 1e-12)
}

# Evaluates one cell's expression, recording a value or an error message and
# any warnings, so that a failure in one cell never stops the others.
contract_eval <- function(expr, env) {
  warnings <- character()
  res <- withCallingHandlers(
    tryCatch(list(ok = TRUE, value = eval(parse(text = expr), env)),
             error = function(e) list(ok = FALSE, value = conditionMessage(e))),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  res$warnings <- unique(warnings)
  res
}
