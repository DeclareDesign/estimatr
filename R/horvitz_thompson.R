#' Horvitz-Thompson Estimator with Inverse Probability Weighting
#'
#' Estimates treatment effects via inverse probability weighting when
#' treatment assignment probabilities are known.  Supports all
#' \pkg{randomizr} designs as well as arbitrary designs when the user
#' supplies per-unit condition probabilities.
#'
#' @param formula A formula \code{Y ~ Z}.
#' @param data A \code{data.frame}.
#' @param condition_prs Treatment probability specification. One of:
#'   (a) an \code{ra_declaration} from \pkg{randomizr} (preferred; enables
#'       design-aware variance estimation for all standard designs),
#'   (b) a named numeric vector of condition probabilities,
#'       e.g. \code{c("0" = 0.4, "1" = 0.6)}, or
#'   (c) a two-column matrix or data frame of per-unit probabilities with
#'       columns named by condition labels.
#'   Options (b) and (c) use the conservative Young's simple-randomization
#'   variance (valid for any design; exact for Bernoulli randomization).
#' @param condition1 Label of the control condition (first sorted condition
#'   by default).
#' @param condition2 Label of the treatment condition (second sorted
#'   condition by default).
#' @param se_type \code{"youngs"} (default) or \code{"none"}.
#' @param ci Logical; whether to compute p-values and confidence intervals.
#' @param alpha Significance level, 0.05 by default.
#'
#' @return An object of class \code{"horvitz_thompson"} with fields
#'   \code{coefficients}, \code{std.error}, \code{statistic},
#'   \code{p.value}, \code{conf.low}, \code{conf.high}, \code{df},
#'   \code{nobs}, \code{vcov}, \code{se_type}, \code{condition1},
#'   \code{condition2}, \code{outcome}, and \code{term}.
#'
#' @export
horvitz_thompson <- function(formula,
                              data,
                              condition_prs = NULL,
                              condition1 = NULL,
                              condition2 = NULL,
                              se_type = "youngs",
                              ci = TRUE,
                              alpha = 0.05) {

  se_type <- match.arg(se_type, c("youngs", "none"))

  # ---- parse formula and data ----

  cpr_q  <- rlang::enquo(condition_prs)
  data_q <- rlang::enquo(data)
  cpr    <- if (rlang::quo_is_missing(cpr_q)) NULL else rlang::eval_tidy(cpr_q)
  data_df <- if (rlang::quo_is_missing(data_q)) NULL else rlang::eval_tidy(data_q)

  m_formula <- rlang::eval_tidy(rlang::enquo(formula))
  stopifnot("`formula` must be a formula" = inherits(m_formula, "formula"))
  mf <- stats::model.frame(m_formula, data = data_df, na.action = stats::na.omit)
  Y  <- as.numeric(stats::model.response(mf))
  Z  <- mf[[all.vars(m_formula[[3L]])[1L]]]

  # ---- resolve conditions ----

  conds <- if (is.factor(Z)) levels(Z) else sort(unique(Z))
  if (is.null(condition1)) condition1 <- conds[1L]
  if (is.null(condition2)) condition2 <- conds[2L]

  in12 <- Z == condition1 | Z == condition2
  Y <- Y[in12]; Z <- Z[in12]
  t2 <- which(Z == condition2); t1 <- which(Z == condition1)
  N_eff <- length(t2) + length(t1)
  if (length(t2) == 0L || length(t1) == 0L)
    stop("No units observed in one of the two conditions.")

  # ---- resolve probabilities and design ----

  prs <- ht_prs(cpr, Z, condition1, condition2, in12, data_df)

  # ---- IPW outcomes and point estimate ----

  Y2 <- Y[t2] / prs$pi2[t2]
  Y1 <- Y[t1] / prs$pi1[t1]
  estimate <- (sum(Y2) - sum(Y1)) / N_eff

  # ---- variance ----

  std.error <- NA_real_
  if (se_type == "youngs") {
    var_est <- ht_youngs(Y2, Y1, N_eff, prs, t2, t1)
    if (isTRUE(var_est >= 0)) std.error <- sqrt(var_est)
  }

  # ---- return ----

  term <- as.character(condition2)
  vcov_mat <- matrix(std.error^2, 1L, 1L, dimnames = list(term, term))
  stat <- estimate / std.error

  ret <- list(
    coefficients = setNames(estimate, term),
    std.error    = setNames(std.error, term),
    statistic    = setNames(if (ci) stat else NA_real_, term),
    p.value      = setNames(if (ci) 2 * stats::pt(abs(stat), Inf, lower.tail = FALSE)
                            else NA_real_, term),
    conf.low     = if (ci) estimate - stats::qnorm(1 - alpha / 2) * std.error else NA_real_,
    conf.high    = if (ci) estimate + stats::qnorm(1 - alpha / 2) * std.error else NA_real_,
    df           = setNames(NA_real_, term),
    alpha        = alpha,
    nobs         = N_eff,
    outcome      = deparse(m_formula[[2L]], nlines = 1L),
    term         = term,
    condition1   = condition1,
    condition2   = condition2,
    se_type      = se_type,
    vcov         = vcov_mat,
    call         = match.call()
  )
  class(ret) <- "horvitz_thompson"
  ret
}

# ---- probability resolver ----

# Returns list(pi1, pi2, design, ...) where pi1/pi2 are per-unit vectors
# for the in-study (condition1 or condition2) units, and design carries
# whatever is needed for variance computation.
ht_prs <- function(cpr, Z, condition1, condition2, in12, data_df) {
  n <- length(Z)  # in-study unit count

  if (inherits(cpr, "ra_declaration")) {
    pm     <- cpr$probabilities_matrix  # N_total × n_conditions
    cnames <- cpr$conditions
    # ra_custom may have empty conditions; fall back to prob_mat column names
    if (length(cnames) == 0L)
      cnames <- sub("^prob_", "", colnames(pm))
    c1_col <- match(as.character(condition1), cnames)
    c2_col <- match(as.character(condition2), cnames)
    if (is.na(c1_col) || is.na(c2_col))
      stop("condition1/condition2 (", condition1, ", ", condition2,
           ") not found in ra_declaration conditions: ",
           paste(cnames, collapse = ", "), ".")
    pi1_full <- pm[in12, c1_col]
    pi2_full <- pm[in12, c2_col]
    return(ht_prs_declaration(cpr, pi1_full, pi2_full, in12))
  }

  # Named scalar vector: e.g. c("0"=0.4, "1"=0.6)
  if (is.numeric(cpr) && !is.null(names(cpr)) && length(cpr) <= 10L) {
    c1_pr <- cpr[as.character(condition1)]
    c2_pr <- cpr[as.character(condition2)]
    if (any(is.na(c(c1_pr, c2_pr))))
      stop("condition1/condition2 not found in condition_prs names.")
    return(list(pi1 = rep_len(c1_pr, n), pi2 = rep_len(c2_pr, n),
                design = "simple"))
  }

  # Per-unit vector (length == n): treatment probability for each unit
  if (is.numeric(cpr) && length(cpr) == n) {
    return(list(pi1 = 1 - cpr, pi2 = cpr, design = "simple"))
  }

  # Two-column matrix or data frame
  if ((is.matrix(cpr) || is.data.frame(cpr)) && ncol(cpr) == 2L) {
    cnames <- colnames(cpr)
    c1_col <- match(as.character(condition1), cnames)
    c2_col <- match(as.character(condition2), cnames)
    if (is.na(c1_col) || is.na(c2_col))
      stop("condition1/condition2 not found in condition_prs column names.")
    return(list(pi1 = cpr[, c1_col], pi2 = cpr[, c2_col], design = "simple"))
  }

  if (is.null(cpr))
    stop("Must supply `condition_prs`: an ra_declaration, named probability ",
         "vector, or per-unit probability matrix.")
  stop("Unrecognised `condition_prs` format.")
}

ht_prs_declaration <- function(decl, pi1, pi2, in12) {
  N_total <- nrow(decl$probabilities_matrix)

  if (inherits(decl, "ra_simple")) {
    return(list(pi1 = pi1, pi2 = pi2, design = "simple"))
  }

  if (inherits(decl, "ra_custom")) {
    perm  <- decl$permutation_matrix          # N_total × K binary
    K     <- ncol(perm)
    # Full 2N × 2N joint prob matrix (control rows 1:N, treated rows N+1:2N)
    full  <- tcrossprod(rbind(1 - perm, perm)) / K
    # Subset to the in-study units
    idx   <- which(in12)
    n     <- length(idx)
    sub   <- full[c(idx, N_total + idx), c(idx, N_total + idx)]
    return(list(pi1 = pi1, pi2 = pi2, design = "custom",
                joint_mat = sub, n = n))
  }

  if (inherits(decl, "ra_blocked_and_clustered")) {
    blocks   <- decl$blocks[in12]
    clusters <- decl$clusters[in12]
    block_info <- ht_block_info(decl, in12, has_clusters = TRUE)
    return(list(pi1 = pi1, pi2 = pi2, design = "blocked_clustered",
                blocks = blocks, clusters = clusters,
                block_info = block_info))
  }

  if (inherits(decl, "ra_blocked")) {
    blocks <- decl$blocks[in12]
    block_info <- ht_block_info(decl, in12, has_clusters = FALSE)
    return(list(pi1 = pi1, pi2 = pi2, design = "blocked",
                blocks = blocks, block_info = block_info))
  }

  if (inherits(decl, "ra_clustered")) {
    clusters  <- decl$clusters[in12]
    is_simple <- isTRUE(decl$simple)
    pi2_cl    <- decl$probabilities_matrix[!duplicated(decl$clusters), 2L]
    K         <- length(pi2_cl)
    pi2_k     <- pi2_cl[1L]     # same for all clusters in complete design
    pi1_k     <- 1 - pi2_k
    return(list(pi1 = pi1, pi2 = pi2, design = "clustered",
                clusters = clusters, K = K,
                pi2_k = pi2_k, pi1_k = pi1_k,
                is_simple = is_simple))
  }

  # ra_complete (plain complete randomization)
  list(pi1 = pi1, pi2 = pi2, design = "complete",
       N_total = N_total, pi2_b = pi2[1L], pi1_b = pi1[1L])
}

# Build per-block (N_b, n2_b, n1_b) table from an ra_declaration.
ht_block_info <- function(decl, in12, has_clusters) {
  pm       <- decl$probabilities_matrix     # N_total × 2
  blk_full <- decl$blocks
  blk_tbl  <- table(blk_full)

  if (has_clusters) {
    uniq_cl  <- !duplicated(blk_full)
    cl_full  <- decl$clusters
    uniq_cl2 <- !duplicated(cl_full)
    cl_blks  <- blk_full[uniq_cl2]
    cl_pi2   <- pm[uniq_cl2, 2L]
    cl_tbl   <- table(cl_blks)
    block_levels <- names(blk_tbl)
    info <- lapply(block_levels, function(b) {
      K_b   <- as.integer(cl_tbl[b])
      pi2_b <- cl_pi2[cl_blks == b][1L]
      N_b   <- as.integer(blk_tbl[b])
      list(N_b = N_b, K_b = K_b, pi2_b = pi2_b, pi1_b = 1 - pi2_b)
    })
  } else {
    block_levels <- names(blk_tbl)
    info <- lapply(block_levels, function(b) {
      N_b   <- as.integer(blk_tbl[b])
      pi2_b <- pm[blk_full == b, 2L][1L]
      list(N_b = N_b, pi2_b = pi2_b, pi1_b = 1 - pi2_b)
    })
  }
  names(info) <- block_levels
  info
}

# ---- joint probability helpers ----

# Pure-R equivalent of estimatr's gen_joint_pr_complete.
# Handles non-integer m = pi * K by mixing floor/ceil probabilities.
gen_jpr_complete <- function(pi, K) {
  m   <- pi * K
  m_f <- floor(m)
  r   <- m %% 1
  nc  <- K - m_f
  pi11 <- r * ((m_f + 1) / K) * (m_f / (K - 1)) +
          (1 - r) * (m_f / K) * (max(m_f - 1, 0) / (K - 1))
  pi00 <- r * (max(nc - 1, 0) / K) * (max(nc - 2, 0) / (K - 1)) +
          (1 - r) * (nc / K) * (max(nc - 1, 0) / (K - 1))
  pi10 <- r * (max(nc - 1, 0) / K) * ((m_f + 1) / (K - 1)) +
          (1 - r) * (nc / K) * (m_f / (K - 1))
  list(pi11 = pi11, pi00 = pi00, pi10 = pi10)
}

# Young's coefficient for the treated-treated block given marginal pi and
# joint pi11. Works for non-integer m (general form of -(K-m)/(K*(m-1))).
ht_c11 <- function(pi2, pi11) if (pi11 > 0) 1 - pi2^2 / pi11 else 0
ht_c00 <- function(pi1, pi00) if (pi00 > 0) 1 - pi1^2 / pi00 else 0
ht_crs <- function(pi2, pi11) if (pi11 < pi2) (pi11 - pi2^2) / (pi11 - pi2) else 0

# ---- Young's inequality variance ----

ht_youngs <- function(Y2, Y1, N_eff, prs, t2, t1) {
  switch(prs$design,
    simple           = ht_var_simple(Y2, Y1, N_eff),
    complete         = ht_var_complete(Y2, Y1, N_eff,
                                       prs$N_total, prs$pi2_b, prs$pi1_b),
    blocked          = ht_var_blocked(Y2, Y1, N_eff, t2, t1,
                                      prs$blocks, prs$block_info),
    clustered        = ht_var_clustered(Y2, Y1, N_eff, t2, t1,
                                        prs$clusters, prs$K,
                                        prs$pi2_k, prs$pi1_k, prs$is_simple),
    blocked_clustered = ht_var_blocked_clustered(Y2, Y1, N_eff, t2, t1,
                                                  prs$blocks, prs$clusters,
                                                  prs$block_info),
    custom           = ht_var_custom(Y2, Y1, N_eff, prs$joint_mat,
                                     prs$pi2[t2], prs$pi1[t1], t2, t1,
                                     prs$n)
  )
}

# Simple (Bernoulli) / conservative bound: no joint probabilities needed.
ht_var_simple <- function(Y2, Y1, N_eff) {
  (sum(Y2^2) + sum(Y1^2)) / N_eff^2
}

# Complete randomization: O(1) scalar formula (no matrix needed).
# Uses gen_jpr_complete to handle non-integer m (pi * K not an integer),
# matching estimatr's weighted-average approach exactly.
ht_var_complete <- function(Y2, Y1, N_eff, K, pi2, pi1) {
  jpr  <- gen_jpr_complete(pi2, K)
  c11  <- ht_c11(pi2, jpr$pi11)
  c00  <- ht_c00(pi1, jpr$pi00)
  crs  <- ht_crs(pi2, jpr$pi11)
  S2 <- sum(Y2); Q2 <- sum(Y2^2)
  S1 <- sum(Y1); Q1 <- sum(Y1^2)
  (Q2 + Q1 + c11 * (S2^2 - Q2) + c00 * (S1^2 - Q1) - 2 * crs * S2 * S1) / N_eff^2
}

# Blocked complete randomization: sum within-block contributions.
ht_var_blocked <- function(Y2, Y1, N_eff, t2, t1, blocks, block_info) {
  b_t2 <- blocks[t2]; b_t1 <- blocks[t1]
  varN2 <- 0
  for (b in names(block_info)) {
    bi <- block_info[[b]]
    Y2_b <- Y2[b_t2 == b]; Y1_b <- Y1[b_t1 == b]
    if (length(Y2_b) == 0L && length(Y1_b) == 0L) next
    jpr  <- gen_jpr_complete(bi$pi2_b, bi$N_b)
    c11  <- ht_c11(bi$pi2_b, jpr$pi11)
    c00  <- ht_c00(bi$pi1_b, jpr$pi00)
    crs  <- ht_crs(bi$pi2_b, jpr$pi11)
    S2 <- sum(Y2_b); Q2 <- sum(Y2_b^2)
    S1 <- sum(Y1_b); Q1 <- sum(Y1_b^2)
    varN2 <- varN2 + Q2 + Q1 +
             c11 * (S2^2 - Q2) + c00 * (S1^2 - Q1) - 2 * crs * S2 * S1
  }
  varN2 / N_eff^2
}

# Clustered: aggregate to cluster-level IPW sums, then apply complete formula
# at cluster level (or simple formula for simple cluster randomization).
ht_var_clustered <- function(Y2, Y1, N_eff, t2, t1, clusters, K, pi2_k, pi1_k, is_simple) {
  cl_t2 <- clusters[t2]; cl_t1 <- clusters[t1]
  Y2_cl <- as.numeric(tapply(Y2, cl_t2, sum))
  Y1_cl <- as.numeric(tapply(Y1, cl_t1, sum))
  if (is_simple)
    return(ht_var_simple(Y2_cl, Y1_cl, N_eff))
  ht_var_complete(Y2_cl, Y1_cl, N_eff, K, pi2_k, pi1_k)
}

# Blocked + clustered: aggregate to cluster level within each block,
# sum block contributions.
ht_var_blocked_clustered <- function(Y2, Y1, N_eff, t2, t1,
                                      blocks, clusters, block_info) {
  b_t2 <- blocks[t2];   cl_t2 <- clusters[t2]
  b_t1 <- blocks[t1];   cl_t1 <- clusters[t1]
  varN2 <- 0
  for (b in names(block_info)) {
    bi <- block_info[[b]]
    Y2_b <- Y2[b_t2 == b]; cl_t2_b <- cl_t2[b_t2 == b]
    Y1_b <- Y1[b_t1 == b]; cl_t1_b <- cl_t1[b_t1 == b]
    if (length(Y2_b) == 0L && length(Y1_b) == 0L) next
    Y2_cl <- if (length(Y2_b)) as.numeric(tapply(Y2_b, cl_t2_b, sum)) else numeric(0)
    Y1_cl <- if (length(Y1_b)) as.numeric(tapply(Y1_b, cl_t1_b, sum)) else numeric(0)
    jpr  <- gen_jpr_complete(bi$pi2_b, bi$K_b)
    c11  <- ht_c11(bi$pi2_b, jpr$pi11)
    c00  <- ht_c00(bi$pi1_b, jpr$pi00)
    crs  <- ht_crs(bi$pi2_b, jpr$pi11)
    S2 <- sum(Y2_cl); Q2 <- sum(Y2_cl^2)
    S1 <- sum(Y1_cl); Q1 <- sum(Y1_cl^2)
    varN2 <- varN2 + Q2 + Q1 +
             c11 * (S2^2 - Q2) + c00 * (S1^2 - Q1) - 2 * crs * S2 * S1
  }
  varN2 / N_eff^2
}

# Custom design: use the joint probability submatrix (from tcrossprod).
# O(n^2) but feasible since ra_custom requires an explicit permutation matrix.
#
# The Young's coefficient for a pair (i,j) is A[i,j] = 1 - pi_i*pi_j / p_ij,
# not (p_ij - pi_i*pi_j).  This is the correct formula (verified against the
# C++ ht_var_partial and ht_covar_partial for complete randomization designs).
ht_var_custom <- function(Y2, Y1, N_eff, joint_mat, pi2, pi1, t2, t1, n) {
  # joint_mat is 2n×2n: rows/cols 1:n = condition1, n+1:2n = condition2
  p22 <- joint_mat[n + t2, n + t2, drop = FALSE]
  p11 <- joint_mat[t1,     t1,     drop = FALSE]
  p21 <- joint_mat[n + t2, t1,     drop = FALSE]

  # Off-diagonal correction matrices: A[i,j] = 1 - pi_i*pi_j / p_ij
  # Diagonal is zero (self-variance already in sum(Y^2)).
  B22 <- 1 - outer(pi2, pi2) / p22; diag(B22) <- 0
  B11 <- 1 - outer(pi1, pi1) / p11; diag(B11) <- 0
  B21 <- 1 - outer(pi2, pi1) / p21   # no diagonal: i and j are different units

  vp22 <- sum(outer(Y2, Y2) * B22)
  vp11 <- sum(outer(Y1, Y1) * B11)
  vc   <- sum(outer(Y2, Y1) * B21)

  (sum(Y2^2) + sum(Y1^2) + vp22 + vp11 - 2 * vc) / N_eff^2
}
