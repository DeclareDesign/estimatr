#' Horvitz-Thompson Estimator with Inverse Probability Weighting
#'
#' Estimates treatment effects via inverse probability weighting when
#' treatment assignment probabilities are known.  Supports all
#' \pkg{randomizr} designs as well as arbitrary designs supplied via a
#' permutation matrix.
#'
#' @details With more than two arms, \code{condition1} and \code{condition2}
#'   select the contrast, and the estimand remains the average treatment
#'   effect over all N units the design covers.  The estimator therefore
#'   divides by N, not by the number of units landing in the two conditions,
#'   and the variance uses the joint assignment probabilities implied by the
#'   arm sizes.  \code{data} must hold one row per unit of the design, in the
#'   design's order, including units assigned to arms outside the contrast.
#'
#' @param formula A formula \code{Y ~ Z}.
#' @param data A \code{data.frame} with one row per unit of the design.
#' @param condition_prs Treatment probability specification. One of:
#'   \itemize{
#'     \item An \code{ra_declaration} from \pkg{randomizr} — \strong{strongly
#'       preferred.}  All standard designs (simple/Bernoulli, complete,
#'       blocked, clustered, blocked-and-clustered, and arbitrary permutation
#'       matrices) are supported, and the variance estimator uses exact
#'       design-aware joint inclusion probabilities.  Any design for which
#'       you know the block structure, cluster structure, marginal treatment
#'       probabilities, and whether randomization is simple or complete can
#'       be expressed as \code{declare_ra(blocks = bl, clusters = cl,
#'       prob = pi, simple = FALSE)} — there is no parametric design that
#'       requires the alternatives below.  For fully custom designs, use
#'       \code{declare_ra(permutation_matrix = perm)}.
#'     \item A named numeric vector of marginal condition probabilities,
#'       e.g. \code{c("0" = 0.4, "1" = 0.6)}.  Uses the conservative
#'       Young's simple-randomization variance bound, which is valid for
#'       any design but exact only for Bernoulli (simple) randomization.
#'       For complete or blocked designs this overstates uncertainty; use
#'       an \code{ra_declaration} to get the tighter design-aware variance.
#'     \item A two-column matrix or data frame of per-unit probabilities with
#'       columns named by condition labels.  Same conservative variance as (b).
#'   }
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
#'   \code{nobs} (the number of units in the design, including any arms
#'   outside the contrast), \code{vcov}, \code{se_type}, \code{condition1},
#'   \code{condition2}, \code{outcome}, and \code{term}.
#'
#' @examples
#' set.seed(40)
#' dat <- data.frame(y = rnorm(100), z = rep(0:1, 50))
#'
#' # A named vector of condition probabilities gives the conservative
#' # simple-randomization bound, valid for any design
#' horvitz_thompson(y ~ z, data = dat, condition_prs = c("0" = 0.5, "1" = 0.5))
#'
#' # Passing the randomization declaration instead is what buys the
#' # design-aware variance, and it is the recommended form
#' if (requireNamespace("randomizr", quietly = TRUE)) {
#'   decl <- randomizr::declare_ra(N = 100, m = 50)
#'   dat$z2 <- randomizr::conduct_ra(decl)
#'   print(horvitz_thompson(y ~ z2, data = dat, condition_prs = decl))
#'
#'   # Blocked and clustered designs need no extra arguments: the declaration
#'   # already carries the structure
#'   bl <- rep(1:4, each = 25)
#'   decl_bl <- randomizr::declare_ra(blocks = bl, prob = 0.4)
#'   dat$z3 <- randomizr::conduct_ra(decl_bl)
#'   print(horvitz_thompson(y ~ z3, data = dat, condition_prs = decl_bl))
#'
#'   # Any two arms of a multi-arm design can be contrasted, with the estimand
#'   # still defined over all N units
#'   decl3 <- randomizr::declare_ra(N = 100, conditions = c("a", "b", "c"))
#'   dat$z4 <- randomizr::conduct_ra(decl3)
#'   print(horvitz_thompson(y ~ z4, data = dat, condition_prs = decl3,
#'                          condition1 = "a", condition2 = "c"))
#' }
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

  # Track which original-data rows survived na.omit so we can index the
  # ra_declaration's probability matrix correctly when rows are dropped.
  na_dropped <- attr(mf, "na.action")
  if (!is.null(data_df) && !is.null(na_dropped)) {
    surv_idx <- setdiff(seq_len(nrow(data_df)), as.integer(na_dropped))
  } else {
    surv_idx <- seq_len(nrow(mf))
  }

  if (inherits(cpr, "ra_declaration") && !is.null(data_df) &&
      nrow(cpr$probabilities_matrix) != nrow(data_df))
    stop("`condition_prs` declares ", nrow(cpr$probabilities_matrix),
         " units but `data` has ", nrow(data_df), " rows. The declaration ",
         "must cover exactly the rows of `data`, in the same order, so that ",
         "each unit is matched to its own assignment probabilities.")

  # The estimand is the average effect over every unit the assignment
  # probabilities refer to, so N is the whole analysis sample and not just
  # the units landing in the two conditions being contrasted. The two
  # coincide when the design has two arms.
  N <- nrow(mf)

  # ---- resolve conditions ----

  conds <- if (is.factor(Z)) levels(Z) else sort(unique(Z))
  if (is.null(condition1)) condition1 <- conds[1L]
  if (is.null(condition2)) condition2 <- conds[2L]

  in12 <- Z == condition1 | Z == condition2
  Y <- Y[in12]; Z <- Z[in12]
  t2 <- which(Z == condition2); t1 <- which(Z == condition1)
  if (length(t2) == 0L || length(t1) == 0L)
    stop("No units observed in one of the two conditions.")

  # Row indices in the original data for the in-study (condition1/2) units
  row_idx <- surv_idx[in12]

  # ---- resolve probabilities and design ----

  prs <- ht_prs(cpr, Z, condition1, condition2, row_idx, data_df)

  # ---- IPW outcomes and point estimate ----

  Y2 <- Y[t2] / prs$pi2[t2]
  Y1 <- Y[t1] / prs$pi1[t1]
  estimate <- (sum(Y2) - sum(Y1)) / N

  # ---- variance ----

  std.error <- NA_real_
  if (se_type == "youngs") {
    var_est <- ht_youngs(Y2, Y1, N, prs, t2, t1)
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
    nobs         = N,
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
# row_idx: integer vector of original-data row indices for the in-study units.
# Using original row indices (not model-frame logical) ensures correct
# probability lookup from pm even when rows were dropped via na.omit.
ht_prs <- function(cpr, Z, condition1, condition2, row_idx, data_df) {
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
    pi1_full <- pm[row_idx, c1_col]
    pi2_full <- pm[row_idx, c2_col]
    return(ht_prs_declaration(cpr, pi1_full, pi2_full, row_idx,
                              c1_col, c2_col, condition1, condition2))
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

ht_prs_declaration <- function(decl, pi1, pi2, row_idx,
                               c1_col, c2_col, condition1, condition2) {
  N_total <- nrow(decl$probabilities_matrix)

  if (inherits(decl, "ra_simple")) {
    return(list(pi1 = pi1, pi2 = pi2, design = "simple"))
  }

  if (inherits(decl, "ra_custom")) {
    perm  <- decl$permutation_matrix          # N_total × K of condition labels
    K     <- ncol(perm)
    ind1  <- (perm == condition1) * 1
    ind2  <- (perm == condition2) * 1
    # Full 2N × 2N joint prob matrix (condition1 rows 1:N, condition2 rows N+1:2N)
    full  <- tcrossprod(rbind(ind1, ind2)) / K
    # Subset to the in-study units using original-data row indices
    idx   <- row_idx
    n     <- length(idx)
    sub   <- full[c(idx, N_total + idx), c(idx, N_total + idx)]
    return(list(pi1 = pi1, pi2 = pi2, design = "custom",
                joint_mat = sub, n = n))
  }

  if (inherits(decl, "ra_blocked_and_clustered")) {
    blocks   <- decl$blocks[row_idx]
    clusters <- decl$clusters[row_idx]
    block_info <- ht_block_info(decl, TRUE, c1_col, c2_col)
    return(list(pi1 = pi1, pi2 = pi2, design = "blocked_clustered",
                blocks = blocks, clusters = clusters,
                block_info = block_info))
  }

  if (inherits(decl, "ra_blocked")) {
    blocks <- decl$blocks[row_idx]
    block_info <- ht_block_info(decl, FALSE, c1_col, c2_col)
    return(list(pi1 = pi1, pi2 = pi2, design = "blocked",
                blocks = blocks, block_info = block_info))
  }

  if (inherits(decl, "ra_clustered")) {
    clusters  <- decl$clusters[row_idx]
    is_simple <- isTRUE(decl$simple)
    uniq_cl   <- !duplicated(decl$clusters)
    K         <- sum(uniq_cl)
    # Same for every cluster under complete cluster randomization
    pi2_k     <- decl$probabilities_matrix[uniq_cl, c2_col][1L]
    pi1_k     <- decl$probabilities_matrix[uniq_cl, c1_col][1L]
    return(list(pi1 = pi1, pi2 = pi2, design = "clustered",
                clusters = clusters, K = K,
                pi2_k = pi2_k, pi1_k = pi1_k,
                is_simple = is_simple))
  }

  # ra_complete (plain complete randomization)
  list(pi1 = pi1, pi2 = pi2, design = "complete",
       N_total = N_total, pi2_b = pi2[1L], pi1_b = pi1[1L])
}

# Build per-block design info from an ra_declaration (full design, not subset).
# Both condition columns are read from the probability matrix rather than
# taking pi1 = 1 - pi2, which only holds when the design has two arms.
ht_block_info <- function(decl, has_clusters, c1_col, c2_col) {
  pm       <- decl$probabilities_matrix     # N_total × n_conditions
  blk_full <- decl$blocks
  blk_tbl  <- table(blk_full)
  block_levels <- names(blk_tbl)

  if (has_clusters) {
    uniq_cl <- !duplicated(decl$clusters)
    cl_blks <- blk_full[uniq_cl]
    cl_pi2  <- pm[uniq_cl, c2_col]
    cl_pi1  <- pm[uniq_cl, c1_col]
    cl_tbl  <- table(cl_blks)
    info <- lapply(block_levels, function(b) {
      in_b <- cl_blks == b
      list(N_b = as.integer(blk_tbl[b]), K_b = as.integer(cl_tbl[b]),
           pi2_b = cl_pi2[in_b][1L], pi1_b = cl_pi1[in_b][1L])
    })
  } else {
    info <- lapply(block_levels, function(b) {
      in_b <- blk_full == b
      list(N_b = as.integer(blk_tbl[b]),
           pi2_b = pm[in_b, c2_col][1L], pi1_b = pm[in_b, c1_col][1L])
    })
  }
  names(info) <- block_levels
  info
}

# ---- joint probability helpers ----

# Pure-R equivalent of estimatr's gen_joint_pr_complete, for two-arm designs.
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
  list(pi11 = pi11, pi00 = pi00)
}

# Young's coefficient for the treated-treated block given marginal pi and
# joint pi11. Works for non-integer m (general form of -(K-m)/(K*(m-1))).
ht_c11 <- function(pi2, pi11) if (pi11 > 0) 1 - pi2^2 / pi11 else 0
ht_c00 <- function(pi1, pi00) if (pi00 > 0) 1 - pi1^2 / pi00 else 0
ht_crs <- function(pi2, pi11) if (pi11 < pi2) (pi11 - pi2^2) / (pi11 - pi2) else 0

# Young's coefficients A_ij = 1 - pi_i pi_j / pi_ij for the three kinds of
# pair (both in condition2, both in condition1, one of each) under complete
# randomization of n units with m2 = pi2*n in condition2 and m1 = pi1*n in
# condition1. Exchangeability gives
#   pi_22 = m2(m2-1)/(n(n-1)), pi_11 = m1(m1-1)/(n(n-1)), pi_21 = m2 m1/(n(n-1))
# so the cross coefficient collapses to 1/n for any complete design.
#
# With two arms m1 = n - m2, and the design may imply a non-integer m2, in
# which case the realized count is floor(m2) or floor(m2)+1. gen_jpr_complete
# averages over that mixture, so two-arm contrasts keep that path. With more
# than two arms m1 is a free quantity and the closed forms below apply.
ht_coefs <- function(pi2, pi1, n) {
  if (abs(pi2 + pi1 - 1) < 1e-10) {
    jpr <- gen_jpr_complete(pi2, n)
    return(list(c22 = ht_c11(pi2, jpr$pi11),
                c11 = ht_c00(pi1, jpr$pi00),
                c21 = ht_crs(pi2, jpr$pi11)))
  }
  m2 <- pi2 * n
  m1 <- pi1 * n
  list(
    c22 = if (m2 > 1) 1 - m2 * (n - 1) / (n * (m2 - 1)) else 0,
    c11 = if (m1 > 1) 1 - m1 * (n - 1) / (n * (m1 - 1)) else 0,
    c21 = 1 / n
  )
}

# Young's bound contribution from one completely randomized group of n units
# (a whole sample, a block, or a set of clusters). Y2_g and Y1_g are the IPW
# outcomes of the group's units in each condition. Returns the numerator; the
# caller divides by N^2 once.
ht_group_contrib <- function(Y2_g, Y1_g, pi2, pi1, n) {
  cf <- ht_coefs(pi2, pi1, n)
  S2 <- sum(Y2_g); Q2 <- sum(Y2_g^2)
  S1 <- sum(Y1_g); Q1 <- sum(Y1_g^2)
  Q2 + Q1 + cf$c22 * (S2^2 - Q2) + cf$c11 * (S1^2 - Q1) - 2 * cf$c21 * S2 * S1
}

# ---- Young's inequality variance ----

ht_youngs <- function(Y2, Y1, N, prs, t2, t1) {
  switch(prs$design,
    simple           = ht_var_simple(Y2, Y1, N),
    complete         = ht_var_complete(Y2, Y1, N,
                                       prs$N_total, prs$pi2_b, prs$pi1_b),
    blocked          = ht_var_blocked(Y2, Y1, N, t2, t1,
                                      prs$blocks, prs$block_info),
    clustered        = ht_var_clustered(Y2, Y1, N, t2, t1,
                                        prs$clusters, prs$K,
                                        prs$pi2_k, prs$pi1_k, prs$is_simple),
    blocked_clustered = ht_var_blocked_clustered(Y2, Y1, N, t2, t1,
                                                  prs$blocks, prs$clusters,
                                                  prs$block_info),
    custom           = ht_var_custom(Y2, Y1, N, prs$joint_mat,
                                     prs$pi2[t2], prs$pi1[t1], t2, t1,
                                     prs$n)
  )
}

# Simple (Bernoulli) / conservative bound: no joint probabilities needed.
ht_var_simple <- function(Y2, Y1, N) {
  (sum(Y2^2) + sum(Y1^2)) / N^2
}

# Complete randomization: O(1) scalar formula (no matrix needed).
ht_var_complete <- function(Y2, Y1, N, n_units, pi2, pi1) {
  ht_group_contrib(Y2, Y1, pi2, pi1, n_units) / N^2
}

# Blocked complete randomization: sum within-block contributions.
ht_var_blocked <- function(Y2, Y1, N, t2, t1, blocks, block_info) {
  b_t2 <- blocks[t2]; b_t1 <- blocks[t1]
  varN2 <- 0
  for (b in names(block_info)) {
    bi <- block_info[[b]]
    Y2_b <- Y2[b_t2 == b]; Y1_b <- Y1[b_t1 == b]
    if (length(Y2_b) == 0L && length(Y1_b) == 0L) next
    varN2 <- varN2 + ht_group_contrib(Y2_b, Y1_b, bi$pi2_b, bi$pi1_b, bi$N_b)
  }
  varN2 / N^2
}

# Clustered: aggregate to cluster-level IPW sums, then apply complete formula
# at cluster level (or simple formula for simple cluster randomization).
ht_var_clustered <- function(Y2, Y1, N, t2, t1, clusters, K, pi2_k, pi1_k, is_simple) {
  cl_t2 <- clusters[t2]; cl_t1 <- clusters[t1]
  Y2_cl <- as.numeric(tapply(Y2, cl_t2, sum))
  Y1_cl <- as.numeric(tapply(Y1, cl_t1, sum))
  if (is_simple)
    return(ht_var_simple(Y2_cl, Y1_cl, N))
  ht_var_complete(Y2_cl, Y1_cl, N, K, pi2_k, pi1_k)
}

# Blocked + clustered: aggregate to cluster level within each block,
# sum block contributions.
ht_var_blocked_clustered <- function(Y2, Y1, N, t2, t1,
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
    varN2 <- varN2 + ht_group_contrib(Y2_cl, Y1_cl, bi$pi2_b, bi$pi1_b, bi$K_b)
  }
  varN2 / N^2
}

# Custom design: use the joint probability submatrix (from tcrossprod).
# O(n^2) but feasible since ra_custom requires an explicit permutation matrix.
#
# The Young's coefficient for a pair (i,j) is A[i,j] = 1 - pi_i*pi_j / p_ij,
# not (p_ij - pi_i*pi_j).  This is the correct formula (verified against the
# C++ ht_var_partial and ht_covar_partial for complete randomization designs).
ht_var_custom <- function(Y2, Y1, N, joint_mat, pi2, pi1, t2, t1, n) {
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

  (sum(Y2^2) + sum(Y1^2) + vp22 + vp11 - 2 * vc) / N^2
}
