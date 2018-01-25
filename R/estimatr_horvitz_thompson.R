#' Horvitz-Thompson estimator for two-armed trials
#'
#' @description Horvitz-Thompson estimators that are unbiased for designs in
#' which the randomization scheme is known
#'
#' @param formula an object of class formula, as in \code{\link{lm}}, such as
#' \code{Y ~ Z} with only one variable on the right-hand side, the treatment.
#' @param data A data.frame.
#' @param condition_prs An optional bare (unquoted) name of the variable with
#' the condition 2 (treatment) probabilities. See details.
#' @param blocks An optional bare (unquoted) name of the block variable. Use
#' for blocked designs only. See details.
#' @param clusters An optional bare (unquoted) name of the variable that
#' corresponds to the clusters in the data; used for cluster randomized
#' designs. For blocked designs, clusters must be within blocks.
#' @param simple logical, optional. Whether the randomization is simple
#' (TRUE) or complete (FALSE). This is ignored if \code{blocks} are specified,
#' as all blocked designs use complete randomization, or either
#' \code{declaration} or \code{condition_pr_mat} are passed. Otherwise, it
#' defaults to \code{TRUE}.
#' @param condition_pr_mat An optional 2n * 2n matrix of marginal and joint
#' probabilities of all units in condition1 and condition2. See details.
#' @param declaration An object of class \code{"ra_declaration"}, from
#' the \code{\link[randomizr]{declare_ra}} function in the \pkg{randomizr}
#' package. This is the third way that one can specify a design for this
#' estimator. Cannot be used along with any of \code{condition_prs},
#' \code{blocks}, \code{clusters}, or \code{condition_pr_mat}. See details.
#' @param subset An optional bare (unquoted) expression specifying a subset of
#' observations to be used.
#' @param se_type can be one of \code{c("youngs", "constant")} and corresponds
#' the estimator of the standard errors. Default estimator uses Young's
#' inequality (and is conservative) while the other uses a constant treatment
#' effects assumption and only works for simple randomized designs at the
#' moment.
#' @param condition1 value in the treatment vector of the condition
#' to be the control. Effects are
#' estimated with \code{condition1} as the control and \code{condition2} as the
#' treatment. If unspecified, \code{condition1} is the "first" condition and
#' \code{condition2} is the "second" according to levels if the treatment is a
#' factor or according to a sortif it is a numeric or character variable (i.e
#' if unspecified and the treatment is 0s and 1s, \code{condition1} will by
#' default be 0 and \code{condition2} will be 1). See the examples for more.
#' @param condition2 value in the treatment vector of the condition to be the
#' treatment. See \code{condition1}.
#' @param ci logical. Whether to compute and return p-values and
#' confidence intervals, TRUE by default.
#' @param alpha The significance level, 0.05 by default.
#' @param return_condition_pr_mat logical. Whether to return the condition
#' probability matrix. Returns NULL if the design is simple randomization,
#' FALSE by default.
#'
#' @details This function implements the Horvitz-Thompson estimator for
#' treatment effects for two-armed trials. This estimator is useful for estimating unbiased
#' treatment effects given any randomization scheme as long as the
#' randomization scheme is known.
#'
#' In short, the Horvitz-Thompson estimator essentially reweights each unit
#' by the probability of it being in its observed condition. Pivotal to the
#' estimation of treatment effects using this estimator are the marginal
#' condition probabilities (i.e., the probability that any one unit is in
#' a particular treatment condition). Pivotal to the estimating the variance
#' variance whenever the design is more complicated than simple randomization,
#' are the
#' joint condition probabilities (i.e., the probabilities that any two units
#' have a particular set of treatment conditions, either the same or
#' different). The estimator we provide here considers the case with two
#' treatment conditions.
#'
#' Users interested in more details can see the
#' \href{http://estimatr.declaredesign.org/articles/technical-notes.html}{technical notes}
#' for more information and references, or see the references below.
#'
#' There are three distinct ways that users can specify the design to the
#' function. The preferred way is to use
#' the \code{\link[randomizr]{declare_ra}} function in the \pkg{randomizr}
#' package. This function takes several arguments, including blocks, clusters,
#' treatment probabilities, whether randomization is simple or not, and more.
#' Passing the outcome of that function, an object of class
#' \code{"ra_declaration"} to the \code{declaration} argument in this function,
#' will lead to a call of the \code{\link{declaration_to_condition_pr_mat}}
#' function which generates the condition probability matrix needed to
#' estimate treatment effects and standard errors. We provide many examples
#' below of how this could be done.
#'
#' The second way is to pass the names of vectors in your \code{data} to
#' \code{condition_prs}, \code{blocks}, and \code{clusters}. You can further
#' specify whether the randomization was simple or complete using the \code{simple}
#' argument. Note that if \code{blocks} are specified the randomization is
#' always treated as complete. From these vectors, the function learns how to
#' build the condition probability matrix that is used in estimation.
#'
#' In the case
#' where \code{condition_prs} is specified, this function assumes those
#' probabilities are the marginal probability that each unit is in condition2
#' and then uses the other arguments (\code{blocks}, \code{clusters},
#' \code{simple}) to learn the rest of the design. If users do not pass
#' \code{condition_prs}, this function learns the probability of being in
#' condition2 from the data. That is, none of these arguments are specified,
#' we assume that there was a simple randomization where the probability
#' of each unit being in condition2 was the average of all units in condition2.
#' Similarly, we learn the block-level probability of treatment within
#' \code{blocks} by looking at the mean number of units in condition2 if
#' \code{condition_prs} is not specified.
#'
#' The third way is to pass a \code{condition_pr_mat} directly. One can
#' see more about this object in the documentation for
#' \code{\link{declaration_to_condition_pr_mat}} and
#' \code{\link{permutations_to_condition_pr_mat}}. Essentially, this 2n * 2n
#' matrix allows users to specify marginal and joint marginal probabilities
#' of units being in conditions 1 and 2 of arbitrary complexity. Users should
#' only use this option if they are certain they know what they are doing.
#'
#' @return Returns an object of class \code{"horvitz_thompson"}.
#'
#' The post-estimation commands functions \code{summary} and \code{\link{tidy}}
#' return results in a \code{data.frame}. To get useful data out of the return,
#' you can use these data frames, you can use the resulting list directly, or
#' you can use the generic accessor functions \code{coef} and
#' \code{confint}.
#'
#' An object of class \code{"horvitz_thompson"} is a list containing at
#' least the following components:
#'
#'   \item{coefficients}{the estimated coefficients}
#'   \item{se}{the estimated standard errors}
#'   \item{df}{the estimated degrees of freedom}
#'   \item{p}{the p-values from from a two-sided z-test using \code{coefficients}, \code{se}, and \code{df}}
#'   \item{ci_lower}{the lower bound of the \code{1 - alpha} percent confidence interval}
#'   \item{ci_upper}{the upper bound of the \code{1 - alpha} percent confidence interval}
#'   \item{coefficient_name}{a character vector of coefficient names}
#'   \item{alpha}{the significance level specified by the user}
#'   \item{N}{the number of observations used}
#'   \item{outcome}{the name of the outcome variable}
#'   \item{condition_pr_mat}{the condition probability matrix if \code{return_condition_pr_mat} is TRUE}
#'
#'
#' @seealso \code{\link[randomizr]{declare_ra}}
#'
#' @references
#' Aronow, Peter M, and Joel A Middleton. 2013. "A Class of Unbiased Estimators of the Average Treatment Effect in Randomized Experiments." Journal of Causal Inference 1 (1): 135-54. \url{https://doi.org/10.1515/jci-2012-0009}.
#'
#' Aronow, Peter M, and Cyrus Samii. 2017. "Estimating Average Causal Effects Under Interference Between Units." Annals of Applied Statistics, forthcoming. \url{https://arxiv.org/abs/1305.6156v3}.
#'
#' Middleton, Joel A, and Peter M Aronow. 2015. "Unbiased Estimation of the Average Treatment Effect in Cluster-Randomized Experiments." Statistics, Politics and Policy 6 (1-2): 39-75. \url{https://doi.org/10.1515/spp-2013-0002}.
#'
#' @examples
#'
#' # Set seed
#' set.seed(42)
#'
#' # Simulate data
#' n <- 10
#' dat <- data.frame(y = rnorm(n))
#'
#' library(randomizr)
#'
#' #----------
#' # 1. Simple random assignment
#' #----------
#' dat$p <- 0.5
#' dat$z <- rbinom(n, size = 1, prob = dat$p)
#'
#' # If you only pass condition_prs, we assume simple random sampling
#' horvitz_thompson(y ~ z, data = dat, condition_prs = p)
#' # Assume constant effects instead
#' horvitz_thompson(y ~ z, data = dat, condition_prs = p, se_type = "constant")
#'
#' # Also can use randomizr to pass a declaration
#' srs_declaration <- declare_ra(N = nrow(dat), prob = 0.5, simple = TRUE)
#' horvitz_thompson(y ~ z, data = dat, declaration = srs_declaration)
#'
#' #----------
#' # 2. Complete random assignment
#' #----------
#'
#' dat$z <- sample(rep(0:1, each = n/2))
#' # Can use a declaration
#' crs_declaration <- declare_ra(N = nrow(dat), prob = 0.5, simple = FALSE)
#' horvitz_thompson(y ~ z, data = dat, declaration = crs_declaration)
#' # Can precompute condition_pr_mat and pass it
#' # (faster for multiple runs with same condition probability matrix)
#' crs_pr_mat <- declaration_to_condition_pr_mat(crs_declaration)
#' horvitz_thompson(y ~ z, data = dat, condition_pr_mat = crs_pr_mat)
#'
#' #----------
#' # 3. Clustered treatment, complete random assigment
#' #-----------
#' # Simulating data
#' dat$cl <- rep(1:4, times = c(2, 2, 3, 3))
#' dat$prob <- 0.5
#' clust_crs_decl <- declare_ra(N = nrow(dat), clusters = dat$cl, prob = 0.5)
#' dat$z <- conduct_ra(clust_crs_decl)
#' # Easiest to specify using declaration
#' ht_cl <- horvitz_thompson(y ~ z, data = dat, declaration = clust_crs_decl)
#' # Also can pass the condition probability and the clusters
#' ht_cl_manual <- horvitz_thompson(
#'   y ~ z,
#'   data = dat,
#'   clusters = cl,
#'   condition_prs = prob,
#'   simple = FALSE
#' )
#' ht_cl
#' ht_cl_manual
#'
#' # Blcoked estimators specified similarly
#'
#' #----------
#' # More complicated assignment
#' #----------
#'
#' # arbitrary permutation matrix
#' possible_treats <- cbind(
#'   c(1, 1, 0, 1, 0, 0, 0, 1, 1, 0),
#'   c(0, 1, 1, 0, 1, 1, 0, 1, 0, 1),
#'   c(1, 0, 1, 1, 1, 1, 1, 0, 0, 0)
#' )
#' arb_pr_mat <- permutations_to_condition_pr_mat(possible_treats)
#' # Simulating a column to be realized treatment
#' dat$z <- possible_treats[, sample(ncol(possible_treats), size = 1)]
#' horvitz_thompson(y ~ z, data = dat, condition_pr_mat = arb_pr_mat)
#'
#' @export
horvitz_thompson <-
  function(formula,
           data,
           blocks,
           clusters,
           simple = NULL,
           condition_prs,
           condition_pr_mat = NULL,
           declaration = NULL,
           subset,
           condition1 = NULL,
           condition2 = NULL,
           se_type = c("youngs", "constant"),
           ci = TRUE,
           alpha = .05,
           return_condition_pr_mat = FALSE) {

    # -----
    # Check arguments
    # -----
    if (length(all.vars(formula[[3]])) > 1) {
      stop(
        "'formula' must have only one variable on the right-hand side: the ",
        "treatment variable"
      )
    }

    se_type <- match.arg(se_type)

    # -----
    # Parse arguments, clean data
    # -----
    # User can either use declaration or the arguments, not both!
    if (!is.null(declaration)) {
      if (ncol(declaration$probabilities_matrix) > 2) {
        stop(
          "Cannot use horvitz_thompson() with a `declaration` with more than ",
          "two treatment arms for now"
        )
      }

      if (!missing(clusters) |
        !missing(condition_prs) |
        !missing(blocks) |
        !is.null(condition_pr_mat)) {
        stop(
          "Cannot use `declaration` with any of `clusters`, `condition_prs`, ",
          "`blocks`, `condition_pr_mat`"
        )
      }

      # Declaration can only be used if it is the same length as the data being passed
      if (nrow(declaration$probabilities_matrix) != nrow(data)) {
        stop(
          "Cannot use `declaration` if the number of observations in the ",
          "`declaration` 'N' is different than the rows in data"
        )
      }

      # Add clusters, blocks, and treatment probabilities to data so they can be cleaned with clean_model_data
      if (!is.null(declaration$clusters)) {
        if (is.null(data[[".clusters_ddinternal"]])) {
          data[[".clusters_ddinternal"]] <- declaration$clusters
          clusters <- ".clusters_ddinternal"
        } else {
          stop(
            "Can't have a variable called '.clusters_ddinternal' in your `data`",
            "as we use that name for clusters from `declaration` objects."
          )
        }
      }

      if (!is.null(declaration$blocks)) {
        if (is.null(data[[".blocks_ddinternal"]])) {
          data[[".blocks_ddinternal"]] <- declaration$blocks
          blocks <- ".blocks_ddinternal"
        } else {
          stop(
            "Can't have a variable called '.blocks_ddinternal' in your `data`",
            "as we use that name for blocks from `declaration` objects."
          )
        }
      }

      if (is.null(data[[".treatment_prob_ddinternal"]])) {
        if (!is.null(condition2)) {
          treatnum <-
            which(declaration$cleaned_arguments$condition_names == condition2)

          if (length(treatnum) == 0) {
            stop(
              "If `condition2` and `declaration` are both specified, ",
              "`condition2` must match the condition_names in `declaration`.",
              "\n`condition2`: ", condition2, "\n`condition_names`: ",
              paste0(
                declaration$cleaned_arguments$condition_names,
                collapse = ", "
              )
            )
          }

          treatment_prob <- declaration$probabilities_matrix[, treatnum]
        } else {
          # assuming treatment is second column
          treatment_prob <- declaration$probabilities_matrix[, 2]
        }
        data[[".treatment_prob_ddinternal"]] <- treatment_prob
        condition_prs <- ".treatment_prob_ddinternal"
      } else {
        stop(
          "Can't have a variable called '.treatment_prob_ddinternal' in your `data`",
          "as we use that name for treatment probabilites."
        )
      }
    }

    ## Clean data
    where <- parent.frame()
    model_data <- eval(substitute(
      clean_model_data(
        formula = formula,
        data = data,
        subset = subset,
        cluster = clusters,
        condition_pr = condition_prs,
        block = blocks,
        where = where
      )
    ))

    ## condition_pr_mat, if supplied, must be same length
    if (!is.null(condition_pr_mat) && (2 * length(model_data$outcome) != nrow(condition_pr_mat))) {
      stop(
        "After cleaning the data, it has ", length(model_data$outcome), " ",
        "while `condition_pr_mat` has ", nrow(condition_pr_mat), ". ",
        "`condition_pr_mat` should have twice the rows"
      )
    }

    data <- data.frame(
      y = model_data$outcome,
      t = model_data$original_treatment,
      stringsAsFactors = FALSE
    )

    # Parse conditions
    condition_names <- parse_conditions(
      treatment = data$t,
      condition1 = condition1,
      condition2 = condition2,
      estimator = "horvitz_thompson"
    )
    condition2 <- condition_names[[2]]
    condition1 <- condition_names[[1]]

    data$clusters <- model_data$cluster
    data$blocks <- model_data$block
    if (!is.null(model_data$condition_pr)) {
      data$condition_probabilities <- model_data$condition_pr
    }

    # ----------
    # Learn design
    # ----------

    # Declaration is passed
    if (!is.null(declaration)) {

      # Use output from clean_model_data to rebuild declaration
      if (nrow(declaration$probabilities_matrix) != length(data$y)) {
        declaration$probabilities_matrix <- cbind(
          1 - data$condition_probabilities,
          data$condition_probabilities
        )
      }

      # If simple, just use condition probabilities shortcut!
      if (declaration$ra_type == "simple") {
        condition_pr_mat <- NULL
      } else {
        # TODO to allow for declaration with multiple arms, get probability matrix
        # and build it like decl$pr_mat <- cbind(decl$pr_mat[, c(cond1, cond2)])
        condition_pr_mat <- declaration_to_condition_pr_mat(declaration)
      }
    } else if (is.null(condition_pr_mat)) {
      # need to learn it if no declaration and not passed
      # check simple arg
      simple <- ifelse(is.null(simple), TRUE, simple)

      if (is.null(data$blocks) && is.null(data$clusters)) {
        # no blocks or clusters
        if (simple) {
          # don't need condition_pr_mat, just the condition_prs
          # if the user passed it, we're fine and can use just the
          # marginal probabilities
          # if user didn't pass, we have to guess
          if (is.null(data$condition_probabilities)) {
            pr_treat <- mean(data$t == condition2)

            message(
              "Assuming simple random assignment with probability of treatment ",
              "equal to the mean number of obs in `condition2`, which is roughly ",
              round(pr_treat, 3)
            )

            data$condition_probabilities <- pr_treat
          }
        } else {
          # If we don't know the prob, learn it
          if (is.null(data$condition_probabilities)) {
            pr_treat <- mean(data$t == condition2)
            message(
              "Learning probability of complete random assignment from data with",
              "prob = ", round(pr_treat, 3)
            )
            condition_pr_mat <- gen_pr_matrix_complete(
              pr = pr_treat,
              n_total = length(data$y)
            )
          } else {
            if (length(unique(data$condition_probabilities)) > 1) {
              stop(
                "Treatment probabilities must be fixed for complete randomized designs"
              )
            }

            condition_pr_mat <- gen_pr_matrix_complete(
              pr = data$condition_probabilities[1],
              n_total = length(data$y)
            )
          }
        }
      } else if (is.null(data$blocks)) {
        # clustered case
        message(
          "`simple = ", simple, ", using ",
          ifelse(simple, "simple", "complete"), " cluster randomization"
        )

        if (is.null(data$condition_probabilities)) {

          # Split by cluster and get complete randomized values within each cluster
          cluster_treats <- get_cluster_treats(data, condition2)
          pr_treat <- mean(cluster_treats$treat_clust)
          message(
            "`condition_prs` not found, estimating probability of treatment ",
            "to be constant at mean of clusters in `condition2` at prob =", pr_treat
          )

          # Some redundancy in following fn
          condition_pr_mat <- gen_pr_matrix_cluster(
            clusters = data$clusters,
            treat_probs = rep(pr_treat, length(data$clusters)),
            simple = simple
          )
        } else {

          # Just to check if cluster has same treatment within
          get_cluster_treats(data, condition2)

          condition_pr_mat <- gen_pr_matrix_cluster(
            clusters = data$clusters,
            treat_probs = data$condition_probabilities,
            simple = simple
          )
        }
      } else {
        # blocked case
        message(
          "Assuming complete random assignment of clusters within blocks. ",
          "User can use `declaration` or `condition_pr_mat` to have full ",
          "control over the design."
        )

        if (is.null(data$condition_probabilities)) {
          message(
            "`condition_prs` not found, estimating probability of treatment ",
            "to be proportion of units or clusters in condition2 in each block"
          )
          condition_pr_mat <- gen_pr_matrix_block(
            blocks = data$blocks,
            clusters = data$clusters,
            t = data$t,
            condition2 = condition2
          )
        } else {
          condition_pr_mat <- gen_pr_matrix_block(
            blocks = data$blocks,
            clusters = data$clusters,
            p2 = data$condition_probabilities
          )
        }
      }
    }

    # Need the marginal condition prs for later
    if (is.null(data$condition_probabilities)) {
      data$condition_probabilities <-
        diag(condition_pr_mat)[(length(data$y) + 1):(2 * length(data$y))]
    }

    # Check some things that must be true, could do this earlier
    # but don't have condition_probabilities then, and unfortunatey
    # this loops over data clusters a second time
    if (!is.null(data$clusters)) {
      if (any(!tapply(
        data$condition_probabilities,
        data$clusters,
        function(x) all(x == x[1])
      ))) {
        stop("`condition_prs` must be constant within `cluster`")
      }
    }

    rm(model_data)

    # -----
    # Estimation
    # -----

    if (is.null(data$blocks)) {
      return_frame <-
        horvitz_thompson_internal(
          condition_pr_mat = condition_pr_mat,
          condition1 = condition1,
          condition2 = condition2,
          data = data,
          se_type = se_type,
          alpha = alpha
        )
    } else {
      clust_per_block <- check_clusters_blocks(data)

      N <- nrow(data)

      data$index <- 1:N

      block_dfs <- split(data, data$blocks)

      block_estimates <- lapply(block_dfs, function(x) {
        horvitz_thompson_internal(
          data = x,
          condition1 = condition1,
          condition2 = condition2,
          condition_pr_mat = condition_pr_mat[c(x$index, N + x$index), c(x$index, N + x$index)],
          se_type = se_type,
          alpha = alpha
        )
      })

      block_estimates <- do.call(rbind, block_estimates)

      N_overall <- with(block_estimates, sum(N))

      n_blocks <- nrow(block_estimates)

      diff <- with(block_estimates, sum(coefficients * N / N_overall))

      se <- with(block_estimates, sqrt(sum(se ^ 2 * (N / N_overall) ^ 2)))

      return_frame <- data.frame(
        coefficients = diff,
        se = se,
        N = N_overall
      )
    }

    return_frame$df <- NA
    return_list <- add_cis_pvals(return_frame, alpha, ci, ttest = FALSE)

    # -----
    # Build and return output
    # -----

    return_list <- dim_like_return(
      return_list,
      alpha = alpha,
      formula = formula,
      conditions = list(condition1, condition2)
    )

    if (return_condition_pr_mat) {
      return_list[["condition_pr_mat"]] <- condition_pr_mat
    }

    attr(return_list, "class") <- "horvitz_thompson"

    return(return_list)
  }

var_ht_total_no_cov <-
  function(y, ps) {
    sum((1 - ps) * ps * y ^ 2)
  }


horvitz_thompson_internal <-
  function(condition_pr_mat = NULL,
           condition1 = NULL,
           condition2 = NULL,
           data,
           pair_matched = FALSE,
           se_type,
           alpha = .05) {

    # TODO, add estimator from Middleton & Aronow 2015 page 51

    t2 <- which(data$t == condition2)
    t1 <- which(data$t == condition1)

    N <- length(t2) + length(t1)

    se <- NA

    collapsed <- !is.null(data$clusters)
    if (collapsed) {
      if (se_type == "constant") {
        stop(
          "`se_type` = 'constant' only supported for simple random designs ",
          "at the moment"
        )
      }
      # used for cluster randomized designs
      k <- length(unique(data$clusters))

      y2_totals <- tapply(data$y[t2], data$clusters[t2], sum)
      y1_totals <- tapply(data$y[t1], data$clusters[t1], sum)

      to_drop <- which(duplicated(data$clusters))
      t2 <- which(data$t[-to_drop] == condition2)
      t1 <- which(data$t[-to_drop] == condition1)

      # reorder totals because tapply above sorts on cluster
      y2_totals <-
        y2_totals[as.character(data$clusters[-to_drop][t2])]
      y1_totals <-
        y1_totals[as.character(data$clusters[-to_drop][t1])]

      condition_pr_mat <-
        condition_pr_mat[-c(to_drop, N + to_drop), -c(to_drop, N + to_drop)]

      ps2 <- diag(condition_pr_mat)[k + t2]
      ps1 <- diag(condition_pr_mat)[t1]

      # for now rescale, with joint pr need squared top alone
      Y2 <- y2_totals / ps2
      Y1 <- y1_totals / ps1

      diff <- (sum(Y2) - sum(Y1)) / N

      se <-
        sqrt(
          sum(Y2 ^ 2) +
            sum(Y1 ^ 2) +
            ht_var_partial(
              Y2,
              condition_pr_mat[(k + t2), (k + t2), drop = FALSE]
            ) +
            ht_var_partial(
              Y1,
              condition_pr_mat[t1, t1, drop = FALSE]
            ) -
            2 * ht_covar_partial(
              Y2,
              Y1,
              condition_pr_mat[(k + t2), t1, drop = FALSE],
              ps2,
              ps1
            )
        ) / N
    } else {
      # All non-clustered designs

      ps2 <- data$condition_probabilities[t2]
      ps1 <- 1 - data$condition_probabilities[t1]

      Y2 <- data$y[t2] / ps2
      Y1 <- data$y[t1] / ps1

      diff <- (sum(Y2) - sum(Y1)) / N

      if (is.null(condition_pr_mat)) {
        # Simple random assignment
        # joint inclusion probabilities = product of marginals
        if (se_type == "constant") {
          # Scale again
          y0 <- ifelse(
            data$t == condition1,
            data$y / (1 - data$condition_probabilities),
            (data$y - diff) / (1 - data$condition_probabilities)
          )
          y1 <- ifelse(
            data$t == condition2,
            data$y / data$condition_probabilities,
            (data$y + diff) / data$condition_probabilities
          )


          se <-
            sqrt(
              var_ht_total_no_cov(y1, data$condition_probabilities) +
                var_ht_total_no_cov(y0, 1 - data$condition_probabilities) +
                # TODO why is it +2 instead of - (looking at old samii/aronow)
                2 * sum(c(data$y[t2], data$y[t1] + diff) * c(data$y[t2] - diff, data$y[t1]))
            ) / N
        } else {
          # Young's inequality
          se <-
            sqrt(sum(Y2 ^ 2) + sum(Y1 ^ 2)) / N
        }
      } else {
        # Complete random assignment
        if (se_type == "constant") {
          stop(
            "`se_type` = 'constant' only supported for simple random designs ",
            "at the moment"
          )
        } else {
          # Young's inequality
          # this is the "clustered" estimator where each unit is a cluster
          # shouldn't apply to clustered designs but may if user passes a
          # condition_pr_mat
          varN2 <-
            sum(Y2 ^ 2) +
            sum(Y1 ^ 2) +
            ht_var_partial(
              Y2,
              condition_pr_mat[(N + t2), (N + t2), drop = FALSE]
            ) +
            ht_var_partial(
              Y1,
              condition_pr_mat[t1, t1, drop = FALSE]
            ) -
            2 * ht_covar_partial(
              Y2,
              Y1,
              condition_pr_mat[(N + t2), t1, drop = FALSE],
              ps2,
              ps1
            )

          if (!is.nan(varN2)) {
            if (varN2 < 0) {
              warning("Variance below 0")
              se <- NA
            } else {
              se <- sqrt(varN2) / N
            }
          } else {
            warning(
              "Variance is NaN. This is likely the result of a complex condition probability matrix"
            )
            se <- NA
          }
        }
      }
    }

    return_frame <-
      data.frame(
        coefficients = diff,
        se = se,
        N = N,
        stringsAsFactors = FALSE
      )

    return(return_frame)
  }
