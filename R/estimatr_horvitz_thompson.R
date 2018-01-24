#' Horvitz-Thompson estimator of treatment effects
#'
#' @param formula An object of class "formula", such as Y ~ Z
#' @param data A data.frame.
#' @param condition_prs An optional bare (unquoted) name of the variable with the condition 2 (treatment) probabilities.
#' @param blocks An optional bare (unquoted) name of the block variable. Use for blocked designs only.
#' @param clusters An optional bare (unquoted) name of the variable that corresponds to the clusters in the data; used for cluster randomized designs. For blocked designs, clusters must be within blocks.
#' @param simple An optional boolean for whether the randomization is simple (TRUE) or complete (FALSE). This is ignored if \code{blocks} are specified, as all blocked designs use complete randomization, or either \code{declaration} or \code{condition_pr_mat} are passed. Otherwise, it defaults to \code{TRUE}.
#' @param condition_pr_mat An optional 2n * 2n matrix of marginal and joint probabilities of all units in condition1 and condition2. Takes precedent over any other means of specifying the design as it can encapsulate everything relevant about the design. See details.
#' @param declaration An object of class \code{"ra_declaration"}, from the \code{\link{randomizr}} package that is an alternative way of specifying the design. Cannot be used along with any of \code{condition_prs}, \code{blocks}, \code{clusters}, or \code{condition_pr_mat}. See details.
#' @param subset An optional bare (unquoted) expression specifying a subset of observations to be used.
#' @param se_type can be one of \code{c("youngs", "constant")} and correspond's to estimating the standard errors using Young's inequality (default, conservative), or the constant effects assumption.
#' @param collapsed A boolean used to collapse clusters to their cluster totals for variance estimation, FALSE by default.
#' @param ci A boolean for whether to compute and return pvalues and confidence intervals, TRUE by default.
#' @param alpha The significance level, 0.05 by default.
#' @param condition1 values of the conditions to be compared. Effects are estimated with condition1 as control and condition2 as treatment. If unspecified, condition1 is the "first" condition and condition2 is the "second" according to r defaults.
#' @param condition2 values of the conditions to be compared. Effects are estimated with condition1 as control and condition2 as treatment. If unspecified, condition1 is the "first" condition and condition2 is the "second" according to r defaults.
#' @param return_condition_pr_mat a boolean for whether to return the condition probability matrix
#'
#' @details This function implements the Horvitz-Thompson estimator for treatment effects. This estimator is useful for estimating unbiased treatment effects given any randomization scheme as long as the randomization scheme is well known. See the \href{http://estimatr.declaredesign.org/articles/technical-notes.html}{technical notes} for more information and references.
#'
#' Horvitz-Thompson relies on knowledge of the probability each unit was in condition2, in condition1, and the joint probability with every other unit for all combinations of treatment conditions. This design information can be passed to the \code{horvitz_thompson} function in three ways.
#'
#' @export
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
#' # Simple random assignment
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
#' # Complete random assignemtn
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
#' # Clustered treatment, complete random assigment
#' # Simulating data
#' dat$cl <- rep(1:4, times = c(2, 2, 3, 3))
#' clust_crs_decl <- declare_ra(N = nrow(dat), clusters = dat$cl, prob = 0.5)
#' dat$z <- conduct_ra(clust_crs_decl)
#' # Regular SE using Young's inequality
#' horvitz_thompson(y ~ z, data = dat, declaration = clust_crs_decl)
#' # SE using collapsed cluster totals and Young's inequality
#' horvitz_thompson(y ~ z, data = dat, declaration = clust_crs_decl, collapsed = TRUE)
#'
horvitz_thompson <-
  function(formula,
           data,
           condition_prs,
           blocks,
           clusters,
           simple = NULL,
           condition_pr_mat = NULL,
           declaration = NULL,
           subset,
           se_type = c('youngs', 'constant'),
           collapsed = FALSE,
           ci = TRUE,
           alpha = .05,
           condition1 = NULL,
           condition2 = NULL,
           return_condition_pr_mat = FALSE) {

    #-----
    # Check arguments
    #-----
    if (length(all.vars(formula[[3]])) > 1) {
      stop(
        "'formula' must have only one variable on the right-hand side: the ",
        "treatment variable."
      )
    }

    se_type <- match.arg(se_type)

    #-----
    # Parse arguments, clean data
    #-----
    # User can either use declaration or the arguments, not both!
    if (!is.null(declaration)) {

      if (ncol(declaration$probabilities_matrix) > 2) {
        stop(
          "Cannot use horvitz_thompson() with a `declaration` with more than ",
          "two treatment arms for now."
        )
      }

      if (!missing(clusters) |
          !missing(condition_prs) |
          !missing(blocks) |
          !is.null(condition_pr_mat)) {
        stop("Cannot use declaration with any of clusters, condition_prs, blocks, condition_pr_mat.")
      }

      # Declaration can only be used if it is the same length as the data being passed
      if (nrow(declaration$probabilities_matrix) != nrow(data)) {
        stop("Cannot use declaration if declaration 'N' is different than the rows in data.")
      }

      # Add clusters, blocks, and treatment probabilities to data so they can be cleaned with clean_model_data
      if (!is.null(declaration$clusters)) {
        if (is.null(data[[".clusters_ddinternal"]])) {
          data[[".clusters_ddinternal"]] <- declaration$clusters
          clusters <- '.clusters_ddinternal'
        } else {
          stop("estimatr stores clusters from declarations in a variable called .clusters_ddinternal in your data. Please remove it and try again.")
        }
      }

      if (!is.null(declaration$blocks)) {
        if (is.null(data[[".blocks_ddinternal"]])) {
          data[[".blocks_ddinternal"]] <- declaration$blocks
          blocks <- '.blocks_ddinternal'
        } else {
          stop("estimatr stores blocks from declarations in a variable called .blocks_ddinternal in your data. Please remove it and try again.")
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
        condition_prs <- '.treatment_prob_ddinternal'
      } else {
        stop(
          "Cannot have a variable in your `data` called ",
          "'.treatment_prob_ddinternal' as this function uses that to store ",
          "the treatment probabilities.\nPlease remove/rename it and try again"
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
    if (!is.null(condition_pr_mat) && (2*length(model_data$outcome) != nrow(condition_pr_mat))) {
      stop(
        "After cleaning the data, it has ", length(model_data$outcome), " ",
        "while condition_pr_mat has ", nrow(condition_pr_mat), ". ",
        "condition_pr_mat should have twice the rows"
      )
    }

    data <- data.frame(y = model_data$outcome,
                       t = model_data$original_treatment)

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

    #----------
    # Learn design
    #----------

    # Declaration is passed
    if (!is.null(declaration)) {

      # Use output from clean_model_data to rebuild declaration
      if (nrow(declaration$probabilities_matrix) != length(data$y)) {
        declaration$probabilities_matrix <- cbind(1 - data$condition_probabilities,
                                                  data$condition_probabilities)
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
              "equal to the mean number of obs in condition2, which is roughly ",
              round(pr_treat, 3)
            )

            data$condition_probabilities <- pr_treat
          }
        } else {
          # If we don't know the prob, learn it
          if (is.null(data$condition_probabilities)) {
            pr_treat <- mean(data$t == condition2)
            message(
              "Learning probability of complete random assignment from data ",
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
            "to be constant at mean of clusters in condition2 at: ", pr_treat
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
            t = data$t
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
        diag(condition_pr_mat)[(length(data$y)+1):(2*length(data$y))]
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

    #-----
    # Estimation
    #-----

    if (is.null(data$blocks)) {
      return_frame <-
        horvitz_thompson_internal(
          condition_pr_mat = condition_pr_mat,
          condition1 = condition1,
          condition2 = condition2,
          data = data,
          collapsed = collapsed,
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
          collapsed = collapsed,
          se_type = se_type,
          alpha = alpha
        )
      })

      block_estimates <- do.call(rbind, block_estimates)

      N_overall <- with(block_estimates, sum(N))

      n_blocks <- nrow(block_estimates)

      diff <- with(block_estimates, sum(est * N/N_overall))

      se <- with(block_estimates, sqrt(sum(se^2 * (N/N_overall)^2)))

      return_frame <- data.frame(
        est = diff,
        se = se,
        N = N_overall
      )

    }

    return_frame$df <- NA
    return_list <- add_cis_pvals(return_frame, alpha, ci, ttest = FALSE)

    #-----
    # Build and return output
    #-----

    return_list <- dim_like_return(return_list,
                                   alpha = alpha,
                                   formula = formula,
                                   conditions = list(condition1, condition2))

    if (return_condition_pr_mat) {
      return_list[["condition_pr_mat"]] <- condition_pr_mat
    }

    attr(return_list, "class") <- "horvitz_thompson"

    return(return_list)
  }

var_ht_total_no_cov <-
  function(y, ps) {
    sum((1 - ps) * ps * y^2)
  }


horvitz_thompson_internal <-
  function(condition_pr_mat = NULL,
           condition1 = NULL,
           condition2 = NULL,
           data,
           pair_matched = FALSE,
           collapsed,
           se_type,
           alpha = .05) {

    #if (length(condition_names) == 1) {
    #  stop("Must have units with both treatment conditions within each block.")
    #}

    # print(table(condition_pr_mat))

    t2 <- which(data$t == condition2)
    t1 <- which(data$t == condition1)

    ps2 <- data$condition_probabilities[t2]
    ps1 <- 1 - data$condition_probabilities[t1]

    Y2 <- data$y[t2] / ps2
    Y1 <- data$y[t1] / ps1

    N <- length(Y2) + length(Y1)

    # Estimator from Middleton & Aronow 2015 page 51
    # TODO figure out why below is not equivalent for constant pr cluster randomized exps
    # if (!is.null(data$clusters) & estimator == 'ma') {
    #
    #   print('ma')
    #   k <- length(unique(data$clusters))
    #   c_2 <- data$clusters[t2]
    #   c_1 <- data$clusters[t1]
    #   k_2 <- length(unique(c_2))
    #   k_1 <- length(unique(c_1))
    #
    #   totals_2 <- tapply(data$y[t2], c_2, sum)
    #   totals_1 <- tapply(data$y[t1], c_1, sum)
    #
    #   diff <- (mean(totals_2) - mean(totals_1)) * k / N
    #
    #   se <-
    #     sqrt(
    #       k^2 / N^2 * (
    #         mean((totals_2 - mean(totals_2))^2) / (k_2 - 1) +
    #         mean((totals_1 - mean(totals_1))^2) / (k_1 - 1)
    #       )
    #     )
    #
    # } else {
    diff <- (sum(Y2) - sum(Y1)) / N

    se <- NA

    if (collapsed) {

      # TODO this may not work with 3 valued treatments as it
      # may not subset everything correctly
      if (is.null(data$clusters)) {
        stop(
          "The collapsed estimator only works if you either pass a ",
          "declaration with clusters or explicitly specify the clusters."
        )
      }

      k <- length(unique(data$clusters))

      y2_totals <- tapply(data$y[t2], data$clusters[t2], sum)
      y1_totals <- tapply(data$y[t1], data$clusters[t1], sum)

      to_drop <- which(duplicated(data$clusters))
      t2 <- which(data$t[-to_drop] == condition2)
      t1 <- which(data$t[-to_drop] == condition1)

      # print(t2)
      # print(t1)
      #
      # print(to_drop)
      # print(y2_totals)
      # print(y1_totals)
      # reorder totals because tapply will sort on cluster
      y2_totals <- y2_totals[as.character(data$clusters[-to_drop][t2])]
      y1_totals <- y1_totals[as.character(data$clusters[-to_drop][t1])]
      # print(y2_totals)
      # print(y1_totals)

      condition_pr_mat <- condition_pr_mat[-c(to_drop, N + to_drop), -c(to_drop, N + to_drop)]
      # print(dimnames(condition_pr_mat))
      # print(k)
      ps2 <- diag(condition_pr_mat)[k + t2]
      ps1 <- diag(condition_pr_mat)[t1]
      # for now rescale, with joint pr need squared top alone
      Y2 <- y2_totals / ps2
      Y1 <- y1_totals / ps1
#
#       print(sum(Y2^2))
#       print(sum(Y1^2))
#
#       print((y2_totals/ps2)^2)
#       print((y1_totals/ps1)^2)
      diff <- (sum(Y2) - sum(Y1)) / N

      se <-
        sqrt(
          sum(Y2^2) +
          sum(Y1^2) +
          ht_var_partial(
            Y2,
            condition_pr_mat[(k + t2), (k + t2), drop = F]
          ) +
          ht_var_partial(
            Y1,
            condition_pr_mat[t1, t1, drop = F]
          ) -
          2 * ht_covar_partial(Y2,
                               Y1,
                               condition_pr_mat[(k + t2), t1, drop = F],
                               ps2,
                               ps1)
        ) / N

    } else if (is.null(condition_pr_mat)) {
      # Simple random assignment
      # joint inclusion probabilities = product of marginals
      if (se_type ==  'constant') {
        # rescaled potential outcomes
        y0 <- ifelse(data$t==condition1,
                     data$y / (1-data$condition_probabilities),
                     (data$y - diff) / (1-data$condition_probabilities))
        y1 <- ifelse(data$t==condition2,
                     data$y / data$condition_probabilities,
                     (data$y + diff) / data$condition_probabilities)

        # print(var_ht_total_no_cov(y1, data$condition_probabilities) +
        #       var_ht_total_no_cov(y0, 1 - data$condition_probabilities))
        #
        # print(-2*sum(c(data$y[t2], data$y[t1] + diff) * c(data$y[t2] - diff, data$y[t1])))

        se <-
          sqrt(
            (
              var_ht_total_no_cov(y1, data$condition_probabilities) +
              var_ht_total_no_cov(y0, 1 - data$condition_probabilities) +
              # TODO why is it +2 instead of - (looking at old samii/aronow)
              2 * sum(c(data$y[t2], data$y[t1] + diff) * c(data$y[t2] - diff, data$y[t1]))

            )
          ) / N
      } else {
        se <-
          sqrt(
            sum(Y2^2) + sum(Y1^2)
          ) / N
      }
    } else {
      # Complete random assignment
      if (se_type ==  'constant') {
        stop(
          "`se_type` = 'constant' only supported for simple random designs ",
          "at the moment"
        )
        # rescaled potential outcomes
        # y0 <- ifelse(data$t==condition1,
        #              data$y / (1-data$condition_probabilities),
        #              (data$y - diff) / (1-data$condition_probabilities))
        # y1 <- ifelse(data$t==condition2,
        #              data$y / data$condition_probabilities,
        #              (data$y + diff) / data$condition_probabilities)
        # se <-
        #   sqrt(
        #     ht_var_total(
        #       c(-y0, y1),
        #       condition_pr_mat
        #     )
        #   ) / N
      } else {
        #print('full youngs')
        #print(condition_pr_mat)
        # Young's inequality

        varN2 <-
          sum(Y2^2) +
            sum(Y1^2) +
            ht_var_partial(
              Y2,
              condition_pr_mat[(N + t2), (N + t2), drop = F]
            ) +
            ht_var_partial(
              Y1,
              condition_pr_mat[t1, t1, drop = F]
            ) -
            2 * ht_covar_partial(Y2,
                                 Y1,
                                 condition_pr_mat[(N + t2), t1, drop = F],
                                 ps2,
                                 ps1)

        if (!is.nan(varN2)) {
          if (varN2 < 0) {
            warning("Variance below 0, consider using constant effects assumption or a different estimator")
            se <- NA
          } else {
            se <- sqrt(varN2) / N
          }

        } else {
          warning("Variance is NaN. This is likely the result of a complex condition probability matrix")
          se <- NA
        }
      }
    }
    # }

    return_frame <-
      data.frame(
        est = diff,
        se = se,
        N = N
      )

    return(return_frame)

  }
