// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// These functions return variance for clustered horvitz thompson designs
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double ht_var(const double & p1p2,
              const double & p1,
              const double & p2,
              const double & y1,
              const double & y2) {
  return (p1p2 - p1 * p2) * y1 * y2;
}

// [[Rcpp::export]]
double ht_var_total(const arma::vec & y,
                    const arma::mat & p) {

  double total_variance = 0.0;

  for (unsigned i = 0; i < y.n_elem; ++i) {
    for (unsigned j = 0; j < y.n_elem; ++j) {
      total_variance +=
        ht_var(
          p(i, j),
          p(i, i),
          p(j, j),
          y(i),
          y(j)
        );
    }
  }

  return total_variance;
}

// If p is symmetric the following is ~30pct faster
// [[Rcpp::export]]
double ht_var_total2(const arma::vec & y,
                     const arma::mat & p) {

  double upper_triangle_variance = 0.0;
  double diag_variance = 0.0;

  for (unsigned i = 0; i < y.n_elem; ++i) {
    for (unsigned j = i; j < y.n_elem; ++j) {
      double temp_var =
        ht_var(
          p(i, j),
          p(i, i),
          p(j, j),
          y(i),
          y(j)
        );

      if (i == j) {
        diag_variance += temp_var;
      } else {
        upper_triangle_variance += temp_var;
      }
    }
  }

  return diag_variance + 2.0 * upper_triangle_variance;
}

// [[Rcpp::export]]
double ht_covar_partial(const arma::vec & y1,
                        const arma::vec & y0,
                        const arma::mat & p10,
                        const arma::vec & p1,
                        const arma::vec & p0) {
  double cov_total = 0.0;

  for (unsigned i = 0; i < y1.n_elem; ++i) {
    for(unsigned j = 0; j < y0.n_elem; ++j) {
      cov_total += y1(i) * y0(j) * (p10(i, j) - p1(i) * p0(j)) / p10(i, j);
    }
  }

  return cov_total;
}

// [[Rcpp::export]]
double ht_var_partial(const arma::vec & y,
                      const arma::mat & p) {
  double var_total = 0.0;

  for (unsigned i = 0; i < y.n_elem; ++i) {
    for(unsigned j = 0; j < y.n_elem; ++j) {
      if(i != j) {
        var_total += y(i) * y(j) * (p(i, j) - p(i,i) * p(j,j)) / p(i, j);

      }
    }
  }

  return var_total;
}

// [[Rcpp::export]]
double ht_covar_total(const arma::vec & y0,
                      const arma::vec & y1,
                      const arma::mat & p00,
                      const arma::mat & p11,
                      const arma::mat & pj) {

  double cov_total = 0.0;

  for (unsigned i = 0; i < y0.n_elem; ++i) {
    for (unsigned j = 0; j < y0.n_elem; ++j) {
      if (i != j) {
        cov_total +=
          ht_var(
            pj(i, j),
            p00(i, i),
            p11(j, j),
            y0(i),
            y1(j)
          );
      }
    }
  }

  return cov_total;
}

// [[Rcpp::export]]
double joint_incl_pr(const double & pi,
                     const double & pj,
                     const double & nleft,
                     const double & ntotal) {
  return (2.0 * nleft * pi * pj) / (2.0 * ntotal - pi - pj);
}

// [[Rcpp::export]]
arma::mat gen_pr_matrix_complete(const arma::vec & prs) {

  double n = prs.n_elem;
  arma::mat mat_11(n, n);
  arma::mat mat_00(n, n);
  arma::mat mat_01(n, n);
  arma::mat mat_10(n, n);

  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) {
      if (i == j) {
        mat_11(i, j) = prs[i];
        mat_00(i, j) = 1 - prs[i];
        mat_01(i, j) = 0;
        mat_10(i, j) = 0;
      } else {
        mat_11(i, j) = joint_incl_pr(prs[i], prs[j], n-1, n);
        mat_00(i, j) = joint_incl_pr(1-prs[i], 1-prs[j], n-1, n);
        mat_01(i, j) = joint_incl_pr(1-prs[i], prs[j], n, n);
        mat_10(i, j) = joint_incl_pr(prs[i], 1-prs[j], n, n);
      }
    }
  }

   arma::mat pr_mat = arma::join_rows(
    arma::join_cols(mat_00, mat_01),
    arma::join_cols(mat_10, mat_11)
  );

  return pr_mat;
}

// old
// [[Rcpp::export]]
double ht_var_total_clusters(const arma::vec & y,
                             const arma::vec & ps,
                             const arma::vec & cluster) {
  double total_variance = 0.0;
  arma::vec levels = unique(cluster);

  // iterate over unique cluster values
  for(arma::vec::const_iterator j = levels.begin();
      j != levels.end();
      ++j){

    arma::uvec cluster_ids = find(cluster == *j);
    arma::vec y_cl = y(cluster_ids);
    double ps_cl = ps(cluster_ids[0]);
    for (unsigned i = 0; i < y_cl.n_elem; ++i) {
      for (unsigned j = 0; j < y_cl.n_elem; ++j) {
        total_variance += (1 - ps_cl) * ps_cl * y_cl[i] * y_cl[j];
      }
    }
  }

  return total_variance;

}
