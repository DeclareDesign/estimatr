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
double ht_covar_total(const arma::vec & y0,
                      const arma::vec & y1,
                      const arma::mat & p00,
                      const arma::mat & p11,
                      const arma::mat & pj) {

  double cov_total = 0.0;

  for (unsigned i = 0; i < y0.n_elem; ++i) {
    for (unsigned j = 0; j < y0.n_elem; ++j) {
      Rcpp::Rcout << i << std::endl << j << std::endl;
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
