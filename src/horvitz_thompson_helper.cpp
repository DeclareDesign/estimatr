// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// These functions return variance for clustered horvitz thompson designs
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace Rcpp;

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
