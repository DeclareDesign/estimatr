// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppEigen.h>
using namespace Rcpp;

// These functions help compute the variance for the Horvitz-Thompson estimators
// TODO use symmetry and matrices to improve speed

// [[Rcpp::export]]
double ht_covar_partial(const Eigen::VectorXd & y1,
                        const Eigen::VectorXd & y0,
                        const Eigen::MatrixXd & p10,
                        const Eigen::VectorXd & p1,
                        const Eigen::VectorXd & p0) {
  double cov_total = 0.0;

  for (int i = 0; i < y1.size(); ++i) {
    for(int j = 0; j < y0.size(); ++j) {
      if(p10(i, j) == 0) {
        cov_total += y1(i) * y0(j) * (p10(i, j) - p1(i) * p0(j));
      } else {
        cov_total += y1(i) * y0(j) * (p10(i, j) - p1(i) * p0(j)) / p10(i, j);
      }
    }
  }

  return cov_total;
}

// [[Rcpp::export]]
double ht_var_partial(const Eigen::VectorXd & y,
                      const Eigen::MatrixXd & p) {
  double var_total = 0.0;

  for (int i = 0; i < y.size(); ++i) {
    for(int j = 0; j < y.size(); ++j) {
      if(i != j) {
        if (p(i, j) == 0) {
          var_total += y(i) * y(j) * (p(i, j) - p(i,i) * p(j,j)) +
            std::pow(y(i), 2) * p(i, i) / 2.0 + std::pow(y(j), 2) * p(j, j) / 2.0;
        } else {
          var_total += y(i) * y(j) * (p(i, j) - p(i,i) * p(j,j)) / p(i, j);
        }

      }
    }
  }

  return var_total;
}
