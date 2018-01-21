// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppEigen.h>
using namespace Rcpp;

// These functions help compute the variance for the Horvitz-Thompson estimators

// [[Rcpp::export]]
double ht_var(const double & p1p2,
              const double & p1,
              const double & p2,
              const double & y1,
              const double & y2) {
  return (p1p2 - p1 * p2) * y1 * y2;
}

// [[Rcpp::export]]
double ht_var_total(const Eigen::VectorXd & y,
                    const Eigen::MatrixXd & p) {

  double total_variance = 0.0;

  for (int i = 0; i < y.size(); ++i) {
    for (int j = 0; j < y.size(); ++j) {
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

// // [[Rcpp::export]]
// double ht_covar_total(const Eigen::VectorXd & y0,
//                       const Eigen::VectorXd & y1,
//                       const Eigen::MatrixXd & p00,
//                       const Eigen::MatrixXd & p11,
//                       const Eigen::MatrixXd & pj) {
//
//   double cov_total = 0.0;
//
//   for (int i = 0; i < y0.size(); ++i) {
//     for (int j = 0; j < y0.size(); ++j) {
//       if (i != j) {
//         cov_total +=
//           ht_var(
//             pj(i, j),
//             p00(i, i),
//             p11(j, j),
//             y0(i),
//             y1(j)
//           );
//       }
//     }
//   }
//
//   return cov_total;
// }

// unused right now
// If p is symmetric the following is ~30pct faster
// // [[Rcpp::export]]
// double ht_var_total2(const Eigen::VectorXd & y,
//                      const Eigen::MatrixXd & p) {
//
//   double upper_triangle_variance = 0.0;
//   double diag_variance = 0.0;
//
//   for (int i = 0; i < y.size(); ++i) {
//     for (int j = i; j < y.size(); ++j) {
//       double temp_var =
//         ht_var(
//           p(i, j),
//           p(i, i),
//           p(j, j),
//           y(i),
//           y(j)
//         );
//
//       if (i == j) {
//         diag_variance += temp_var;
//       } else {
//         upper_triangle_variance += temp_var;
//       }
//     }
//   }
//
//   return diag_variance + 2.0 * upper_triangle_variance;
// }
