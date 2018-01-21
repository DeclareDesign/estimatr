// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppEigen.h>
using namespace Rcpp;

// These functions help build condition probability matrices for Horvitz-Thompson

// [[Rcpp::export]]
double joint_incl_pr(const double & pi,
                     const double & pj,
                     const int & same,
                     const double & ntotal) {
  return pi * ((pj * ntotal - same) / (ntotal - 1));
}

// [[Rcpp::export]]
Eigen::MatrixXd gen_pr_matrix_complete(const Eigen::VectorXd & prs) {

  int n = prs.size();
  Eigen::MatrixXd pr_mat(2*n, 2*n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) {
        pr_mat(i, j) = 1 - prs(i);
        pr_mat(i + n, j) = 0;
        pr_mat(i, j + n) = 0;
        pr_mat(i + n, j + n) = prs(i);
      } else {
        pr_mat(i, j) = joint_incl_pr(1-prs(i), 1-prs(j), 1, n);
        pr_mat(i + n, j) = joint_incl_pr(prs(i), 1-prs(j), 0, n);
        pr_mat(i, j + n) = joint_incl_pr(1-prs(i), prs(j), 0, n);
        pr_mat(i + n, j + n) = joint_incl_pr(prs(i), prs(j), 1, n);
      }
    }
  }

  return pr_mat;
}
