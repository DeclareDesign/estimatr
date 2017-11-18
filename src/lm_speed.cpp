// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

// These functions implement
// Using Heteroscedasticity Consistent Standard Errors in the Linear Regression Model
// J. Scott Long and Laurie H. Ervin
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;

// Much of what follows is modified from RcppEigen Vignette by Douglas Bates and Dirk Eddelbuettel
// https://cran.r-project.org/web/packages/RcppEigen/vignettes/RcppEigen-Introduction.pdf
Eigen::MatrixXd AtA(const Eigen::Map<Eigen::MatrixXd>& A) {
  int n(A.cols());
  return Eigen::MatrixXd(n,n).setZero().selfadjointView<Eigen::Lower>()
                             .rankUpdate(A.adjoint());
}

//' @export
// [[Rcpp::export]]
List lm_ei_test(Eigen::Map<Eigen::MatrixXd>& X,
                const Eigen::Map<Eigen::VectorXd>& y,
                const bool& chol,
                const bool& trychol) {

  const int n(X.rows()), p(X.cols());
  const Eigen::MatrixXd Xt(X.transpose());
  Eigen::MatrixXd XtX_inv, R_inv, Vcov_hat;
  Eigen::VectorXd beta_hat(Eigen::VectorXd::Constant(p, ::NA_REAL));
  Eigen::VectorXd dof(p);

  if (chol) {
    const Eigen::LLT<Eigen::MatrixXd> llt(Xt*X);
    beta_hat = llt.solve(X.adjoint() * y);
    R_inv = llt.matrixL().solve(Eigen::MatrixXd::Identity(p, p));
    XtX_inv = R_inv.transpose() * R_inv;
  } else {
    try {
      if (trychol) {
        const Eigen::LLT<Eigen::MatrixXd> llt(Xt*X);

        // Catch case where X is rank-deficient
        if(llt.info() == Eigen::NumericalIssue) {
          throw std::runtime_error("Possibly non semi-positive definite matrix!");
        } else {
          beta_hat = llt.solve(X.adjoint() * y);
          R_inv = llt.matrixL().solve(Eigen::MatrixXd::Identity(p, p));
          XtX_inv = R_inv.transpose() * R_inv;
        }
      } else {
        throw std::runtime_error("move to QR");
      }

    } catch (const std::exception& e) {

      Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR(X);
      const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType Pmat(PQR.colsPermutation());

      const int r(PQR.rank());

      R_inv = PQR.matrixQR().topLeftCorner(r, r).triangularView<Eigen::Upper>().solve(Eigen::MatrixXd::Identity(r, r));
      Eigen::VectorXd effects(PQR.householderQ().adjoint() * y);
      beta_hat.head(r) = R_inv * effects.head(r);
      // Rcout << beta_hat << std::endl << std::endl;
      beta_hat = Pmat * beta_hat;

      // Get all column indices
      Eigen::ArrayXi Pmat_indices = Pmat.indices();
      // Get the order for the columns you are keeping
      Eigen::ArrayXi Pmat_keep = Pmat_indices.head(r);
      // Get the indices for columns you are discarding
      Eigen::ArrayXi Pmat_toss = Pmat_indices.tail(p - r);

      // this code takes the indices you are keeping, and, while preserving order, keeps them in the range [0, r-1]
      // For each one, see how many dropped indices are smaller, and subtract that difference
      // Ex: p = 4, r = 2
      // Pmat_indices = {3, 1, 0, 2}
      // Pmat_keep = {3, 1}
      // Pmat_toss = {0, 2}
      // Now we go through each keeper, count how many in toss are smaller, and then modify accordingly
      // 3 - 2 and 1 - 1
      // Pmat_keep = {1, 0}
      for(Eigen::Index i=0; i<r; ++i)
      {
        Pmat_keep(i) = Pmat_keep(i) - (Pmat_toss < Pmat_keep(i)).count();
      }
      // Rcout << Pmat_indices << std::endl << std::endl;
      // Rcout << Pmat_keep << std::endl << std::endl;
      //
      // Rcout << Pmat * Eigen::MatrixXd::Identity(p, p) << std::endl << std::endl;

      Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P = Eigen::PermutationWrapper<Eigen::ArrayXi>(Pmat_keep);
      // Rcout << P * Eigen::MatrixXd::Identity(r, r) << std::endl << std::endl;

      R_inv = P * R_inv * P;
      XtX_inv = R_inv * R_inv.transpose();
    }
  }

  return List::create(_["beta_hat"]= beta_hat,
                      _["Vcov_hat"]= Vcov_hat,
                      _["XtX_inv"]= XtX_inv,
                      _["dof"]= dof);
}

//' @export
// [[Rcpp::export]]
List lm_old(arma::mat & X,
                      const arma::vec & y) {

  // t(X)
  arma::mat Xt = arma::trans(X);

  // xtx inverse
  arma::mat XtX_inv = arma::inv(Xt * X);

  // beta_hat <- solve(XtX) %*% t(X_mat) %*% Y
  arma::colvec beta_hat = XtX_inv*Xt*y;

  arma::mat Vcov_hat;

  arma::colvec dof(X.n_cols);
  dof.fill(-99);


  return List::create(_["beta_hat"]= beta_hat,
                      _["XtX_inv"]= XtX_inv,
                      _["Vcov_hat"]= Vcov_hat,
                      _["dof"]= dof);
}
