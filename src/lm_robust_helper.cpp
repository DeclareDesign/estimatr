// [[Rcpp::depends(RcppArmadillo)]]

// These functions implement
// Using Heteroscedasticity Consistent Standard Errors in the Linear Regression Model
// J. Scott Long and Laurie H. Ervin

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List lm_robust_helper(const arma::vec & y, const arma::mat & X, const String type) {

  // xt
  arma::mat Xt = arma::trans(X);

  //XtX <- t(X_mat) %*% X_mat
  arma::mat XtX = Xt*X;

  // xtx inverse
  arma::mat XtX_inv = arma::inv(XtX);

  // beta_hat <- solve(XtX) %*% t(X_mat) %*% Y
  arma::colvec beta_hat = XtX_inv*Xt*y;

  // residuals
  arma::colvec ei = y - X*beta_hat;

  // hat values
  arma::colvec ei2 = arma::pow(ei, 2);

  double n = X.n_rows;
  double k = X.n_cols;

  arma::sp_mat D_hat(n, n);
  arma::mat Vcov_hat;

  if(type == "HC0"){
    D_hat.diag() = ei2;
    Vcov_hat = XtX_inv * Xt * D_hat * X * XtX_inv;

  } else if (type == "HC1") {
    D_hat.diag() = ei2;
    Vcov_hat = n/(n-k) * XtX_inv * Xt * D_hat * X * XtX_inv;
  } else if (type == "HC2") {
    int n = X.n_rows;
    arma::colvec hii(n);

    for(int i = 0; i < n; i++){
      arma::rowvec Xi = X.row(i);
      hii[i] = arma::as_scalar(Xi * XtX_inv * arma::trans(Xi));
    }
    D_hat.diag() = ei2 /(1 - hii);
    Vcov_hat = XtX_inv * Xt * D_hat * X * XtX_inv;
  } else if (type == "HC3") {
    int n = X.n_rows;
    arma::colvec hii(n);

    for(int i = 0; i < n; i++){
      arma::rowvec Xi = X.row(i);
      hii[i] = arma::as_scalar(Xi * XtX_inv * arma::trans(Xi));
    }
    D_hat.diag() = ei2 /pow(1 - hii, 2);
    Vcov_hat = XtX_inv * Xt * D_hat * X * XtX_inv;
  } else if (type == "classical") {
    // the next line due to Dirk Eddelbuettel
    // http://dirk.eddelbuettel.com/code/rcpp.armadillo.html
    double s2 = std::inner_product(ei.begin(), ei.end(), ei.begin(), 0.0)/(n - k);
    Vcov_hat = s2 * XtX_inv;
  }

  return List::create(_["beta_hat"]= beta_hat, _["Vcov_hat"]= Vcov_hat);
}
