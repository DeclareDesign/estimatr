// [[Rcpp::depends(RcppArmadillo)]]

// These functions implement
// Using Heteroscedasticity Consistent Standard Errors in the Linear Regression Model
// J. Scott Long and Laurie H. Ervin

#include <RcppArmadillo.h>
using namespace Rcpp;

// Faster than x * diagmat(d)
// [[Rcpp::export]]
arma::mat mult_diag(const arma::mat& x, const arma::vec& d) {

  arma::mat out(x.n_rows, x.n_cols);
  for (int j = 0; j < x.n_cols; ++j) {
    out.col(j) = x.col(j) * d(j);
  }

  return out;
}

// (X)^{-1/2}
// [[Rcpp::export]]
arma::mat mat_sq_inv(const arma::mat & X) {
  arma::vec eigvals;
  arma::mat eigvecs;
  arma::eig_sym(eigvals, eigvecs, X);
  arma::vec d2 = 1.0 / arma::sqrt(eigvals);
  return(mult_diag(eigvecs, d2) * arma::trans(eigvecs));
}

// [[Rcpp::export]]
List lm_robust_helper(const arma::vec & y,
                      const arma::mat & X,
                      const Rcpp::Nullable<Rcpp::NumericVector> & cluster,
                      const bool & ci,
                      const String type) {

  // xt
  arma::mat Xt = arma::trans(X);

  //XtX <- t(X_mat) %*% X_mat
  arma::mat XtX = Xt*X;

  // xtx inverse
  arma::mat XtX_inv = arma::inv(XtX);

  // beta_hat <- solve(XtX) %*% t(X_mat) %*% Y
  arma::colvec beta_hat = XtX_inv*Xt*y;

  arma::mat Vcov_hat;
  arma::colvec df(X.n_cols);

  if(type != "none"){

    // residuals
    arma::colvec ei = y - X*beta_hat;

    // hat values
    arma::colvec ei2 = arma::pow(ei, 2);

    double n = X.n_rows;
    double k = X.n_cols;

    arma::sp_mat D_hat(n, n);

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
    } else if ((type == "BM" | type == "stata") & cluster.isNotNull()) {


      // Code adapted from Michal Kolesar
      // https://github.com/kolesarm/Robust-Small-Sample-Standard-Errors

      // Get unique cluster values
      arma::vec clusters = Rcpp::as<arma::vec>(cluster);
      arma::vec levels = unique(clusters);
      double J = levels.n_elem;

      if (type == "BM") {

        arma::mat tutX(J, k);
        arma::cube Gs(n, k, J);

        // iterator used to fill tutX
        int clusternum = 0;

        // iterate over unique cluster values
        for(arma::vec::const_iterator j = levels.begin();
            j != levels.end();
            ++j){

          // Vcov matrix
          arma::uvec cluster_ids = find(clusters == *j);
          int cluster_size = cluster_ids.n_elem;
          arma::mat XtX_inv_tXj = XtX_inv * Xt.cols(cluster_ids);
          // (I_j - Xj %*% (X'X)^{-1} %*% Xj') ^ {-1/2} %*% Xj
          arma::mat minP_Xt = mat_sq_inv(
            arma::eye(cluster_size, cluster_size) - X.rows(cluster_ids) * XtX_inv_tXj
          ) * X.rows(cluster_ids);

          tutX.row(clusternum) = arma::trans(ei(cluster_ids)) * minP_Xt;

          if(ci) {
            // Adjusted degrees of freedom
            arma::mat ss = arma::zeros(n, cluster_size);

            for(int i = 0; i < cluster_size; i++){
              ss(cluster_ids[i], i) = 1;
            }

            Gs.slice(clusternum) = (ss - X * XtX_inv_tXj) * minP_Xt * XtX_inv;
          }

          clusternum++;
        }


        Vcov_hat = XtX_inv * arma::trans(tutX) * tutX * XtX_inv;

        if (ci) {
          for(int p = 0; p < k; p++){
            arma::mat G = Gs.subcube(arma::span::all, arma::span(p), arma::span::all);
            arma::mat GG = arma::trans(G) * G;
            df[p] = std::pow(arma::trace(GG), 2) / arma::accu(arma::pow(GG, 2));
          }
        }
      }

      if (type == "stata") {

        arma::mat XteetX(k, k);
        XteetX.fill(0.0);

        // iterate over unique cluster values
        for(arma::vec::const_iterator j = levels.begin();
            j != levels.end();
            ++j){
          arma::uvec cluster_ids = find(clusters == *j);

          XteetX += Xt.cols(cluster_ids) * ei.elem(cluster_ids) *
            arma::trans(ei.elem(cluster_ids)) * X.rows(cluster_ids);
        }

        Vcov_hat = ((J * (n - 1)) / ((J - 1) * (n - k))) * XtX_inv * XteetX * XtX_inv;

        df.fill(J - 1);

      }

    }
  }

  return List::create(_["beta_hat"]= beta_hat,
                      _["Vcov_hat"]= Vcov_hat,
                      _["df"]= df);
}



