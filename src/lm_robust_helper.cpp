// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// These functions implement
// Using Heteroscedasticity Consistent Standard Errors in the Linear Regression Model
// J. Scott Long and Laurie H. Ervin
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace Rcpp;


//----------
// Helper functions
//----------

// This does a matrix multiplication of a matrix X with a diagonal matrix where
// the diagonal is D.
// This is faster than x * arma::diagmat(d)
// [[Rcpp::export]]
arma::mat mult_diag(const arma::mat& x, const arma::vec& d) {

  arma::mat out(x.n_rows, x.n_cols);
  for (unsigned j = 0; j < x.n_cols; ++j) {
    out.col(j) = x.col(j) * d(j);
  }

  return out;
}

// Get's the inverse square root of a matrix (X)^{-1/2}
// Used in computing the BM standard errors (I - P_ss^&{-1/2}
// [[Rcpp::export]]
arma::mat mat_sqrt_inv(const arma::mat & X, const bool & tol) {
  arma::vec eigvals;
  arma::mat eigvecs;
  arma::eig_sym(eigvals, eigvecs, X);
  for (unsigned i = 0; i < eigvals.n_elem; ++i) {
    if (tol) {
      if (eigvals(i) > std::pow(10.0, -12.0)) {
        eigvals(i) = 1.0 / std::sqrt(eigvals(i));
      } else {
        eigvals(i) = 0;
      }
    } else {
      eigvals(i) = 1.0 / std::sqrt(eigvals(i));
    }
  }
  return(mult_diag(eigvecs, eigvals) * arma::trans(eigvecs));
}



//----------
// Main LM fit and SE function
//----------
// [[Rcpp::export]]
List lm_robust_helper(const arma::vec & y,
                      arma::mat & X,
                      const Rcpp::Nullable<Rcpp::NumericMatrix> & Xunweighted,
                      const Rcpp::Nullable<Rcpp::NumericVector> & weight,
                      const double & weight_mean,
                      const Rcpp::Nullable<Rcpp::NumericVector> & cluster,
                      const bool & ci,
                      const String type,
                      const std::vector<bool> & which_covs) {

  arma::mat Xoriginal(X.n_rows, X.n_cols);
  arma::vec weights;
  if(Xunweighted.isNotNull()) {
    weights = Rcpp::as<arma::vec>(weight);
    Xoriginal = Rcpp::as<arma::mat>(Xunweighted);
  }

  // t(X)
  arma::mat Xt = arma::trans(X);

  // xtx inverse
  arma::mat XtX_inv = arma::inv(Xt * X);

  // beta_hat <- solve(XtX) %*% t(X_mat) %*% Y
  arma::colvec beta_hat = XtX_inv*Xt*y;

  arma::mat Vcov_hat;
  arma::colvec dof(X.n_cols);
  dof.fill(-99);

  if(type != "none"){

    // residuals
    arma::colvec ei = y - X * beta_hat;

    double n = X.n_rows;
    double k = X.n_cols;

    if (type == "classical") {

      // the next line due to Dirk Eddelbuettel
      // http://dirk.eddelbuettel.com/code/rcpp.armadillo.html
      double s2 = std::inner_product(ei.begin(), ei.end(), ei.begin(), 0.0)/(n - k);
      Vcov_hat = s2 * XtX_inv;

    } else if ( (type == "HC0") | (type == "HC1") | (type == "HC2") | (type == "HC3") ) {

      // hat values
      arma::colvec ei2 = arma::pow(ei, 2);
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
      }

    } else if ( (type == "stata") | (type == "CR2") ) {

      // Code adapted from Michal Kolesar
      // https://github.com/kolesarm/Robust-Small-Sample-Standard-Errors

      // Get unique cluster values
      arma::vec clusters = Rcpp::as<arma::vec>(cluster);
      arma::vec levels = unique(clusters);
      double J = levels.n_elem;

      if (type == "CR2") {

        // Instead of having X and Xt be weighted by sqrt(W), this has them weighted by W
        Xt = Xt.each_row() % arma::sqrt(arma::trans(weights));
        X = X.each_col() % arma::sqrt(weights);

        arma::mat tutX(J, k);


        // used for the dof corrction
        arma::cube H1s(k, k, J);
        arma::cube H2s(k, k, J);
        arma::cube H3s(k, k, J);
        arma::mat P_diags(k, J);

        arma::mat M_U_ct = arma::trans(arma::chol(XtX_inv));

        // Unweighted residuals
        ei = y / sqrt(weights) - Xoriginal * beta_hat;

        int clusternum = 0;

        arma::mat MUWTWUM = XtX_inv * Xt * arma::trans(Xt) * XtX_inv;
        arma::mat Omega_ct = arma::trans(arma::chol(MUWTWUM));

        // iterate over unique cluster values
        for(arma::vec::const_iterator j = levels.begin();
            j != levels.end();
            ++j){

          arma::uvec cluster_ids = find(clusters == *j);
          int cluster_size = cluster_ids.n_elem;
          arma::mat P(cluster_size, k);

          arma::mat D = arma::eye(cluster_size, cluster_size);
          arma::mat H = Xoriginal.rows(cluster_ids) * XtX_inv * Xt.cols(cluster_ids);
          arma::mat I_min_H = D - H;

          // Code from clubSandwich
          // uwTwu <- Map(function(uw, th) uw %*% th %*% t(uw),
          //             uw = UW_list, th = Theta_list)
          // MUWTWUM <- M_U %*% Reduce("+", uwTwu) %*% M_U

          //(thet - h %*% thet - thet %*% t(h) + u %*% MUWTWUM %*% t(u))

          // A' W R in clubSand notation
          arma::mat At_WX = mat_sqrt_inv(
            I_min_H - arma::trans(H) + Xoriginal.rows(cluster_ids) * MUWTWUM * arma::trans(Xoriginal.rows(cluster_ids)),
            true
          ) * X.rows(cluster_ids) ;

          if (ci) {
            arma::mat ME = (XtX_inv / weight_mean) * arma::trans(At_WX);

            P_diags.col(clusternum) = arma::sum(arma::pow(ME, 2), 1);

            arma::mat MEU = ME * Xoriginal.rows(cluster_ids);

            H1s.slice(clusternum) = MEU * M_U_ct;
            H2s.slice(clusternum) = ME * X.rows(cluster_ids) * M_U_ct;
            H3s.slice(clusternum) = MEU * Omega_ct;
          }

          // t(ei) %*% (I - P_ss)^{-1/2} %*% Xj
          // each ro  w is the contribution of the cluster to the meat
          // Below use  t(tutX) %*% tutX to sum contributions across clusters

          tutX.row(clusternum) = arma::trans(ei(cluster_ids)) * At_WX;

          clusternum++;
        }

        Vcov_hat = XtX_inv * (arma::trans(tutX) * tutX) * XtX_inv;

        if (ci) {
          for(int p = 0; p < k; p++){
            // only compute for covars that we need the DoF for
            if (which_covs[p]) {

              // cast subcube as mat
              arma::mat H1 = H1s.subcube(arma::span(p), arma::span::all, arma::span::all);
              arma::mat H2 = H2s.subcube(arma::span(p), arma::span::all, arma::span::all);
              arma::mat H3 = H3s.subcube(arma::span(p), arma::span::all, arma::span::all);

              arma::mat uf = arma::trans(H1) * H2;

              arma::mat P_array = arma::trans(H3)*H3 - uf - arma::trans(uf);
              P_array.diag() += P_diags.row(p);

              dof[p] = std::pow(arma::trace(P_array), 2) / arma::accu(arma::pow(P_array, 2));
            }
          }
        }
      } else if (type == "stata") {

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

        dof.fill(J - 1);

      }

    }
  }

  return List::create(_["beta_hat"]= beta_hat,
                      _["Vcov_hat"]= Vcov_hat,
                      _["dof"]= dof);
}



