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
arma::mat mat_sqrt_inv(const arma::mat & X) {
  arma::vec eigvals;
  arma::mat eigvecs;
  arma::eig_sym(eigvals, eigvecs, X);
  return(mult_diag(eigvecs, 1.0 / arma::sqrt(eigvals)) * arma::trans(eigvecs));
}

//----------
// Main LM fit and SE function
//----------
// [[Rcpp::export]]
List lm_robust_helper(const arma::vec & y,
                      const arma::mat & X,
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
    Rcpp::Rcout << "weights" << std::endl;
    weights = Rcpp::as<arma::vec>(weight);
    Xoriginal = Rcpp::as<arma::mat>(Xunweighted);
  }

  // t(X)
  arma::mat Xt = arma::trans(X);

  //XtX <- t(X_mat) %*% X_mat
  arma::mat XtX = Xt*X;

  // xtx inverse
  arma::mat XtX_inv = arma::inv(XtX);

  // beta_hat <- solve(XtX) %*% t(X_mat) %*% Y
  arma::colvec beta_hat = XtX_inv*Xt*y;

  arma::mat Vcov_hat;
  arma::colvec dof(X.n_cols);
  dof.fill(-99);

  if(type != "none"){

    Rcpp::Rcout << y << std::endl;

    Rcpp::Rcout << y * sqrt(weight_mean)<< std::endl;
    Rcpp::Rcout << X << std::endl;

    Rcpp::Rcout << X * sqrt(weight_mean) << std::endl;
    // residuals
    arma::colvec ei = y - X * beta_hat;
    Rcpp::Rcout << ei << std::endl;

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
    } else if (( (type == "BM") | (type == "stata") | (type == "CR2") ) & cluster.isNotNull()) {

      // Code adapted from Michal Kolesar
      // https://github.com/kolesarm/Robust-Small-Sample-Standard-Errors

      // Get unique cluster values
      arma::vec clusters = Rcpp::as<arma::vec>(cluster);
      arma::vec levels = unique(clusters);
      double J = levels.n_elem;

      if (type == "CR2") {

        // if weights
        Xt = Xt.each_row() % arma::sqrt(arma::trans(weights));

        arma::mat tutX(J, k);

        // Cube stores the n by k G_s matrix for each cluster
        // used for the BM dof corrction
        arma::cube Gs(n, k, J);

        ei = y / sqrt(weights) - (X.each_col() / sqrt(weights)) * beta_hat;
        Rcpp::Rcout << ei << std::endl;
        // iterator used to fill tutX
        int clusternum = 0;

        Rcpp::Rcout << XtX_inv << std::endl;
        arma::mat MUWTWUM = XtX_inv * Xt * arma::trans(Xt) * XtX_inv;

        // iterate over unique cluster values
        for(arma::vec::const_iterator j = levels.begin();
            j != levels.end();
            ++j){

          arma::uvec cluster_ids = find(clusters == *j);
          int cluster_size = cluster_ids.n_elem;
        //    H_jj <- Map(function(u, uw) u %*% M_U %*% uw,
        //u = U_list, uw = UW_list)

          arma::mat D = arma::eye(cluster_size, cluster_size);
          arma::mat H = Xoriginal.rows(cluster_ids) * XtX_inv * Xt.cols(cluster_ids);
          arma::mat I_min_H = D - H;
         // uwTwu <- Map(function(uw, th) uw %*% th %*% t(uw),
          //             uw = UW_list, th = Theta_list)
           // MUWTWUM <- M_U %*% Reduce("+", uwTwu) %*% M_U

          //(thet - h %*% thet - thet %*% t(h) + u %*% MUWTWUM %*% t(u))

          arma::mat A = D * arma::sqrtmat_sympd(arma::pinv(
            I_min_H - arma::trans(H) + Xoriginal.rows(cluster_ids) * MUWTWUM * arma::trans(Xoriginal.rows(cluster_ids))
          )) * D * arma::trans(Xt.cols(cluster_ids)) ;

          //arma::mat A = D * arma::sqrtmat_sympd(arma::pinv(
          //  D * (I_min_H) * D * D
          //)) * D * X.rows(cluster_ids) ;

          Rcpp::Rcout << A << std::endl;


          // t(ei) %*% (I - P_ss)^{-1/2} %*% Xj
          // each ro  w is the contribution of the cluster to the meat
          // Below use  t(tutX) %*% tutX to sum contributions across clusters

          tutX.row(clusternum) = arma::trans(ei(cluster_ids))  * A;

          clusternum++;
        }

        Rcpp::Rcout << tutX << std::endl;
        Rcpp::Rcout << weight_mean << std::endl;
        Vcov_hat = XtX_inv * (arma::trans(tutX) * tutX) * XtX_inv;

        Rcpp::Rcout <<  (arma::trans(tutX) * tutX)<< std::endl;

        if (ci) {
          for(int p = 0; p < k; p++){
            // only compute for covars that we need the DoF for
            if (which_covs[p]) {
              arma::mat G = Gs.subcube(arma::span::all, arma::span(p), arma::span::all);
              arma::mat GG = arma::trans(G) * G;
              dof[p] = std::pow(arma::trace(GG), 2) / arma::accu(arma::pow(GG, 2));
            }
          }
        }
      }

      if (type == "BM") {

        arma::mat tutX(J, k);

        // Cube stores the n by k G_s matrix for each cluster
        // used for the BM dof corrction
        arma::cube Gs(n, k, J);

        // iterator used to fill tutX
        int clusternum = 0;

        // iterate over unique cluster values
        for(arma::vec::const_iterator j = levels.begin();
            j != levels.end();
            ++j){

          arma::uvec cluster_ids = find(clusters == *j);
          int cluster_size = cluster_ids.n_elem;

          arma::mat XtX_inv_tXj = XtX_inv * Xt.cols(cluster_ids);

          // (I_j - Xj %*% (X'X)^{-1} %*% Xj') ^ {-1/2} %*% Xj
          // i.e. (I_j - P_ss) ^ {-1/2} %*% Xj
          arma::mat minP_Xt = mat_sqrt_inv(
            arma::eye(cluster_size, cluster_size) -
              X.rows(cluster_ids) * XtX_inv_tXj
          ) * X.rows(cluster_ids);

          // t(ei) %*% (I - P_ss)^{-1/2} %*% Xj
          // each row is the contribution of the cluster to the meat
          // Below use  t(tutX) %*% tutX to sum contributions across clusters
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
            // only compute for covars that we need the DoF for
            if (which_covs[p]) {
              arma::mat G = Gs.subcube(arma::span::all, arma::span(p), arma::span::all);
              arma::mat GG = arma::trans(G) * G;
              dof[p] = std::pow(arma::trace(GG), 2) / arma::accu(arma::pow(GG, 2));
            }
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

        dof.fill(J - 1);

      }

    }
  }

  return List::create(_["beta_hat"]= beta_hat,
                      _["Vcov_hat"]= Vcov_hat,
                      _["dof"]= dof);
}



