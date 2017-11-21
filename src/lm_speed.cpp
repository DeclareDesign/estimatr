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
Eigen::MatrixXd mult_diag(const Eigen::MatrixXd& x, const Eigen::ArrayXd& d) {

  Eigen::MatrixXd out(x.rows(), x.cols());
  for (unsigned j = 0; j < x.cols(); ++j) {
    out.col(j) = x.col(j) * d(j);
  }

  return out;
}

//' @export
// [[Rcpp::export]]
List lm_ei_test(Eigen::Map<Eigen::MatrixXd>& X,
                const Eigen::Map<Eigen::VectorXd>& y,
                const Rcpp::Nullable<Rcpp::NumericMatrix> & Xunweighted,
                const Rcpp::Nullable<Rcpp::NumericVector> & weight,
                const double & weight_mean,
                const Rcpp::Nullable<Rcpp::NumericVector> & cluster,
                const bool & ci,
                const String type,
                const std::vector<bool> & which_covs,
                const bool& chol,
                const bool& trychol) {

  const int n(X.rows()), p(X.cols());
  const Eigen::MatrixXd Xt(X.transpose());
  Eigen::MatrixXd XtX_inv, R_inv, Vcov_hat;
  Eigen::VectorXd beta_hat(Eigen::VectorXd::Constant(p, ::NA_REAL));
  Eigen::VectorXd dof(p);

  bool full_rank = true;

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
          full_rank = false;
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

  if(type != "none"){

    // residuals
    Eigen::VectorXd ei = y - X * beta_hat;

    if (type == "classical") {

      // the next line due to Dirk Eddelbuettel
      // http://dirk.eddelbuettel.com/code/rcpp.armadillo.html
      double s2 = ei.dot(ei)/(n - p);
      Vcov_hat = s2 * XtX_inv;

    } else if ( (type == "HC0") | (type == "HC1") | (type == "HC2") | (type == "HC3") ) {

      // hat values
      Eigen::ArrayXd ei2 = ei.array().pow(2);

      if(type == "HC0"){
        Vcov_hat = XtX_inv * mult_diag(Xt, ei2) * X * XtX_inv;

      } else if (type == "HC1") {
        Vcov_hat = n/(n-p) * XtX_inv * mult_diag(Xt, ei2) * X * XtX_inv;
      } else if (type == "HC2") {
        Eigen::ArrayXd hii(n);

        for(int i = 0; i < n; i++){
          Eigen::VectorXd Xi = X.row(i);
          hii(i) = 1.0 - (Xi.transpose() * XtX_inv * Xi);
        }
        Vcov_hat = XtX_inv * mult_diag(Xt, ei2 / hii) * X * XtX_inv;
      } else if (type == "HC3") {
        Eigen::ArrayXd hii(n);

        for(int i = 0; i < n; i++){
          Eigen::VectorXd Xi = X.row(i);
          hii(i) = std::pow(1.0 - Xi.transpose() * XtX_inv * Xi, 2);
        }
        Vcov_hat = XtX_inv * mult_diag(Xt, ei2 / hii) * X * XtX_inv;
      }

//    } //else if ( (type == "stata") | (type == "CR2") ) {
//
//       // Code adapted from Michal Kolesar
//       // https://github.com/kolesarm/Robust-Small-Sample-Standard-Errors
//
//       // Get unique cluster values
//       arma::vec clusters = Rcpp::as<arma::vec>(cluster);
//       arma::vec levels = unique(clusters);
//       double J = levels.n_elem;
//
//       //Rcpp::Rcout << 'here' << std::endl;
//
//       if (type == "CR2") {
//
//         arma::mat Xoriginal(X.n_rows, X.n_cols);
//         arma::vec weights;
//         bool weighted;
//         if(Xunweighted.isNotNull()) {
//           weights = Rcpp::as<arma::vec>(weight);
//           Xoriginal = Rcpp::as<arma::mat>(Xunweighted);
//           weighted = true;
//         } else {
//           Xoriginal = X;
//           weighted = false;
//         }
//
//         if (weighted) {
//           // Instead of having X and Xt be weighted by sqrt(W), this has them weighted by W
//           Xt = Xt.each_row() % arma::trans(weights);
//           X = X.each_col() % weights;
//           // Unweighted residuals
//           ei = y / weights - Xoriginal * beta_hat;
//         }
//
//         arma::mat tutX(J, k);
//
//         // used for the dof corrction
//         arma::cube H1s(k, k, J);
//         arma::cube H2s(k, k, J);
//         arma::cube H3s(k, k, J);
//         arma::mat P_diags(k, J);
//
//         arma::mat M_U_ct = arma::chol(XtX_inv, "lower");
//
//         int clusternum = 0;
//
//         arma::mat MUWTWUM = XtX_inv * Xt * arma::trans(Xt) * XtX_inv;
//         arma::mat Omega_ct = arma::chol(MUWTWUM, "lower");
//
//         // iterate over unique cluster values
//         for(arma::vec::const_iterator j = levels.begin();
//             j != levels.end();
//             ++j){
//
//           arma::uvec cluster_ids = find(clusters == *j);
//           int cluster_size = cluster_ids.n_elem;
//
//           arma::mat D = arma::eye(cluster_size, cluster_size);
//           arma::mat H = Xoriginal.rows(cluster_ids) * XtX_inv * Xt.cols(cluster_ids);
//           arma::mat I_min_H = D - H;
//
//           // Code from clubSandwich
//           // uwTwu <- Map(function(uw, th) uw %*% th %*% t(uw),
//           //             uw = UW_list, th = Theta_list)
//           // MUWTWUM <- M_U %*% Reduce("+", uwTwu) %*% M_U
//
//           //(thet - h %*% thet - thet %*% t(h) + u %*% MUWTWUM %*% t(u))
//
//           // A' W R in clubSand notation
//           arma::mat At_WX = mat_sqrt_inv(
//             I_min_H - arma::trans(H) + Xoriginal.rows(cluster_ids) * MUWTWUM * arma::trans(Xoriginal.rows(cluster_ids)),
//             true
//           ) * X.rows(cluster_ids) ;
//
//           if (ci) {
//
//             arma::mat ME(k, cluster_size);
//             if (weighted) {
//               ME = (XtX_inv / weight_mean) * arma::trans(At_WX);
//             } else {
//               ME = XtX_inv * arma::trans(At_WX);
//             }
//
//             P_diags.col(clusternum) = arma::sum(arma::pow(ME, 2), 1);
//
//             arma::mat MEU = ME * Xoriginal.rows(cluster_ids);
//
//             H1s.slice(clusternum) = MEU * M_U_ct;
//             H2s.slice(clusternum) = ME * X.rows(cluster_ids) * M_U_ct;
//             H3s.slice(clusternum) = MEU * Omega_ct;
//           }
//
//           // t(ei) %*% (I - P_ss)^{-1/2} %*% Xj
//           // each ro  w is the contribution of the cluster to the meat
//           // Below use  t(tutX) %*% tutX to sum contributions across clusters
//
//           tutX.row(clusternum) = arma::trans(ei(cluster_ids)) * At_WX;
//
//           clusternum++;
//         }
//
//         Vcov_hat = XtX_inv * (arma::trans(tutX) * tutX) * XtX_inv;
//
//         if (ci) {
//           for(int p = 0; p < k; p++){
//             // only compute for covars that we need the DoF for
//             if (which_covs[p]) {
//
//               // cast subcube as mat
//               arma::mat H1 = H1s.subcube(arma::span(p), arma::span::all, arma::span::all);
//               arma::mat H2 = H2s.subcube(arma::span(p), arma::span::all, arma::span::all);
//               arma::mat H3 = H3s.subcube(arma::span(p), arma::span::all, arma::span::all);
//
//               arma::mat uf = arma::trans(H1) * H2;
//
//               arma::mat P_array = arma::trans(H3)*H3 - uf - arma::trans(uf);
//               P_array.diag() += P_diags.row(p);
//
//               dof[p] = std::pow(arma::trace(P_array), 2) / arma::accu(arma::pow(P_array, 2));
//             }
//           }
//         }
//       } else if (type == "stata") {
//
//         arma::mat XteetX(k, k);
//         XteetX.fill(0.0);
//
//         // iterate over unique cluster values
//         for(arma::vec::const_iterator j = levels.begin();
//             j != levels.end();
//             ++j){
//           arma::uvec cluster_ids = find(clusters == *j);
//
//           XteetX += Xt.cols(cluster_ids) * ei.elem(cluster_ids) *
//             arma::trans(ei.elem(cluster_ids)) * X.rows(cluster_ids);
//         }
//
//         Vcov_hat = ((J * (n - 1)) / ((J - 1) * (n - k))) * XtX_inv * XteetX * XtX_inv;
//
//         dof.fill(J - 1);
//
//       }
//
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
