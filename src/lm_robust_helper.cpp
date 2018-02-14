// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppEigen.h>
using namespace Rcpp;

// Much of what follows is modified from RcppEigen Vignette by Douglas Bates and Dirk Eddelbuettel
// https://cran.r-project.org/web/packages/RcppEigen/vignettes/RcppEigen-Introduction.pdf
// [[Rcpp::export]]
Eigen::MatrixXd AtA(const Eigen::MatrixXd& A) {
  int n(A.cols());
  return Eigen::MatrixXd(n,n).setZero().selfadjointView<Eigen::Lower>()
                             .rankUpdate(A.adjoint());
}

// [[Rcpp::export]]
List lm_solver(Eigen::Map<Eigen::MatrixXd>& Xfull,
               const Eigen::Map<Eigen::VectorXd>& y,
               const Rcpp::Nullable<Rcpp::NumericMatrix> & Xunweighted,
               const Rcpp::Nullable<Rcpp::NumericVector> & weight,
               const double & weight_mean,
               const Rcpp::Nullable<Rcpp::IntegerVector> & cluster,
               const int & J,
               const bool & ci,
               const String type,
               const std::vector<bool> & which_covs,
               const bool& try_cholesky) {

  const int n(Xfull.rows()), p(Xfull.cols());
  int r = p;
  Eigen::MatrixXd Xfullt(Xfull.transpose());
  Eigen::MatrixXd XtX_inv, R_inv, Vcov_hat;
  Eigen::VectorXd beta_out(Eigen::VectorXd::Constant(p, ::NA_REAL));


  bool do_qr = !try_cholesky;
  if (try_cholesky) {
    const Eigen::LLT<Eigen::MatrixXd> llt(Xfullt*Xfull);

    // Catch case where Xfull is rank-deficient
    if (llt.info() == Eigen::NumericalIssue) {
      do_qr = true;
    } else{
      beta_out = llt.solve(Xfull.adjoint() * y);
      R_inv = llt.matrixL().solve(Eigen::MatrixXd::Identity(p, p));
      XtX_inv = R_inv.transpose() * R_inv;
    }

  }

  if (do_qr) {
    // Rcout << Xfull << std::endl;

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR(Xfull);
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType Pmat(PQR.colsPermutation());

    r = PQR.rank();

    R_inv = PQR.matrixQR().topLeftCorner(r, r).triangularView<Eigen::Upper>().solve(Eigen::MatrixXd::Identity(r, r));

    // Rcout << "R_inv:" << std::endl;
    // Rcout << R_inv << std::endl;

    Eigen::VectorXd effects(PQR.householderQ().adjoint() * y);

    // Rcout << "effects:" << std::endl;
    // Rcout << effects << std::endl;

    beta_out.head(r) = R_inv * effects.head(r);
    // Rcout << "beta_out:" << std::endl;
    // Rcout << beta_out << std::endl;
    beta_out = Pmat * beta_out;

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
    // Rcout << "Pmat * I" << std::endl;
    // Rcout << Pmat * Eigen::MatrixXd::Identity(r, r) << std::endl << std::endl;

    R_inv = P * R_inv * P;
    XtX_inv = R_inv * R_inv.transpose();
    // Rcout << "XtX_inv: " << std::endl;
    // Rcout << XtX_inv << std::endl;
  }

  // TODO See whether beta was not fit (i.e. dependent columns of X)
  //Eigen::ArrayBase<Derived> = beta_out.array().isNaN();

  // Two copies of the data because both weighted and unweighted needed for CR2
  Eigen::MatrixXd Xoriginal(n, r);
  Eigen::MatrixXd Xoriginalfull(n, p);
  Eigen::MatrixXd X(n, r);
  Eigen::VectorXd beta_hat(r);
  Eigen::ArrayXd weights(n);
  bool weighted = Xunweighted.isNotNull();
  if (weighted) {
    weights = Rcpp::as<Eigen::Map<Eigen::ArrayXd> >(weight);
    Xoriginalfull = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Xunweighted);
  }
  // Reorder objects if rank-deficient
  if (r < p) {
    // Rcout << "Rank-deficient!" << std::endl;
    // Rcout << "beta_out:" << beta_out << std::endl;

    int j = 0;
    for (int i = 0; i < p; i++) {
      // Rcout << "isnan: " << std::isnan(beta_out(i)) << std::endl;
      if (!std::isnan(beta_out(i))) {
        // Rcout << "not NA: " << i << std::endl;
        beta_hat(j) = beta_out(i);
        X.col(j) = Xfull.col(i);
        if (weighted) {
          Xoriginal.col(j) = Xoriginalfull.col(i);
        } else {
          Xoriginal.col(j) = Xfull.col(i);
        }
        j++;
      }
    }
  } else {
    beta_hat = beta_out;
    X = Xfull;
    if (weighted) {
      Xoriginal = Xoriginalfull;
    } else {
      Xoriginal = Xfull;
    }
  }

  // Rcout << "X:" << X << std::endl << std::endl;
  // Rcout << "Xoriginal:" << Xoriginal << std::endl << std::endl;
  // Rcout << "beta_hat:" << beta_hat << std::endl << std::endl;

  Eigen::VectorXd dof = Eigen::VectorXd::Constant(r, -99.0);
  Eigen::VectorXd ei = Eigen::VectorXd::Constant(n, -99.0);
  double s2 = -99.0;

  // Standard error calculations
  if (type != "none"){

    // residuals
    ei = y - X * beta_hat;
    s2 = ei.dot(ei)/(n - r);

    // Rcout << ei << std::endl;

    if (type == "classical") {

      Vcov_hat = s2 * XtX_inv;

    } else if ( (type == "HC0") | (type == "HC1") | (type == "HC2") | (type == "HC3") ) {

      // hat values
      Eigen::ArrayXd ei2 = ei.array().pow(2);

      if(type == "HC0"){

        Vcov_hat = XtX_inv * (X.transpose() * ei2.matrix().asDiagonal()) * X * XtX_inv;

      } else if (type == "HC1") {

        Vcov_hat =
          (double)n / ((double)n - (double)r) *
          XtX_inv * (X.transpose() * ei2.matrix().asDiagonal()) * X * XtX_inv;

      } else if (type == "HC2") {
        Eigen::ArrayXd hii(n);

        for (int i = 0; i < n; i++) {
          Eigen::VectorXd Xi = X.row(i);
          hii(i) = ei2(i) / (1.0 - (Xi.transpose() * XtX_inv * Xi));
        }
        Vcov_hat = XtX_inv * (X.transpose() * hii.matrix().asDiagonal()) * X * XtX_inv;
      } else if (type == "HC3") {
        Eigen::ArrayXd hii(n);

        for (int i = 0; i < n; i++) {
          Eigen::VectorXd Xi = X.row(i);
          hii(i) = ei2(i) / (std::pow(1.0 - Xi.transpose() * XtX_inv * Xi, 2));
        }
        Vcov_hat = XtX_inv * (X.transpose() * hii.matrix().asDiagonal()) * X * XtX_inv;
      }

    } else if ( (type == "stata") || (type == "CR2") || (type == "CR0") ) {

      Eigen::Map<Eigen::ArrayXi> clusters = Rcpp::as<Eigen::Map<Eigen::ArrayXi> >(cluster);

      // Following code implements the clubSandwich package CR2 estimatr found here:
      // https://github.com/jepusto/clubSandwich
      // Much of the code is also a translation/adaptation from the R code in that package
      if (type == "CR2") {

        Eigen::VectorXd cr2_eis;
        if (weighted) {
          // Instead of having X and Xt be weighted by sqrt(W), this has them weighted by W
          X.array().colwise() *= weights;
          // Unweighted residuals
          cr2_eis = y.array() / weights - (Xoriginal * beta_hat).array();
        } else {
          cr2_eis = ei;
        }

        Eigen::MatrixXd tutX(J, r);

        // used for the dof corrction
        Eigen::MatrixXd H1s(r, r*J);
        Eigen::MatrixXd H2s(r, r*J);
        Eigen::MatrixXd H3s(r, r*J);
        Eigen::MatrixXd P_diags(r, J);

        Eigen::MatrixXd M_U_ct = XtX_inv.llt().matrixL();

        Eigen::MatrixXd MUWTWUM = XtX_inv * X.transpose() * X * XtX_inv;
        Eigen::MatrixXd Omega_ct = MUWTWUM.llt().matrixL();

        // Rcout << Omega_ct << std::endl;

        // Rcout << "SIZES" << std::endl;
        // Rcout << "X: " << X.rows() << "x" << X.cols() << std::endl;
        // Rcout << "Xoriginal: " << Xoriginal.rows() << "x" << Xoriginal.cols() << std::endl;
        // Rcout << "XtX_inv: " << XtX_inv.rows() << "x" << XtX_inv.cols() << std::endl;

        double current_cluster = clusters(0);
        int j = 0;
        int start_pos = 0;
        int len = 1;

        // iterate over unique cluster values
        for(int i = 1; i <= n; ++i){

          //Rcout << "clusters(i): " << clusters(i) << std::endl;
          if ((i == n) || (clusters(i) != current_cluster)) {
            // Rcout << "Current cluster: " << current_cluster << std::endl;
            // Rcout << "Starting position: " << start_pos << std::endl;
            // Rcout << "len: " << len << std::endl << std::endl;
            // Rcout <<  X.transpose().block(0, start_pos, r, len) << std::endl << std::endl;

            // Rcout <<  X.transpose().block(0, start_pos, r, len) << std::endl << std::endl;

            // TODO H should be symmetric, shouldn't need to transpose
            Eigen::MatrixXd H = Xoriginal.block(start_pos, 0, len, r) * XtX_inv * X.block(start_pos, 0, len, r).transpose();

            // Rcout << "H: " << H << std::endl;

            // Code from clubSandwich
            // uwTwu <- Map(function(uw, th) uw %*% th %*% t(uw),
            //             uw = UW_list, th = Theta_list)
            // MUWTWUM <- M_U %*% Reduce("+", uwTwu) %*% M_U

            //(thet - h %*% thet - thet %*% t(h) + u %*% MUWTWUM %*% t(u))

            // A' W R in clubSand notation
            // Rcout << "At_WX: " << (Eigen::MatrixXd::Identity(len, len) - H) - H.transpose() + Xoriginal.block(start_pos, 0, len, r) * MUWTWUM * Xoriginal.block(start_pos, 0, len, r).transpose() << std::endl;
            // Rcout << "MUWTWUM: " << MUWTWUM << std::endl;

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> At_WX(
                (Eigen::MatrixXd::Identity(len, len) - H) - H +
                  Xoriginal.block(start_pos, 0, len, r) * MUWTWUM * Xoriginal.block(start_pos, 0, len, r).transpose()
            );

            Eigen::VectorXd eigvals = At_WX.eigenvalues();
            for (int m = 0; m < eigvals.size(); ++m) {
              if (eigvals(m) > std::pow(10.0, -12.0)) {
                eigvals(m) = 1.0 / std::sqrt(eigvals(m));
              } else {
                eigvals(m) = 0;
              }
            }

            // Rcout << "Xblock: " << X.block(start_pos, 0, len, r) << std::endl;
            // Rcout << "At_WX_inv: " << At_WX.operatorInverseSqrt() << std::endl;

            // TODO mimic tol in arma implementation
            //Eigen::MatrixXd At_WX_inv = At_WX.operatorInverseSqrt() * X.block(start_pos, 0, len, r);
            Eigen::MatrixXd At_WX_inv =
              ((At_WX.eigenvectors() * eigvals.asDiagonal()) * At_WX.eigenvectors().transpose())
              * X.block(start_pos, 0, len, r);

            if (ci) {

              Eigen::MatrixXd ME(r, len);
              if (weighted) {
                ME = (XtX_inv / weight_mean) * At_WX_inv.transpose();
              } else {
                ME = XtX_inv * At_WX_inv.transpose();
              }

              P_diags.col(j) = ME.array().pow(2).rowwise().sum();

              Eigen::MatrixXd MEU = ME * Xoriginal.block(start_pos, 0, len, r);

              int p_pos = j*r;
              // Rcout << "p_pos: " << p_pos << std::endl;
              H1s.block(0, p_pos, r, r) = MEU * M_U_ct;
              H2s.block(0, p_pos, r, r) = ME * X.block(start_pos, 0, len, r) * M_U_ct;
              H3s.block(0, p_pos, r, r) = MEU * Omega_ct;
            }

            // t(cr2_eis) %*% (I - P_ss)^{-1/2} %*% Xj
            // each ro  w is the contribution of the cluster to the meat
            // Below use  t(tutX) %*% tutX to sum contributions across clusters

            // Rcout << "At_WX_inv: " << At_WX_inv << std::endl;
            // Rcout << "cr2_eis: " << cr2_eis.segment(start_pos, len).transpose() << std::endl;
            tutX.row(j) = cr2_eis.segment(start_pos, len).transpose() * At_WX_inv;

            if (i < n) {
              current_cluster = clusters(i);
              len = 1;
              start_pos = i;
              j++;
            }


          } else {
            len++;
            continue;
          }

        }


        // Rcout << "XtX_inv: " << std::endl << XtX_inv << std::endl;
        // Rcout << "tutX: " << std::endl << tutX << std::endl;

        Vcov_hat = XtX_inv * (tutX.transpose() * tutX) * XtX_inv;

        // Rcout << "H1s" << std::endl << std::endl;
        // Rcout << H1s << std::endl << std::endl;
        // Rcout << "H2s" << std::endl << std::endl;
        // Rcout << H2s << std::endl << std::endl;
        // Rcout << "H3s" << std::endl << std::endl;
        // Rcout << H3s << std::endl << std::endl;

        if (ci) {
          dof.fill(-99);
          int k = 0;
          for(int j = 0; j < p; j++){
            // only compute for covars that we need the DoF for
            if (!std::isnan(beta_out(j))) {
              if (which_covs[j]) {

                Eigen::MatrixXd H1t = H1s.row(k);
                Eigen::MatrixXd H2t = H2s.row(k);
                Eigen::MatrixXd H3t = H3s.row(k);
                // Rcout << H1t << std::endl<< std::endl;

                H1t.resize(r, J);
                H2t.resize(r, J);
                H3t.resize(r, J);

                // Rcout << H1t << std::endl<< std::endl;

                Eigen::MatrixXd uf = H1t.transpose() * H2t;
                Eigen::MatrixXd P_row = P_diags.row(k).asDiagonal();
                Eigen::MatrixXd P_array = (H3t.transpose()*H3t - uf - uf.transpose()) + P_row;
                // Rcout << "P_array: " << P_array << std::endl;

                // Rcout << std::pow(P_array.trace(), 2) << " and " << P_array.array().pow(2).sum() << std::endl;
                dof(k) = std::pow(P_array.trace(), 2) / P_array.array().pow(2).sum();
              }
              k++;
            }
          }
        }
      } else if ((type == "stata") || (type == "CR0")) {

        Eigen::MatrixXd XteetX = Eigen::MatrixXd::Zero(r, r);
        double current_cluster = clusters(0);
        int start_pos = 0;
        int len = 1;

        // iterate over unique cluster values
        for(int i = 1; i <= n; ++i){

          if ((i == n) || (clusters(i) != current_cluster)) {
            // Rcout << current_cluster << std::endl;
            // Rcout << start_pos << std::endl;
            // Rcout << len << std::endl << std::endl;
            // Rcout <<  X.transpose().block(0, start_pos, r, len) << std::endl << std::endl;
            // Rcout << "XteetX: " << AtA(ei.segment(start_pos, len).transpose() * X.block(start_pos, 0, len, r)) << std::endl;
            XteetX += AtA(ei.segment(start_pos, len).transpose() * X.block(start_pos, 0, len, r));

            if (i < n) {
              current_cluster = clusters(i);
              len = 1;
              start_pos = i;
            }

          } else {
            len++;
            continue;
          }
          //Eigen::VectorXi cluster_ids = (clusters == levels[i]).cast<int>();

          //Eigen::VectorXd ei_clus = extract_vec(ei, cluster_ids);
          //Eigen::MatrixXd X_clus = extract_mat_rows(X, cluster_ids);

        }
        // Rcout << "XteetX: " << XteetX << std::endl;
        //
        // Rcout << "J: " << J << std::endl;
        // Rcout << "n: " << n << std::endl;
        // Rcout << "r: " << r << std::endl;
        // Rcout << "corr: " << ((J * (n - 1)) / ((J - 1) * (n - r))) << std::endl;
        // Rcout << "sandwich: " << (XtX_inv * XteetX * XtX_inv) << std::endl;

        if (type == "stata") {
          Vcov_hat =
            (((double)J * (n - 1)) / (((double)J - 1) * (n - r))) *
            (XtX_inv * XteetX * XtX_inv);
        } else {
          Vcov_hat = XtX_inv * XteetX * XtX_inv;
        }

        dof.fill(J - 1);

        //Rcpp::Rcout << 'here' << std::endl;
      }
    }
  }

  // Rcout << beta_out << std::endl;

  return List::create(_["beta_hat"]= beta_out,
                      _["Vcov_hat"]= Vcov_hat,
                      _["dof"]= dof,
                      _["res_var"]= s2,
                      _["XtX_inv"]= XtX_inv,
                      _["residuals"]= ei);
}
