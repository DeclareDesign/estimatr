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
Eigen::MatrixXd Kr(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
  Eigen::MatrixXd AB(A.rows() * B.rows(), A.cols() * B.cols());

  for (int i = 0; i < A.rows(); i++) {
    for (int j = 0; j < A.cols(); j++) {
      AB.block(i*B.rows(), j*B.cols(), B.rows(), B.cols()) = A(i, j) * B;
    }
  }
  return AB;
}

// [[Rcpp::export]]
List lm_solver(const Eigen::Map<Eigen::MatrixXd>& X,
               const Eigen::Map<Eigen::MatrixXd>& y,
               const bool& try_cholesky) {

  const int p(X.cols()), ny(y.cols());
  int r = p;
  Eigen::MatrixXd XtX_inv, R_inv;
  Eigen::MatrixXd beta_out(Eigen::MatrixXd::Constant(p, ny, ::NA_REAL));

  //Rcpp::Rcout << y << std::endl;
  bool do_qr = !try_cholesky;
  if (try_cholesky) {
    const Eigen::LLT<Eigen::MatrixXd> llt(X.transpose() * X);

    // Catch case where X is rank-deficient
    if (llt.info() == Eigen::NumericalIssue) {
      do_qr = true;
    } else{
      beta_out = llt.solve(X.adjoint() * y);
      R_inv = llt.matrixL().solve(Eigen::MatrixXd::Identity(p, p));
      XtX_inv = R_inv.transpose() * R_inv;
    }
  }

  if (do_qr) {
    // Rcout << X << std::endl;

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR(X);
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType Pmat(PQR.colsPermutation());

    r = PQR.rank();

    R_inv = PQR.matrixQR().topLeftCorner(r, r).triangularView<Eigen::Upper>().solve(Eigen::MatrixXd::Identity(r, r));

    // Rcout << "R_inv:" << std::endl;
    // Rcout << R_inv << std::endl;

    Eigen::MatrixXd effects(PQR.householderQ().adjoint() * y);

    // Rcout << "effects:" << std::endl;
    // Rcout << effects << std::endl;

    beta_out.topRows(r) = R_inv * effects.topRows(r);
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

  return List::create(
    _["beta_hat"]= beta_out,
    _["XtX_inv"]= XtX_inv
  );
}

// [[Rcpp::export]]
List lm_variance(const Eigen::Map<Eigen::MatrixXd>& X,
                 const Eigen::Map<Eigen::MatrixXd>& XtX_inv,
                 const Eigen::Map<Eigen::MatrixXd>& ei,
                 const Rcpp::Nullable<Rcpp::IntegerVector> & cluster,
                 const int & J,
                 const bool & ci,
                 const String type,
                 const std::vector<bool> & which_covs) {

  const int n(X.rows()), r(X.cols()), ny(ei.cols());
  const double clustered = ((type == "stata") || (type == "CR0"));
  const int npars = r * ny;
  int sandwich_size = n;
  if (clustered) {
    sandwich_size = J;
  }

  // Rcout << "X:" << X << std::endl << std::endl;
  // Rcout << "ny: " << ny << std::endl;
  // Rcout << "ei: " << std::endl << ei << std::endl;
  Eigen::MatrixXd Vcov_hat;
  Eigen::VectorXd dof = Eigen::VectorXd::Constant(npars, -99.0);
  Eigen::MatrixXd s2 = Eigen::MatrixXd::Constant(ny, ny, -99.0);

  // Standard error calculations
  if (type == "classical") {
    // Classical
    s2 = AtA(ei)/((double)n - (double)r);
    Vcov_hat = Kr(s2, XtX_inv);

    dof.fill(n - r);

  } else {
    // Robust
    Eigen::MatrixXd temp_omega = ei.array().pow(2);

    s2 = temp_omega.colwise().sum()/((double)n - (double)r);

    dof.fill(n - r);

    Eigen::MatrixXd bread(npars, npars);
    Eigen::MatrixXd half_meat(sandwich_size, npars);
    if (ny == 1) {
      bread = XtX_inv;
    } else {
      bread = Kr(Eigen::MatrixXd::Identity(ny, ny), XtX_inv);
    }

    if ( !clustered ) {
      // Robust, no clusters

      if (type == "HC2") {

        for (int i = 0; i < n; i++) {
          Eigen::VectorXd Xi = X.row(i);
          temp_omega.row(i) = temp_omega.row(i) / (1.0 - (Xi.transpose() * XtX_inv * Xi));
        }

      } else if (type == "HC3") {

        for (int i = 0; i < n; i++) {
          Eigen::VectorXd Xi = X.row(i);
          temp_omega.row(i) = temp_omega.row(i) / (std::pow(1.0 - Xi.transpose() * XtX_inv * Xi, 2));
        }

      }

      // Rcout << "temp_omega: " << std::endl << temp_omega << std::endl;
      for (int m = 0; m < ny; m++) {
        half_meat.block(0, r*m, n, r) = X.array().colwise() * temp_omega.col(m).array().sqrt();
      }

    } else {
      // Robust, clustered

      dof.fill(J - 1);
      Eigen::Map<Eigen::ArrayXi> clusters = Rcpp::as<Eigen::Map<Eigen::ArrayXi> >(cluster);

      double current_cluster = clusters(0);
      int clust_num = 0;
      int start_pos = 0;
      int len = 1;

      // iterate over unique cluster values
      for (int i = 1; i <= n; ++i){

        if ((i == n) || (clusters(i) != current_cluster)) {
          // Rcout << current_cluster << std::endl;
          // Rcout << start_pos << std::endl;
          // Rcout << len << std::endl << std::endl;
          // Rcout <<  X.transpose().block(0, start_pos, r, len) << std::endl << std::endl;
          // Rcout << "XteetX: " << AtA(ei.segment(start_pos, len).transpose() * X.block(start_pos, 0, len, r)) << std::endl;

          if (ny > 1) {

            // Stack residuals for this cluster from each model
            Eigen::MatrixXd ei_block = ei.block(start_pos, 0, len, ny);
            Eigen::Map<const Eigen::MatrixXd> ei_long(ei_block.data(), 1, len*ny);
            // Rcout << "clust_num: " << clust_num << std::endl;
            // Rcout << "ei_long:" << std::endl << ei_long << std::endl;
            half_meat.block(clust_num, 0, 1, npars) =
              ei_long *
              Kr(Eigen::MatrixXd::Identity(ny, ny), X.block(start_pos, 0, len, r));

          } else {
            // Rcout << "clust_num: " << clust_num << std::endl;
            // Rcout << "ei:" << std::endl << ei.block(start_pos, 0, len, 1).transpose() << std::endl;
            half_meat.row(clust_num) =
              ei.block(start_pos, 0, len, 1).transpose() *
              X.block(start_pos, 0, len, r);
          }


          if (i < n) {
            current_cluster = clusters(i);
            len = 1;
            start_pos = i;
          }

          clust_num++;

        } else {
          len++;
          continue;
        }
      }
    }

    // Rcout << "bread: " << std::endl<< bread << std::endl;
    // Rcout << "half_meat:" << std::endl << half_meat << std::endl;
    Vcov_hat = bread * (half_meat.transpose() * half_meat) * bread;

  }

  if (type == "HC1") {

    Vcov_hat =
      Vcov_hat *
      (double)n / ((double)n - (double)r);

  } else if (type == "stata") {

    // Rcout << "correction: " << (((double)J * (n - 1)) / (((double)J - 1) * (n - r))) << std::endl;
    Vcov_hat =
      Vcov_hat *
      (((double)J * (n - 1)) / (((double)J - 1) * (n - r)));
  }

  return List::create(_["Vcov_hat"]= Vcov_hat,
                      _["dof"]= dof,
                      _["res_var"]= s2);
}


//  if (weighted) {
// Instead of having X and Xt be weighted by sqrt(W), this has them weighted by W
// X.array().colwise() *= weights;
// }
// [[Rcpp::export]]
List lm_variance_cr2(const Eigen::Map<Eigen::MatrixXd>& X,
                     const Rcpp::Nullable<Rcpp::NumericMatrix> & Xunweighted,
                     const Eigen::Map<Eigen::MatrixXd>& XtX_inv,
                     const Eigen::Map<Eigen::MatrixXd>& ei,
                     const double weight_mean,
                     const Eigen::Map<Eigen::ArrayXi>& clusters,
                     const int & J,
                     const bool & ci,
                     const std::vector<bool> & which_covs) {


  const int n(X.rows()), r(X.cols()), ny(ei.cols());
  const int npars = r * ny;

  // Rcout << "X:" << X << std::endl << std::endl;
  // Rcout << "beta_hat:" << beta_hat << std::endl << std::endl;

  Eigen::MatrixXd Vcov_hat;
  Eigen::VectorXd dof = Eigen::VectorXd::Constant(npars, -99.0);

  Eigen::MatrixXd Xoriginal(n, r);
  if (Xunweighted.isNotNull()) {
    Xoriginal = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Xunweighted);
  } else {
    Xoriginal = X;
  }

  // Following code implements the clubSandwich package CR2 estimator found here:
  // https://github.com/jepusto/clubSandwich, using the R code as a template for
  // implementation

  Eigen::MatrixXd bread(npars, npars);
  Eigen::MatrixXd half_meat(J, npars);
  if (ny == 1) {
    bread = XtX_inv;
  } else {
    bread = Kr(Eigen::MatrixXd::Identity(ny, ny), XtX_inv);
  }

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
  int clust_num = 0;
  int start_pos = 0;
  int len = 1;

  // Rcout << "ei: " << std::endl << ei << std::endl;

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
      Eigen::MatrixXd H =
        Xoriginal.block(start_pos, 0, len, r) *
        XtX_inv *
        X.block(start_pos, 0, len, r).transpose();

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
          (Eigen::MatrixXd::Identity(len, len) - H) - H.transpose() +
            Xoriginal.block(start_pos, 0, len, r) *
            MUWTWUM *
            Xoriginal.block(start_pos, 0, len, r).transpose()
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
        At_WX.eigenvectors() *
        eigvals.asDiagonal() *
        At_WX.eigenvectors().transpose() *
        X.block(start_pos, 0, len, r);

      if (ci) {

        Eigen::MatrixXd ME(r, len);
        if (weight_mean != 1) {
          ME = (XtX_inv / weight_mean) * At_WX_inv.transpose();
        } else {
          ME = XtX_inv * At_WX_inv.transpose();
        }

        P_diags.col(clust_num) = ME.array().pow(2).rowwise().sum();

        Eigen::MatrixXd MEU = ME * Xoriginal.block(start_pos, 0, len, r);

        int p_pos = clust_num*r;
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
      if (ny > 1) {

        // Stack residuals for this cluster from each model
        // Rcout << "len: " << len << std::endl;
        Eigen::MatrixXd ei_block = ei.block(start_pos, 0, len, ny);
        Eigen::Map<const Eigen::MatrixXd> ei_long(ei_block.data(), 1, len*ny);
        // Rcout << "clust_num: " << clust_num << std::endl;
        // Rcout << "ei.block" << ei.block(start_pos, 0, len, ny) << std::endl;
        // Rcout << "ei_long:" << std::endl << ei_long << std::endl;
        half_meat.block(clust_num, 0, 1, npars) =
          ei_long *
          Kr(Eigen::MatrixXd::Identity(ny, ny), At_WX_inv);

      } else {
        // Rcout << "clust_num: " << clust_num << std::endl;
        // Rcout << "ei:" << std::endl << ei.block(start_pos, 0, len, 1).transpose() << std::endl;
        half_meat.row(clust_num) =
          ei.block(start_pos, 0, len, 1).transpose() *
          At_WX_inv;
      }

      if (i < n) {
        current_cluster = clusters(i);
        len = 1;
        start_pos = i;
        clust_num++;
      }

    } else {
      len++;
      continue;
    }

  }

  // Rcout << "bread: " << std::endl << bread << std::endl;
  // Rcout << "half_meat: " << std::endl << half_meat << std::endl;
  // Rcout << "meat: " << std::endl << (half_meat.transpose() * half_meat) << std::endl;

  Vcov_hat = bread * (half_meat.transpose() * half_meat) * bread;

  // Rcout << "H1s" << std::endl << std::endl;
  // Rcout << H1s << std::endl << std::endl;
  // Rcout << "H2s" << std::endl << std::endl;
  // Rcout << H2s << std::endl << std::endl;
  // Rcout << "H3s" << std::endl << std::endl;
  // Rcout << H3s << std::endl << std::endl;

  if (ci) {
    for(int j = 0; j < r; j++){
      if (which_covs[j]) {

        Eigen::MatrixXd H1t = H1s.row(j);
        Eigen::MatrixXd H2t = H2s.row(j);
        Eigen::MatrixXd H3t = H3s.row(j);
        // Rcout << H1t << std::endl<< std::endl;

        H1t.resize(r, J);
        H2t.resize(r, J);
        H3t.resize(r, J);

        // Rcout << H1t << std::endl<< std::endl;

        Eigen::MatrixXd uf = H1t.transpose() * H2t;
        Eigen::MatrixXd P_row = P_diags.row(j).asDiagonal();
        Eigen::MatrixXd P_array = (H3t.transpose()*H3t - uf - uf.transpose()) + P_row;
        // Rcout << "P_array: " << P_array << std::endl;

        // Rcout << std::pow(P_array.trace(), 2) << " and " << P_array.array().pow(2).sum() << std::endl;
        dof(j) = std::pow(P_array.trace(), 2) / P_array.array().pow(2).sum();
      }
    }
  }

  return List::create(
    _["Vcov_hat"]= Vcov_hat,
    _["dof"]= dof
  );
}
