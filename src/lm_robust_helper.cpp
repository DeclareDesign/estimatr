// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export]]
Eigen::MatrixXd eigenAve(const Eigen::ArrayXd& x,
                         const Eigen::VectorXi& fe,
                         const Eigen::VectorXd& weights) {

  std::unordered_map<int, Eigen::Array2d> aves;
  Eigen::ArrayXd avevec(x.rows());

  for (Eigen::Index i=0; i<fe.rows(); i++) {
    Eigen::Array2d dat;
    dat(0) = weights(i) * x(i);
    dat(1) = weights(i);
    if (aves.find(fe(i)) != aves.end()) {
      aves[fe(i)] += dat;
    } else {
      aves[fe(i)] = dat;
    }
  }

  for (Eigen::Index i=0; i<fe.rows(); i++) {
    avevec(i) = x(i) - aves[fe(i)](0)/aves[fe(i)](1);
  }
  return avevec;
}

// Modifiedrom user raahlb on Stack Overflow
// https://stackoverflow.com/a/46303314/4257918
// [[Rcpp::export]]
void removeColumn(Eigen::Map<Eigen::MatrixXd>& matrix, unsigned int colToRemove) {
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;

  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);

}

// [[Rcpp::export]]
List demeanMat(const Eigen::MatrixXd& Y,
               const Eigen::MatrixXd& X,
               const Eigen::MatrixXi& fes,
               const Eigen::VectorXd& weights,
               const bool& has_int,
               const double& eps) {

  int start_col = 0 + has_int;
  // Rcout << start_col << std::endl;

  int n = X.rows();
  int p = X.cols();
  int ny = Y.cols();
  // Drop integer
  Eigen::MatrixXd newX(n, p - start_col);
  Eigen::MatrixXd newY(n, ny);

  // Rcout << X.rows() << std::endl;
  // Rcout << X.cols() << std::endl;
  // Iterate over columns of X, starting at 1 if there is an intercept
  // and then do Y
  for (Eigen::Index i = start_col; i < (p + ny); ++i) {
    Eigen::ArrayXd oldcol(n);
    Eigen::ArrayXd newcol(n);
    if (i < p) {
      oldcol = X.col(i).array() - 1.0;
      newcol = X.col(i).array();
    } else {
      oldcol = Y.col(i-p).array() - 1.0;
      newcol = Y.col(i-p).array();
    }

    while (std::sqrt((oldcol - newcol).pow(2).sum()) >= eps) {
      oldcol = newcol;
      for (Eigen::Index j = 0; j < fes.cols(); ++j) {
        newcol = eigenAve(newcol.matrix(), fes.col(j), weights);
      }
      // Rcout << "oldcol" << std::endl << oldcol << std::endl;
      // Rcout << "newcol" << std::endl << newcol << std::endl;
      // Rcout << std::sqrt((oldcol - newcol).pow(2).sum()) << std::endl;
    }
    if (i < p) {
      newX.col(i - start_col) = newcol;
    } else {
      newY.col(i - p) = newcol;
    }
  }

  return List::create(
    _["newY"]= newY,
    _["newX"]= newX
  );
}

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

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>, Eigen::ArrayXi> XtX_QR(const Eigen::MatrixXd& X) {


  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR(X);
  const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType Pmat(PQR.colsPermutation());

  int r = PQR.rank();
  int p = X.cols();

  Eigen::MatrixXd R_inv = PQR.matrixQR().topLeftCorner(r, r).triangularView<Eigen::Upper>().solve(Eigen::MatrixXd::Identity(r, r));

  // Get all column indices
  Eigen::ArrayXi Pmat_indices = Pmat.indices();
  // Get the order for the columns you are keeping
  Eigen::ArrayXi Pmat_keep = Pmat_indices.head(r);
  // Get the indices for columns you are discarding
  Eigen::ArrayXi Pmat_toss = Pmat_indices.tail(p - r);

  for(Eigen::Index i=0; i<r; ++i)
  {
    Pmat_keep(i) = Pmat_keep(i) - (Pmat_toss < Pmat_keep(i)).count();
  }

  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P = Eigen::PermutationWrapper<Eigen::ArrayXi>(Pmat_keep);

  Eigen::MatrixXd P_R_inv_P = P * R_inv * P;

  Eigen::MatrixXd XtX_inv = P_R_inv_P * P_R_inv_P.transpose();

  return std::make_tuple(XtX_inv, R_inv, PQR, Pmat_toss);
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
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR;
    std::tie(XtX_inv, R_inv, PQR, std::ignore) = XtX_QR(X);

    r = PQR.rank();

    Eigen::MatrixXd effects(PQR.householderQ().adjoint() * y);

    // Rcout << "effects:" << std::endl;
    // Rcout << effects << std::endl;

    beta_out.topRows(r) = R_inv * effects.topRows(r);
    // Rcout << "beta_out:" << std::endl;
    // Rcout << beta_out << std::endl;
    beta_out = PQR.colsPermutation() * beta_out;
  }

  return List::create(
    _["beta_hat"]= beta_out,
    _["XtX_inv"]= XtX_inv
  );
}

// [[Rcpp::export]]
List lm_variance(Eigen::Map<Eigen::MatrixXd>& X,
                 const Eigen::Map<Eigen::MatrixXd>& XtX_inv,
                 const Eigen::Map<Eigen::MatrixXd>& ei,
                 const Rcpp::Nullable<Rcpp::IntegerVector> & cluster,
                 const int& J,
                 const bool& ci,
                 const String type,
                 const std::vector<bool> & which_covs,
                 const int& fe_rank) {

  const int n(X.rows()), r(XtX_inv.cols()), ny(ei.cols());
  //Rcout << "fe_rank:" << fe_rank << std::endl;
  int r_fe = r + fe_rank;
  const bool clustered = ((type == "stata") || (type == "CR0"));
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
    s2 = AtA(ei)/((double)n - (double)r_fe);
    Vcov_hat = Kr(s2, XtX_inv);

    dof.fill(n - r_fe);

  } else {
    // Robust
    Eigen::MatrixXd temp_omega = ei.array().pow(2);

    s2 = temp_omega.colwise().sum()/((double)n - (double)r_fe);

    Eigen::MatrixXd bread(npars, npars);
    Eigen::MatrixXd half_meat(sandwich_size, npars);
    if (ny == 1) {
      bread = XtX_inv;
    } else {
      bread = Kr(Eigen::MatrixXd::Identity(ny, ny), XtX_inv);
    }
    // Rcout << "temp_omega1: " << std::endl << temp_omega << std::endl;

    if ( !clustered ) {
      // Robust, no clusters

      if ((type == "HC2") || (type == "HC3")) {

        Eigen::MatrixXd meatXtX_inv;
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR;
        Eigen::ArrayXi dropped;
        if (X.cols() > r) {
          std::tie(meatXtX_inv, std::ignore, PQR, dropped) = XtX_QR(X);
          for (Eigen::Index i=0; i<dropped.size(); i++) {
            if (dropped(i) < X.cols())
            X.block(0, dropped(i), n, X.cols() - dropped(i) - 1) = X.rightCols(X.cols() - dropped(i) - 1);
          }

          r_fe -= dropped.size();

        } else {
          meatXtX_inv = XtX_inv;
        }

        // Rcout << "meatXtX_inv:" << std::endl << meatXtX_inv << std::endl;

        for (int i = 0; i < n; i++) {
          Eigen::VectorXd Xi = X.leftCols(r_fe).row(i);
          if (type == "HC2") {
            temp_omega.row(i) = temp_omega.row(i) / (1.0 - (Xi.transpose() * meatXtX_inv * Xi));
          } else if (type == "HC3") {
            temp_omega.row(i) = temp_omega.row(i) / (std::pow(1.0 - Xi.transpose() * meatXtX_inv * Xi, 2));
          }
        }

      }
      // Rcout << "temp_omega2: " << std::endl << temp_omega << std::endl;
      for (int m = 0; m < ny; m++) {
        half_meat.block(0, r*m, n, r) = X.leftCols(r).array().colwise() * temp_omega.col(m).array().sqrt();
      }

      dof.fill(n - r_fe);

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
      (double)n / ((double)n - (double)r_fe);

    // } else if (type == "HC2") {
    //
    //
    //   Vcov_hat =
    //     Vcov_hat *
    //     (double)n / ((double)n - (double)fe_rank);

  } else if (type == "stata") {

    // Rcout << "correction: " << (((double)J * (n - 1)) / (((double)J - 1) * (n - r))) << std::endl;
    Vcov_hat =
      Vcov_hat *
      (((double)J * (n - 1)) / (((double)J - 1) * (n - r_fe)));
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
                     const std::vector<bool> & which_covs,
                     const int& fe_rank) {


  const int n(X.rows()), r(XtX_inv.cols()), ny(ei.cols());
  const int r_fe = r + fe_rank;
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

  Eigen::MatrixXd meatXtX_inv;
  if (X.cols() > r) {
    std::tie(meatXtX_inv, std::ignore, std::ignore, std::ignore) = XtX_QR(X);
  } else {
    meatXtX_inv = XtX_inv;
  }

  // used for the dof corrction
  Eigen::MatrixXd H1s(r_fe, r_fe*J);
  Eigen::MatrixXd H2s(r_fe, r_fe*J);
  Eigen::MatrixXd H3s(r_fe, r_fe*J);
  Eigen::MatrixXd P_diags(r_fe, J);

  Eigen::MatrixXd M_U_ct = meatXtX_inv.llt().matrixL();
  Eigen::MatrixXd MUWTWUM = meatXtX_inv * X.transpose() * X * meatXtX_inv;
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

      // H is not symmetric if weighted CR2
      Eigen::MatrixXd H =
        Xoriginal.block(start_pos, 0, len, r_fe) *
        meatXtX_inv *
        X.block(start_pos, 0, len, r_fe).transpose();

      // Rcout << "H: " << H << std::endl;

      // Code from clubSandwich
      // uwTwu <- Map(function(uw, th) uw %*% th %*% t(uw),
      //             uw = UW_list, th = Theta_list)
      // MUWTWUM <- M_U %*% Reduce("+", uwTwu) %*% M_U

      //(thet - h %*% thet - thet %*% t(h) + u %*% MUWTWUM %*% t(u))

      // A' W R in clubSand notation
      // Rcout << "At_WX: " << (Eigen::MatrixXd::Identity(len, len) - H) - H.transpose() + Xoriginal.block(start_pos, 0, len, r) * MUWTWUM * Xoriginal.block(start_pos, 0, len, r).transpose() << std::endl;
      // Rcout << "MUWTWUM: " << MUWTWUM << std::endl;

      // If no FEs
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> At_WX(
          (Eigen::MatrixXd::Identity(len, len) - H) - H.transpose() +
            Xoriginal.block(start_pos, 0, len, r_fe) *
            MUWTWUM *
            Xoriginal.block(start_pos, 0, len, r_fe).transpose()
      );

      // Rcout << "At_WX RETRIEVED" << std::endl;

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
        X.block(start_pos, 0, len, r_fe);

      // Rcout << "At_WX INVERTED" << std::endl;

      if (ci) {

        Eigen::MatrixXd ME(r_fe, len);
        if (weight_mean != 1) {
          ME = (meatXtX_inv / weight_mean) * At_WX_inv.transpose();
        } else {
          ME = meatXtX_inv * At_WX_inv.transpose();
        }

        P_diags.col(clust_num) = ME.array().pow(2).rowwise().sum();

        Eigen::MatrixXd MEU = ME * Xoriginal.block(start_pos, 0, len, r_fe);

        int p_pos = clust_num*r_fe;
        // Rcout << "p_pos: " << p_pos << std::endl;
        H1s.block(0, p_pos, r_fe, r_fe) = MEU * M_U_ct;
        H2s.block(0, p_pos, r_fe, r_fe) = ME * X.block(start_pos, 0, len, r_fe) * M_U_ct;
        H3s.block(0, p_pos, r_fe, r_fe) = MEU * Omega_ct;
      }

      // t(cr2_eis) %*% (I - P_ss)^{-1/2} %*% Xj
      // each ro  w is the contribution of the cluster to the meat
      // Below use  t(tutX) %*% tutX to sum contributions across clusters
      // Rcout << "At_WX dim:" << At_WX_inv.rows() << "x" << At_WX_inv.cols() << std::endl;

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
          Kr(Eigen::MatrixXd::Identity(ny, ny), At_WX_inv.leftCols(r));

      } else {
        // Rcout << "clust_num: " << clust_num << std::endl;
        // Rcout << "ei:" << std::endl << ei.block(start_pos, 0, len, 1).transpose() << std::endl;
        half_meat.row(clust_num) =
          ei.block(start_pos, 0, len, 1).transpose() *
          At_WX_inv.leftCols(r);
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

  // Rcout << "DONE W SEs" << std::endl;
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
        // Rcout << "H1t size:" << std::endl << H1t.rows() << "x" << H1t.cols() << std::endl<< std::endl;

        H1t.resize(r_fe, J);
        H2t.resize(r_fe, J);
        H3t.resize(r_fe, J);

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
