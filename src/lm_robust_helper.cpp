// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppEigen.h>
using namespace Rcpp;


// [[Rcpp::export]]
Eigen::ArrayXXd demeanMat2(const Eigen::MatrixXd& what,
                          const Rcpp::IntegerMatrix& fes,
                          const Rcpp::NumericVector& weights,
                          const int& start_col,
                          const double& eps) {

  int n = what.rows();
  int p = what.cols() - start_col;

  Eigen::ArrayXXd out(n,p);
  Eigen::ArrayXd oldcol(n);

  out =  what.block(0, start_col, n , p);

  int nlevels = max(fes) + 1; // add one because factor codes are 1:N
  Eigen::ArrayXd group_sums(nlevels);
  Eigen::ArrayXXd group_weights(nlevels, fes.cols());

  group_weights.setZero();
  for (int j = 0; j < fes.cols(); ++j) {
    for (int i=0; i<n; i++) {
      int g = fes(i, j);
      group_weights(g,j) = group_weights(g,j) + weights(i);
    }
  }



  for (int k=0; k<p; k++) {

    double delta = 9999;

    do {
      oldcol = out.col(k);

      for (Eigen::Index j = 0; j < fes.cols(); ++j) {

        group_sums.setZero();
        for (int i=0; i<n; i++) {
          int g = fes(i, j);
          group_sums(g) = group_sums(g) + weights(i) * out(i, k);
        }
        group_sums = group_sums / group_weights.col(j);

        // demean *in-place*
        for (int i=0; i<n; i++) {
          int g = fes(i, j);
          out(i, k) = out(i, k) - group_sums(g,0);
        }

      }

      delta = std::sqrt((oldcol - out.col(k)).pow(2).sum());


      // Rcout << delta << std::endl;
      // Rcout << "**********************************" << std::endl;
      // Rcout << group_sums << std::endl;
      // Rcout << "**********************************" << std::endl;
      // Rcout << "old" << oldcol << std::endl;
      // Rcout << "**********************************" << std::endl;
      // Rcout << "outk" << out.col(k) << std::endl;
      // Rcout << "**********************************" << std::endl;
      // return out;

    } while (delta >= eps);

  }

  return out;
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

// Gets padded UtU matrix (where U = cbind(X, FE_dummies))
Eigen::MatrixXd getMeatXtX(Eigen::Map<Eigen::MatrixXd>& X,
                           const Eigen::MatrixXd& XtX_inv) {
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

  R_inv = P * R_inv * P;

  Eigen::MatrixXd meatXtX_inv = R_inv * R_inv.transpose();

  for (Eigen::Index i=0; i<Pmat_toss.size(); i++) {
    if (Pmat_toss(i) < X.cols())
      X.block(0, Pmat_toss(i), X.rows(), X.cols() - Pmat_toss(i) - 1) = X.rightCols(X.cols() - Pmat_toss(i) - 1);
  }

  return meatXtX_inv;
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
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR(X);
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType Pmat(PQR.colsPermutation());

    r = PQR.rank();

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
    Eigen::MatrixXd effects(PQR.householderQ().adjoint() * y);

    // Rcout << "effects:" << std::endl;
    // Rcout << effects << std::endl;

    beta_out.topRows(r) = R_inv * effects.topRows(r);
    // Rcout << "beta_out:" << std::endl;
    // Rcout << beta_out << std::endl;
    beta_out = PQR.colsPermutation() * beta_out;

    R_inv = P * R_inv * P;

    XtX_inv = R_inv * R_inv.transpose();

  }

  return List::create(
    _["beta_hat"]= beta_out,
    _["XtX_inv"]= XtX_inv
  );
}

// [[Rcpp::export]]
List lm_variance(Eigen::Map<Eigen::MatrixXd>& X,
                 const Rcpp::Nullable<Rcpp::NumericMatrix> & Xunweighted,
                 const Eigen::Map<Eigen::MatrixXd>& XtX_inv,
                 const Eigen::Map<Eigen::MatrixXd>& ei,
                 const double weight_mean,
                 const Rcpp::Nullable<Rcpp::IntegerVector> & cluster,
                 const int& J,
                 const bool& ci,
                 const String se_type,
                 const std::vector<bool> & which_covs,
                 const int& fe_rank) {

  const int n(X.rows()), r(XtX_inv.cols()), ny(ei.cols());
  // Rcout << "X:" << std::endl << X << std::endl;
  int r_fe = r + fe_rank;
  const bool clustered = ((se_type == "stata") || (se_type == "CR0") || (se_type == "CR2"));
  const int npars = r * ny;
  int sandwich_size = n;
  if (clustered) {
    sandwich_size = J;
  }

  // For CR2
  Eigen::MatrixXd Xoriginal;
  Eigen::MatrixXd H1s;
  Eigen::MatrixXd H2s;
  Eigen::MatrixXd H3s;
  Eigen::MatrixXd P_diags;
  Eigen::MatrixXd M_U_ct;
  Eigen::MatrixXd MUWTWUM;
  Eigen::MatrixXd Omega_ct;
  Eigen::MatrixXd At_WX_inv;

  Eigen::MatrixXd Vcov_hat;
  Eigen::VectorXd dof = Eigen::VectorXd::Constant(npars, -99.0);
  Eigen::VectorXd res_var = Eigen::VectorXd::Constant(ny, -99.0);

    // Standard error calculations
  if (se_type == "classical") {
    // Classical
    Eigen::MatrixXd s2 = AtA(ei)/((double)n - (double)r_fe);
    Vcov_hat = Kr(s2, XtX_inv);
    res_var = s2.diagonal();

  } else {
    // Robust
    Eigen::MatrixXd temp_omega = ei.array().pow(2);

    res_var = temp_omega.colwise().sum()/((double)n - (double)r_fe);

    Eigen::MatrixXd bread(npars, npars);
    Eigen::MatrixXd half_meat(sandwich_size, npars);
    if (ny == 1) {
      bread = XtX_inv;
    } else {
      bread = Kr(Eigen::MatrixXd::Identity(ny, ny), XtX_inv);
    }

    Eigen::MatrixXd meatXtX_inv;
    if ((se_type == "HC2") || (se_type == "HC3") || (se_type == "CR2")) {
      if (X.cols() > r) {
        meatXtX_inv = getMeatXtX(X, XtX_inv);
        r_fe = meatXtX_inv.cols();
      } else {
        meatXtX_inv = XtX_inv;
      }
    }
    // Rcout << "meatXtX_inv:" << std::endl << meatXtX_inv << std::endl;

    if ( !clustered ) {
      // Rcout << "temp_omega:" << std::endl << temp_omega << std::endl;

      if ((se_type == "HC2") || (se_type == "HC3")) {

        Eigen::ArrayXd new_omega(ny);
        for (int i = 0; i < n; i++) {
          Eigen::VectorXd Xi = X.leftCols(r_fe).row(i);
          // Rcout << i << ":" << Xi << std::endl;

          if (se_type == "HC2") {
            new_omega = temp_omega.row(i) / (1.0 - (Xi.transpose() * meatXtX_inv * Xi));
          } else if (se_type == "HC3") {
            new_omega = temp_omega.row(i) / (std::pow(1.0 - Xi.transpose() * meatXtX_inv * Xi, 2));
          }
          // Perfect fits cause instability, but we can place 0s for those
          // observations and the rest of the estimation works
          new_omega = new_omega.unaryExpr([](double v) {return std::isfinite(v)? v : 0.0;});
          temp_omega.row(i) = new_omega;
        }
      }
      // Rcout << "temp_omega:" << std::endl << temp_omega << std::endl;


      for (int m = 0; m < ny; m++) {
        if (ny > 1) {
          // Preserve signs for off-diagonal vcov blocks in mlm
          half_meat.block(0, r*m, n, r) =  X.leftCols(r).array().colwise() * (ei.col(m).array().sign() * temp_omega.col(m).array().sqrt());
        } else {
          half_meat.block(0, r*m, n, r) =  X.leftCols(r).array().colwise() * temp_omega.col(m).array().sqrt();
        }
      }

    } else {
      // clustered

      if (se_type == "CR2") {
        Xoriginal.resize(n, r);
        if (Xunweighted.isNotNull()) {
          Xoriginal = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Xunweighted);
        } else {
          Xoriginal = X;
        }

        H1s.resize(r_fe, r_fe*J);
        H2s.resize(r_fe, r_fe*J);
        H3s.resize(r_fe, r_fe*J);
        P_diags.resize(r_fe, J);

        M_U_ct = meatXtX_inv.llt().matrixL();
        MUWTWUM = meatXtX_inv * X.leftCols(r_fe).transpose() * X.leftCols(r_fe) * meatXtX_inv;
        Omega_ct = MUWTWUM.llt().matrixL();
      }

      Eigen::Map<Eigen::ArrayXi> clusters = Rcpp::as<Eigen::Map<Eigen::ArrayXi> >(cluster);

      double current_cluster = clusters(0);
      int clust_num = 0;
      int start_pos = 0;
      int len = 1;

      // iterate over unique cluster values
      for (int i = 1; i <= n; ++i){

        if ((i == n) || (clusters(i) != current_cluster)) {

          if (se_type == "CR2") {

            // H is not symmetric if weighted CR2
            Eigen::MatrixXd H =
              Xoriginal.block(start_pos, 0, len, r_fe) *
              meatXtX_inv *
              X.block(start_pos, 0, len, r_fe).transpose();

            // Code from clubSandwich
            // uwTwu <- Map(function(uw, th) uw %*% th %*% t(uw),
            //             uw = UW_list, th = Theta_list)
            // MUWTWUM <- M_U %*% Reduce("+", uwTwu) %*% M_U

            //(thet - h %*% thet - thet %*% t(h) + u %*% MUWTWUM %*% t(u))

            // A' W R in clubSand notation

            // If no FEs
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> At_WX(
                (Eigen::MatrixXd::Identity(len, len) - H) - H.transpose() +
                  Xoriginal.block(start_pos, 0, len, r_fe) *
                  MUWTWUM *
                  Xoriginal.block(start_pos, 0, len, r_fe).transpose()
            );

            Eigen::VectorXd eigvals = At_WX.eigenvalues();
            for (int m = 0; m < eigvals.size(); ++m) {
              if (eigvals(m) > std::pow(10.0, -12.0)) {
                eigvals(m) = 1.0 / std::sqrt(eigvals(m));
              } else {
                eigvals(m) = 0;
              }
            }

            At_WX_inv =
              At_WX.eigenvectors() *
              eigvals.asDiagonal() *
              At_WX.eigenvectors().transpose() *
              X.block(start_pos, 0, len, r_fe);

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
          }

          if (ny > 1) {

            // Stack residuals for this cluster from each model
            // Rcout << "len: " << len << std::endl;
            Eigen::MatrixXd ei_block = ei.block(start_pos, 0, len, ny);
            Eigen::Map<const Eigen::MatrixXd> ei_long(ei_block.data(), 1, len*ny);

            if (se_type == "CR2") {
              half_meat.block(clust_num, 0, 1, npars) =
                ei_long *
                Kr(Eigen::MatrixXd::Identity(ny, ny), At_WX_inv.leftCols(r));
            } else {
              half_meat.block(clust_num, 0, 1, npars) =
                ei_long *
                Kr(Eigen::MatrixXd::Identity(ny, ny), X.block(start_pos, 0, len, r));
            }

          } else {

            if (se_type == "CR2") {
              half_meat.row(clust_num) =
                ei.block(start_pos, 0, len, 1).transpose() *
                At_WX_inv.leftCols(r);
            } else {
              half_meat.row(clust_num) =
                ei.block(start_pos, 0, len, 1).transpose() *
                X.block(start_pos, 0, len, r);
            }

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

    // Rcout << "bread:" << std::endl << bread << std::endl << std::endl;
    // Rcout << "half_meat:" << std::endl << half_meat << std::endl << std::endl;
    // Rcout << "meat:" << std::endl << (half_meat.transpose() * half_meat) << std::endl << std::endl;

    Vcov_hat = bread * (half_meat.transpose() * half_meat) * bread;

  }

  if (se_type == "HC1") {

    Vcov_hat =
      Vcov_hat *
      (double)n / ((double)n - (double)r_fe);

  } else if (se_type == "stata") {

    // Rcout << "correction: " << (((double)J * (n - 1)) / (((double)J - 1) * (n - r))) << std::endl;
    Vcov_hat =
      Vcov_hat *
      (((double)J * (n - 1)) / (((double)J - 1) * (n - r_fe)));
  }

  // Degrees of freedom
  if (ci) {
    if ( !clustered ) {
      dof.fill(n - r_fe); // regular
    } else if (se_type != "CR2") {
      dof.fill(J - 1); // clustered
    } else {
      for (int j = 0; j < r; j++) {
        if (which_covs[j]) {

          Eigen::MatrixXd H1t = H1s.row(j);
          Eigen::MatrixXd H2t = H2s.row(j);
          Eigen::MatrixXd H3t = H3s.row(j);

          H1t.resize(r_fe, J);
          H2t.resize(r_fe, J);
          H3t.resize(r_fe, J);

          Eigen::MatrixXd uf = H1t.transpose() * H2t;
          Eigen::MatrixXd P_row = P_diags.row(j).asDiagonal();
          Eigen::MatrixXd P_array = (H3t.transpose()*H3t - uf - uf.transpose()) + P_row;

          double dof_j = std::pow(P_array.trace(), 2) / P_array.array().pow(2).sum();
          for (int outcome_ix = 0; outcome_ix < ny; outcome_ix++) {
            dof(j + outcome_ix * r) = dof_j;
          }
        }
      }
    }
  }

  return List::create(_["Vcov_hat"]= Vcov_hat,
                      _["dof"]= dof,
                      _["res_var"]= res_var);
}
