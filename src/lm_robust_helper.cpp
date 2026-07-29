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

// Gets padded UtU matrix (where U = cbind(X, FE_dummies))
Eigen::MatrixXd getMeatXtX(Eigen::Map<Eigen::MatrixXd>& X,
                           const Eigen::MatrixXd& XtX_inv) {
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR(X);
  const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType Pmat(PQR.colsPermutation());

  int r = PQR.rank();
  int p = X.cols();

  Eigen::MatrixXd R_inv = PQR.matrixQR().topLeftCorner(r, r).triangularView<Eigen::Upper>().solve(Eigen::MatrixXd::Identity(r, r));

  Eigen::ArrayXi Pmat_indices = Pmat.indices();
  Eigen::ArrayXi Pmat_keep = Pmat_indices.head(r);
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

  bool do_qr = !try_cholesky;
  if (try_cholesky) {
    const Eigen::LLT<Eigen::MatrixXd> llt(X.transpose() * X);

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
    // Eigen's default rank threshold is about epsilon * ncol relative to the
    // largest pivot, which is tight enough that an exactly collinear column
    // can survive as a pivot of order 1e-14 and produce coefficients of order
    // 1e11 instead of NA (estimatr #351, #395).  stats::lm() uses LINPACK
    // dqrdc2 with tol = 1e-7; matching that makes rank detection agree with
    // lm() on degenerate designs.
    PQR.setThreshold(1e-7);
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType Pmat(PQR.colsPermutation());

    r = PQR.rank();

    Eigen::MatrixXd R_inv = PQR.matrixQR().topLeftCorner(r, r).triangularView<Eigen::Upper>().solve(Eigen::MatrixXd::Identity(r, r));

    Eigen::ArrayXi Pmat_indices = Pmat.indices();
    Eigen::ArrayXi Pmat_keep = Pmat_indices.head(r);
    Eigen::ArrayXi Pmat_toss = Pmat_indices.tail(p - r);

    for(Eigen::Index i=0; i<r; ++i)
    {
      Pmat_keep(i) = Pmat_keep(i) - (Pmat_toss < Pmat_keep(i)).count();
    }

    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P = Eigen::PermutationWrapper<Eigen::ArrayXi>(Pmat_keep);
    Eigen::MatrixXd effects(PQR.householderQ().adjoint() * y);

    beta_out.topRows(r) = R_inv * effects.topRows(r);
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
  int r_fe = r + fe_rank;
  const bool clustered = ((se_type == "stata") || (se_type == "CR0") || (se_type == "CR2"));
  const int npars = r * ny;
  int sandwich_size = n;
  if (clustered) {
    sandwich_size = J;
  }

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

  if (se_type == "classical") {
    Eigen::MatrixXd s2 = AtA(ei)/((double)n - (double)r_fe);
    Vcov_hat = Kr(s2, XtX_inv);
    res_var = s2.diagonal();

  } else {
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

    if ( !clustered ) {

      if ((se_type == "HC2") || (se_type == "HC3")) {

        Eigen::ArrayXd new_omega(ny);
        for (int i = 0; i < n; i++) {
          Eigen::VectorXd Xi = X.leftCols(r_fe).row(i);

          if (se_type == "HC2") {
            new_omega = temp_omega.row(i) / (1.0 - (Xi.transpose() * meatXtX_inv * Xi));
          } else if (se_type == "HC3") {
            new_omega = temp_omega.row(i) / (std::pow(1.0 - Xi.transpose() * meatXtX_inv * Xi, 2));
          }
          new_omega = new_omega.unaryExpr([](double v) {return std::isfinite(v)? v : 0.0;});
          temp_omega.row(i) = new_omega;
        }
      }

      for (int m = 0; m < ny; m++) {
        if (ny > 1) {
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

      for (int i = 1; i <= n; ++i){

        if ((i == n) || (clusters(i) != current_cluster)) {

          if (se_type == "CR2") {

            Eigen::MatrixXd H =
              Xoriginal.block(start_pos, 0, len, r_fe) *
              meatXtX_inv *
              X.block(start_pos, 0, len, r_fe).transpose();

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
              H1s.block(0, p_pos, r_fe, r_fe) = MEU * M_U_ct;
              H2s.block(0, p_pos, r_fe, r_fe) = ME * X.block(start_pos, 0, len, r_fe) * M_U_ct;
              H3s.block(0, p_pos, r_fe, r_fe) = MEU * Omega_ct;
            }
          }

          if (ny > 1) {

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

    Vcov_hat = bread * (half_meat.transpose() * half_meat) * bread;

  }

  if (se_type == "HC1") {

    Vcov_hat =
      Vcov_hat *
      (double)n / ((double)n - (double)r_fe);

  } else if (se_type == "stata") {

    Vcov_hat =
      Vcov_hat *
      (((double)J * (n - 1)) / (((double)J - 1) * (n - r_fe)));
  }

  if (ci) {
    if ( !clustered ) {
      dof.fill(n - r_fe);
    } else if (se_type != "CR2") {
      dof.fill(J - 1);
    } else {
      // Avoid O(J^2) P_array by computing trace and Frobenius norm from
      // r_fe×J matrices directly. Complexity: O(r * r_fe^2 * J) vs O(r * J^2).
      // For typical cases (r_fe=3, J=100): ~11x faster for this section.
      for (int j = 0; j < r; j++) {
        if (which_covs[j]) {

          Eigen::MatrixXd H1t = H1s.row(j);
          Eigen::MatrixXd H2t = H2s.row(j);
          Eigen::MatrixXd H3t = H3s.row(j);

          H1t.resize(r_fe, J);  // r_fe × J
          H2t.resize(r_fe, J);
          H3t.resize(r_fe, J);

          Eigen::RowVectorXd p = P_diags.row(j);  // 1 × J

          // r_fe × r_fe products — cheap
          Eigen::MatrixXd G3  = H3t * H3t.transpose();  // symmetric
          Eigen::MatrixXd P31 = H3t * H1t.transpose();
          Eigen::MatrixXd P32 = H3t * H2t.transpose();
          Eigen::MatrixXd G21 = H2t * H1t.transpose();
          Eigen::MatrixXd G11 = H1t * H1t.transpose();  // symmetric
          Eigen::MatrixXd G22 = H2t * H2t.transpose();  // symmetric

          // Column-wise dot products — O(r_fe * J)
          Eigen::RowVectorXd col_sq_A3    = H3t.colwise().squaredNorm();
          Eigen::RowVectorXd col_dot_A1A2 = (H1t.cwiseProduct(H2t)).colwise().sum();

          // trace(P_array) without forming J×J matrix
          double trace_P = H3t.squaredNorm()
                         - 2.0 * H1t.cwiseProduct(H2t).sum()
                         + p.sum();

          // ||P_array||_F^2 without forming J×J matrix
          // P = S - U - U^T + D  (S=H3t^T H3t, U=H1t^T H2t, D=diag(p))
          // ||P||^2 = ||S||^2 - 2<S,Q> + 2<S,D> + ||Q||^2 - 2<Q,D> + ||D||^2
          // where Q = U + U^T (symmetric)
          double sq_norm_P =
              G3.squaredNorm()                                     // ||S||^2 (G3 symmetric → ||G3||_F^2 = trace(G3^2) = trace(S^2))
            - 4.0 * P31.cwiseProduct(P32).sum()                   // -2<S,Q> = -4 trace(S U)
            + 2.0 * col_sq_A3.cwiseProduct(p).sum()               // 2<S,D>
            + 2.0 * G21.cwiseProduct(G21.transpose()).sum()        // 2 trace(U^2)  } ||Q||^2
            + 2.0 * G11.cwiseProduct(G22).sum()                   // 2 ||U||_F^2   }
            - 4.0 * col_dot_A1A2.cwiseProduct(p).sum()            // -2<Q,D>
            + p.squaredNorm();                                     // ||D||^2

          double dof_j = (sq_norm_P > 0.0)
                       ? trace_P * trace_P / sq_norm_P
                       : 0.0;
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

// ---------------------------------------------------------------------------
// demean_cpp: alternating-projections FE demeaning in C++
//
// mat        - N × P matrix to demean (modified in place, returned)
// fe_codes   - list of 1-indexed integer group vectors, one per FE variable
// weights    - length-N weight vector (pass numeric(0) for unweighted)
// eps        - convergence threshold (relative to 1 + max|mat|)
// max_iter   - maximum number of full sweeps over all FE variables
//
// For one-way FE the algorithm converges in exactly 1 iteration.
// For multi-way FE it cycles through the FE variables until the maximum
// absolute change across all cells falls below eps * (1 + max|mat|).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Eigen::MatrixXd demean_cpp(Eigen::MatrixXd mat,
                            Rcpp::List       fe_codes_list,
                            Rcpp::NumericVector weights,
                            double eps      = 1e-8,
                            int    max_iter = 100) {
  const int n = mat.rows();
  const int p = mat.cols();
  const int n_fe = fe_codes_list.size();

  // Unpack FE group codes and group counts
  std::vector<std::vector<int>> fe(n_fe);
  std::vector<int> n_grp(n_fe);
  for (int k = 0; k < n_fe; ++k) {
    Rcpp::IntegerVector gv = Rcpp::as<Rcpp::IntegerVector>(fe_codes_list[k]);
    fe[k].resize(n);
    int mx = 0;
    for (int i = 0; i < n; ++i) {
      fe[k][i] = gv[i] - 1;   // 0-indexed
      if (fe[k][i] > mx) mx = fe[k][i];
    }
    n_grp[k] = mx + 1;
  }

  // Weight vector (all-ones when unweighted)
  Eigen::VectorXd w(n);
  if (weights.size() == n) {
    for (int i = 0; i < n; ++i) w(i) = weights[i];
  } else {
    w.setOnes();
  }

  // Pre-allocate group-sum buffers (reused across iterations and FE variables)
  int max_grp = *std::max_element(n_grp.begin(), n_grp.end());
  Eigen::VectorXd w_sum(max_grp);
  Eigen::MatrixXd wx_sum(max_grp, p);

  for (int iter = 0; iter < max_iter; ++iter) {
    double max_delta = 0.0;

    for (int k = 0; k < n_fe; ++k) {
      const std::vector<int>& g = fe[k];
      const int ng = n_grp[k];

      // Accumulate weighted sums
      w_sum.head(ng).setZero();
      wx_sum.topRows(ng).setZero();
      for (int i = 0; i < n; ++i) {
        const int gi = g[i];
        w_sum(gi) += w(i);
        wx_sum.row(gi) += w(i) * mat.row(i);
      }

      // Subtract weighted group means, tracking max change
      for (int i = 0; i < n; ++i) {
        const int gi = g[i];
        Eigen::RowVectorXd delta = wx_sum.row(gi) / w_sum(gi);
        double row_max = delta.cwiseAbs().maxCoeff();
        if (row_max > max_delta) max_delta = row_max;
        mat.row(i) -= delta;
      }
    }

    // Convergence: change small relative to current scale
    double scale = 1.0 + mat.cwiseAbs().maxCoeff();
    if (max_delta < eps * scale) break;
  }

  return mat;
}
