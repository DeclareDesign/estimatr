// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
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
  // The same 1e-7 threshold lm_solver() uses, and for the same reason: Eigen's
  // default is tight enough that an exactly collinear column can survive as a
  // pivot of order 1e-14. The two must agree, or the meat is read off a rank
  // the coefficients were not fitted at.
  PQR.setThreshold(1e-7);
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

  // Compacting X by removing the tossed columns one at a time is only correct
  // in descending index order: each left shift moves every column to the right
  // of the removed one, so a later removal at a HIGHER index would then name
  // the wrong column. The QR hands back its permutation in pivot order, which
  // is descending only by accident. With one redundant column there is nothing
  // to order, which is why every rank-deficient-by-one probe agreed and CR2
  // with fixed effects was wrong only when two or more columns went.
  std::sort(Pmat_toss.data(), Pmat_toss.data() + Pmat_toss.size(),
            std::greater<int>());

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

// `fe_leverage` is the per-observation leverage contributed by absorbed
// one-way fixed effects. Under a single FE factor the hat value of the full
// [dummies | X] design splits exactly, as
//     h_ii = h_ii(demeaned X) + w_i / (sum of w over i's group)
// so HC2 and HC3 need only this vector rather than a full dummy hat matrix.
// It is R_NilValue when there are no fixed effects, or more than one factor,
// where the split does not hold.
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
                 const int& fe_rank,
                 const Rcpp::Nullable<Rcpp::NumericVector> & fe_leverage) {

  const int n(X.rows()), r(XtX_inv.cols()), ny(ei.cols());
  // Two different things, which used to share one variable. `r_fe` is the RANK
  // consumed by the fit, which sets the residual degrees of freedom and the
  // HC1/stata scale factors; absorbed fixed effects consume it whether or not
  // the design matrix carries their columns. `meat_cols` is how many columns
  // of X the hat values are actually read off. They coincide only when the
  // dummies have been expanded into X.
  int r_fe = r + fe_rank;
  int meat_cols = r;
  // Rcpp::String comparison is not free, and these were being re-evaluated
  // once per observation and once per cluster inside the loops below.
  const bool cr2 = (se_type == "CR2");
  const bool hc2 = (se_type == "HC2");
  const bool hc3 = (se_type == "HC3");
  const bool clustered = ((se_type == "stata") || (se_type == "CR0") || cr2);
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
  // Reported back so R can warn on the condition itself rather than on a NaN,
  // which is no longer the symptom once the denominator is guarded.
  int n_leverage_above_one = 0;

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
    if (hc2 || hc3 || cr2) {
      if (X.cols() > r) {
        // The dummies were expanded into X, so the QR reveals the true rank
        // and the meat is read off the full design.
        meatXtX_inv = getMeatXtX(X, XtX_inv);
        meat_cols = meatXtX_inv.cols();
        r_fe = meat_cols;
      } else {
        // The meat is the plain r-by-r XtX_inv and X carries only its r
        // demeaned columns, so the hat values are read off those. `r_fe` keeps
        // the absorbed rank: setting it to r here dropped the fixed effects
        // from the residual degrees of freedom, which moved every p-value and
        // confidence interval on a one-way FE fit at the HC2 default.
        meatXtX_inv = XtX_inv;
        meat_cols = r;
      }
    }

    if ( !clustered ) {

      if (hc2 || hc3) {

        // h_ii is the row-wise quadratic form X_i' M X_i. Taking it as a
        // matrix product plus a row-wise dot product replaces n separate
        // matrix-vector products, each of which also heap-allocated a row.
        //
        // In row blocks, not in one product: with `fixed_effects` expanded
        // into dummies meat_cols is the number of FE levels, and one n-by-that
        // temporary would be hundreds of megabytes. The block is sized to hold
        // about 8 MB, so the ordinary case (a handful of covariates) is
        // still one pass, and the wide case stays bounded.
        Eigen::VectorXd hii(n);
        const int block_rows = std::max(1, std::min(n, (int)(1048576 / meat_cols)));
        for (int start = 0; start < n; start += block_rows) {
          const int len = std::min(block_rows, n - start);
          hii.segment(start, len) =
            (X.block(start, 0, len, meat_cols) * meatXtX_inv)
              .cwiseProduct(X.block(start, 0, len, meat_cols))
              .rowwise().sum();
        }

        if (fe_leverage.isNotNull()) {
          Rcpp::NumericVector fe_lev(fe_leverage);
          for (int i = 0; i < n; i++) hii(i) += fe_lev[i];
        }

        Eigen::ArrayXd denom = 1.0 - hii.array();

        // A hat value is a projection diagonal and cannot exceed 1. Where the
        // computed one does, the observation is fitted exactly up to rounding
        // and its contribution is the same 0/0 that leverage of exactly 1
        // resolves to 0 below. Left alone the two estimators fail differently
        // and neither failure is informative: HC2 divides by a small negative
        // number, half_meat then takes the square root of it, and every
        // standard error in the fit is NaN however small the offending term;
        // HC3 squares the denominator, which cancels the sign, so it returns a
        // finite number carrying a spurious positive term and says nothing.
        // Setting the denominator to 0 sends both through the isfinite trap.
        n_leverage_above_one = (denom < 0.0).count();
        denom = (denom <= 0.0).select(0.0, denom);
        if (hc3) denom = denom.square();

        for (int m = 0; m < ny; m++) {
          temp_omega.col(m) = (temp_omega.col(m).array() / denom)
            .unaryExpr([](double v) {return std::isfinite(v)? v : 0.0;});
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

      if (cr2) {
        Xoriginal.resize(n, r);
        if (Xunweighted.isNotNull()) {
          Xoriginal = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Xunweighted);
        } else {
          Xoriginal = X;
        }

        H1s.resize(meat_cols, meat_cols*J);
        H2s.resize(meat_cols, meat_cols*J);
        H3s.resize(meat_cols, meat_cols*J);
        P_diags.resize(meat_cols, J);

        M_U_ct = meatXtX_inv.llt().matrixL();
        MUWTWUM = meatXtX_inv * X.leftCols(meat_cols).transpose() * X.leftCols(meat_cols) * meatXtX_inv;
        Omega_ct = MUWTWUM.llt().matrixL();
      }

      Eigen::Map<Eigen::ArrayXi> clusters = Rcpp::as<Eigen::Map<Eigen::ArrayXi> >(cluster);

      double current_cluster = clusters(0);
      int clust_num = 0;
      int start_pos = 0;
      int len = 1;

      for (int i = 1; i <= n; ++i){

        if ((i == n) || (clusters(i) != current_cluster)) {

          if (cr2) {

            Eigen::MatrixXd H =
              Xoriginal.block(start_pos, 0, len, meat_cols) *
              meatXtX_inv *
              X.block(start_pos, 0, len, meat_cols).transpose();

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> At_WX(
                (Eigen::MatrixXd::Identity(len, len) - H) - H.transpose() +
                  Xoriginal.block(start_pos, 0, len, meat_cols) *
                  MUWTWUM *
                  Xoriginal.block(start_pos, 0, len, meat_cols).transpose()
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
              X.block(start_pos, 0, len, meat_cols);

            if (ci) {

              Eigen::MatrixXd ME(meat_cols, len);
              if (weight_mean != 1) {
                ME = (meatXtX_inv / weight_mean) * At_WX_inv.transpose();
              } else {
                ME = meatXtX_inv * At_WX_inv.transpose();
              }

              P_diags.col(clust_num) = ME.array().pow(2).rowwise().sum();

              Eigen::MatrixXd MEU = ME * Xoriginal.block(start_pos, 0, len, meat_cols);

              int p_pos = clust_num*meat_cols;
              H1s.block(0, p_pos, meat_cols, meat_cols) = MEU * M_U_ct;
              H2s.block(0, p_pos, meat_cols, meat_cols) = ME * X.block(start_pos, 0, len, meat_cols) * M_U_ct;
              H3s.block(0, p_pos, meat_cols, meat_cols) = MEU * Omega_ct;
            }
          }

          if (ny > 1) {

            Eigen::MatrixXd ei_block = ei.block(start_pos, 0, len, ny);
            Eigen::Map<const Eigen::MatrixXd> ei_long(ei_block.data(), 1, len*ny);

            if (cr2) {
              half_meat.block(clust_num, 0, 1, npars) =
                ei_long *
                Kr(Eigen::MatrixXd::Identity(ny, ny), At_WX_inv.leftCols(r));
            } else {
              half_meat.block(clust_num, 0, 1, npars) =
                ei_long *
                Kr(Eigen::MatrixXd::Identity(ny, ny), X.block(start_pos, 0, len, r));
            }

          } else {

            if (cr2) {
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
    } else if (!cr2) {
      dof.fill(J - 1);
    } else {
      // Avoid O(J^2) P_array by computing trace and Frobenius norm from
      // meat_cols×J matrices directly. Complexity: O(r * meat_cols^2 * J)
      // vs O(r * J^2). For typical cases (meat_cols=3, J=100): ~11x faster.
      for (int j = 0; j < r; j++) {
        if (which_covs[j]) {

          Eigen::MatrixXd H1t = H1s.row(j);
          Eigen::MatrixXd H2t = H2s.row(j);
          Eigen::MatrixXd H3t = H3s.row(j);

          H1t.resize(meat_cols, J);  // meat_cols × J
          H2t.resize(meat_cols, J);
          H3t.resize(meat_cols, J);

          Eigen::RowVectorXd p = P_diags.row(j);  // 1 × J

          // meat_cols × meat_cols products — cheap
          Eigen::MatrixXd G3  = H3t * H3t.transpose();  // symmetric
          Eigen::MatrixXd P31 = H3t * H1t.transpose();
          Eigen::MatrixXd P32 = H3t * H2t.transpose();
          Eigen::MatrixXd G21 = H2t * H1t.transpose();
          Eigen::MatrixXd G11 = H1t * H1t.transpose();  // symmetric
          Eigen::MatrixXd G22 = H2t * H2t.transpose();  // symmetric

          // Column-wise dot products — O(meat_cols * J)
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
                      _["res_var"]= res_var,
                      _["n_leverage_above_one"]= n_leverage_above_one);
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
Rcpp::NumericMatrix demean_cpp(Eigen::MatrixXd mat,
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

  // Weight vector. Unweighted fits never read it, so it is not built.
  const bool unweighted = (weights.size() != n);
  Eigen::VectorXd w;
  if (!unweighted) {
    w.resize(n);
    for (int i = 0; i < n; ++i) w(i) = weights[i];
  }

  // Pre-allocate group-sum buffers (reused across iterations and FE variables)
  int max_grp = *std::max_element(n_grp.begin(), n_grp.end());
  Eigen::MatrixXd wx_sum(max_grp, p);

  // The group weight sums depend only on the codes and the weights, so they
  // are the same on every sweep. They used to be rebuilt inside the iteration
  // loop, which is an O(n) pass per FE variable per sweep, thrown away and
  // recomputed identically the next time round.
  std::vector<Eigen::VectorXd> w_sums(n_fe);
  for (int k = 0; k < n_fe; ++k) {
    const std::vector<int>& g = fe[k];
    w_sums[k] = Eigen::VectorXd::Zero(n_grp[k]);
    if (unweighted) {
      for (int i = 0; i < n; ++i) w_sums[k](g[i]) += 1.0;
    } else {
      for (int i = 0; i < n; ++i) w_sums[k](g[i]) += w(i);
    }
  }

  // max|mat| after the final sweep, accumulated by the subtraction pass below
  // rather than by a separate full pass over the matrix.
  double max_abs = 0.0;
  // Reported back so the caller can say the sweeps ran out. Alternating
  // projections converge geometrically at a rate set by how well connected the
  // factors are, so a weakly connected design can still be moving when the cap
  // is reached, and the answer is then simply wrong rather than approximate.
  int iters_used = 0;
  bool converged = (n_fe <= 1);

  for (int iter = 0; iter < max_iter; ++iter) {
    double max_delta = 0.0;
    iters_used = iter + 1;

    for (int k = 0; k < n_fe; ++k) {
      const std::vector<int>& g = fe[k];
      const int ng = n_grp[k];
      const Eigen::VectorXd& w_sum = w_sums[k];

      // Both `mat` and `wx_sum` are column-major, so everything below walks
      // one column at a time. Touching a row at a time instead strides across
      // the whole matrix on every element, and measured an order of magnitude
      // slower here; it also built a heap-allocated row vector per
      // observation, of which there are only `ng` distinct values.
      wx_sum.topRows(ng).setZero();

      const std::size_t wstride = (std::size_t) wx_sum.rows();

      // An unweighted fit multiplied every element by a 1.0 it had just read
      // out of a vector of ones. The results are bit-identical either way;
      // this is the innermost loop of the whole demeaning.
      for (int c = 0; c < p; ++c) {
        const double* mc = mat.data() + (std::size_t) c * n;
        double* wc = wx_sum.data() + (std::size_t) c * wstride;
        if (unweighted) {
          for (int i = 0; i < n; ++i) wc[g[i]] += mc[i];
        } else {
          for (int i = 0; i < n; ++i) wc[g[i]] += w(i) * mc[i];
        }
      }

      // Group sums become group means once per group rather than once per
      // observation. A group with no observations would divide by zero, and is
      // skipped here exactly as the subtraction below never reaches it.
      for (int c = 0; c < p; ++c) {
        double* wc = wx_sum.data() + (std::size_t) c * wstride;
        for (int j = 0; j < ng; ++j) {
          if (w_sum(j) == 0.0) continue;
          wc[j] /= w_sum(j);
          const double v = std::abs(wc[j]);
          if (v > max_delta) max_delta = v;
        }
      }

      // The last FE variable's subtraction touches every cell of `mat`, so the
      // running maximum it leaves behind IS max|mat| for the finished sweep,
      // which is what the convergence test below needs. Computing it here
      // costs nothing; `mat.cwiseAbs().maxCoeff()` was a second full pass.
      max_abs = 0.0;
      for (int c = 0; c < p; ++c) {
        double* mc = mat.data() + (std::size_t) c * n;
        const double* wc = wx_sum.data() + (std::size_t) c * wstride;
        for (int i = 0; i < n; ++i) {
          mc[i] -= wc[g[i]];
          const double v = std::abs(mc[i]);
          if (v > max_abs) max_abs = v;
        }
      }
    }

    // Convergence: change small relative to current scale
    double scale = 1.0 + max_abs;
    if (max_delta < eps * scale) { converged = true; break; }
  }

  Rcpp::NumericMatrix out(Rcpp::wrap(mat));
  out.attr("iterations") = iters_used;
  out.attr("converged") = converged;
  return out;
}

// Weighted cross-tabulation of two integer code vectors into a dense matrix.
//
// `i1` runs 1..n1 and `i2` runs 1..n2; an observation whose code is below 1 in
// either vector is skipped, which is how the reference level of a
// contrast-coded factor drops out. An empty `w` means unweighted, so the cell
// counts are wanted and no unit weight vector has to be materialised.
//
// This replaces a composite index plus `rowsum()`. `rowsum()` hashes all n
// observations to rediscover groups that are already known to be 1..n1 and
// 1..n2, and then labels its result rows with the group values AS STRINGS, so
// the caller converted them back with `as.integer(rownames(.))`. At n =
// 100,000 over 1,000 by 30 groups that pair cost 7.3 ms against 0.04 ms for a
// single pass, and fe_leverage() does one per factor pair.
// [[Rcpp::export]]
Rcpp::NumericMatrix xtab_cpp(const Rcpp::IntegerVector& i1,
                             const Rcpp::IntegerVector& i2,
                             const int n1,
                             const int n2,
                             const Rcpp::NumericVector& w) {
  const R_xlen_t n = i1.size();
  const bool unweighted = (w.size() == 0);
  Rcpp::NumericMatrix out(n1, n2);

  for (R_xlen_t t = 0; t < n; ++t) {
    const int r = i1[t];
    const int c = i2[t];
    if (r < 1 || c < 1 || r > n1 || c > n2) continue;
    out(r - 1, c - 1) += unweighted ? 1.0 : w[t];
  }

  return out;
}
