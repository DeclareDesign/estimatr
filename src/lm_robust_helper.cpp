// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
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

// stats::lm() detects rank with LINPACK dqrdc2, which compares each column's
// remaining norm against ITS OWN original norm, so rescaling a column cannot
// change the rank it reports. Eigen's setThreshold() compares every pivot
// against the LARGEST pivot in the matrix, which is not scale invariant: one
// column in large units pushes the others below the threshold and they are
// dropped as collinear on a design that is full rank. Normalizing the columns
// before the QR makes the criterion per column, which is what dqrdc2 does and
// what the 1e-7 threshold was meant to reproduce. The caller undoes the
// scaling on whatever it reads off the factorization.
Eigen::VectorXd columnScales(const Eigen::Ref<const Eigen::MatrixXd>& X) {
  Eigen::VectorXd scales(X.cols());
  for (Eigen::Index j = 0; j < X.cols(); ++j) {
    const double norm_j = X.col(j).norm();
    scales(j) = (norm_j > 0.0) ? 1.0 / norm_j : 1.0;
  }
  return scales;
}

// The scales of the kept columns, in ascending original-column order, which is
// the order R_inv is permuted back into.
Eigen::VectorXd keptScales(const Eigen::VectorXd& scales,
                           const Eigen::ArrayXi& Pmat_toss,
                           const int p,
                           const int r) {
  Eigen::ArrayXi is_tossed = Eigen::ArrayXi::Zero(p);
  for (Eigen::Index i = 0; i < Pmat_toss.size(); ++i) {
    is_tossed(Pmat_toss(i)) = 1;
  }
  Eigen::VectorXd kept(r);
  Eigen::Index k = 0;
  for (int j = 0; j < p; ++j) {
    if (!is_tossed(j)) kept(k++) = scales(j);
  }
  return kept;
}

// Gets padded UtU matrix (where U = cbind(X, FE_dummies))
Eigen::MatrixXd getMeatXtX(Eigen::Map<Eigen::MatrixXd>& X,
                           const Eigen::MatrixXd& XtX_inv) {
  // Read off X before the compaction below rewrites it.
  const Eigen::VectorXd scales = columnScales(X);
  const Eigen::MatrixXd X_scaled = X * scales.asDiagonal();
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR(X_scaled);
  // The same criterion lm_solver() uses, normalization included, and for the
  // same reason: Eigen's default is tight enough that an exactly collinear
  // column can survive as a pivot of order 1e-14. The two must agree, or the
  // meat is read off a rank the coefficients were not fitted at.
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

  // R_inv came off the normalized design, so it inverts D * XtX * D rather
  // than XtX, where D is diagonal in the column scales.
  const Eigen::VectorXd kept = keptScales(scales, Pmat_toss, p, r);
  Eigen::MatrixXd meatXtX_inv =
    kept.asDiagonal() * (R_inv * R_inv.transpose()) * kept.asDiagonal();

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
    // Normalized for the reason the QR below is, and with a second payoff
    // here. With unit-norm columns the Gram matrix has a unit diagonal, so
    // each L_ii is the norm of that column's residual after the earlier
    // columns are projected out, measured against its own original norm.
    // That is exactly dqrdc2's rank test, and it is what info() does not
    // give: Eigen's LLT reports success on a numerically singular Gram
    // matrix, so testing info() alone let try_cholesky = TRUE return a
    // coefficient for every column of a rank-deficient design, splitting
    // arbitrarily between collinear ones, where the QR path and lm() return
    // NA. A design that is merely ill conditioned falls back to the QR and
    // pays its cost, which is the safe direction to be wrong in.
    const Eigen::VectorXd scales = columnScales(X);
    const Eigen::MatrixXd X_scaled = X * scales.asDiagonal();
    const Eigen::LLT<Eigen::MatrixXd> llt(X_scaled.transpose() * X_scaled);

    if (llt.info() == Eigen::NumericalIssue ||
        llt.matrixLLT().diagonal().minCoeff() < 1e-7) {
      do_qr = true;
    } else {
      beta_out = scales.asDiagonal() * llt.solve(X_scaled.adjoint() * y);
      R_inv = llt.matrixL().solve(Eigen::MatrixXd::Identity(p, p));
      XtX_inv =
        scales.asDiagonal() * (R_inv.transpose() * R_inv) * scales.asDiagonal();
    }
  }

  if (do_qr) {
    const Eigen::VectorXd scales = columnScales(X);
    const Eigen::MatrixXd X_scaled = X * scales.asDiagonal();
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR(X_scaled);
    // Eigen's default rank threshold is about epsilon * ncol relative to the
    // largest pivot, which is tight enough that an exactly collinear column
    // can survive as a pivot of order 1e-14 and produce coefficients of order
    // 1e11 instead of NA (estimatr #351, #395).  stats::lm() uses LINPACK
    // dqrdc2 with tol = 1e-7; matching that takes the normalization above as
    // well as the threshold, since dqrdc2's test is per column.
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

    // The fit is of the normalized design, so each coefficient carries its own
    // column's scale. Applied here, in pivot order, rather than to the whole of
    // beta_out, which would put the dropped columns' NA through an arithmetic
    // operation that need not preserve the payload.
    Eigen::MatrixXd beta_scaled = R_inv * effects.topRows(r);
    for (Eigen::Index i = 0; i < r; ++i) {
      beta_scaled.row(i) *= scales(Pmat_indices(i));
    }
    beta_out.topRows(r) = beta_scaled;
    beta_out = PQR.colsPermutation() * beta_out;

    R_inv = P * R_inv * P;

    const Eigen::VectorXd kept = keptScales(scales, Pmat_toss, p, r);
    XtX_inv = kept.asDiagonal() * (R_inv * R_inv.transpose()) * kept.asDiagonal();

  }

  return List::create(
    _["beta_hat"]= beta_out,
    _["XtX_inv"]= XtX_inv
  );
}

// Satterthwaite degrees of freedom for one row of the CR2 components: a single
// coefficient, or a linear combination of coefficients. Each row of H1s, H2s
// and H3s is linear in the corresponding row of the per-cluster matrix the
// components are built from, and each column of P_diags is that row's squared
// norm, so a combination's components are the combination applied to that
// matrix before squaring and the formula is the same. Avoids the O(J^2) P array
// by computing its trace and Frobenius norm from meat_cols-by-J matrices:
// O(meat_cols^2 * J) against O(J^2).
static double cr2_satterthwaite(const Eigen::MatrixXd& H1s,
                                const Eigen::MatrixXd& H2s,
                                const Eigen::MatrixXd& H3s,
                                const Eigen::MatrixXd& P_diags,
                                const int j,
                                const int meat_cols,
                                const int J) {
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

  return (sq_norm_P > 0.0) ? trace_P * trace_P / sq_norm_P : 0.0;
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
                 const Rcpp::Nullable<Rcpp::NumericVector> & fe_leverage,
                 const int& n_eff,
                 const Rcpp::Nullable<Rcpp::NumericMatrix> & hypotheses = R_NilValue) {

  const int n(X.rows()), r(XtX_inv.cols()), ny(ei.cols());
  // `n` sizes the loops and the matrices; `n_use` counts observations for the
  // degrees of freedom and the HC1/stata scale factors. They differ when rows
  // carry zero weight: such a row contributes nothing to the fit and is not an
  // observation, which is how `lm()` counts it too. -1 means "no distinction".
  const int n_use = (n_eff > 0) ? n_eff : n;
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
  int n_leverage_near_one = 0;

  // Linear combinations of coefficients whose CR2 Satterthwaite degrees of
  // freedom lh_robust() needs. A combination has its own, which is neither any
  // one coefficient's nor bounded by them.
  const bool has_hypotheses = hypotheses.isNotNull() && cr2 && ci;
  int n_hypotheses = 0;
  Eigen::MatrixXd C_hyp, H1c, H2c, H3c, P_hyp;
  if (hypotheses.isNotNull()) {
    n_hypotheses = Rcpp::NumericMatrix(hypotheses).nrow();
  }
  Eigen::VectorXd hypothesis_dof = Eigen::VectorXd::Constant(n_hypotheses, -99.0);
  if (!has_hypotheses) n_hypotheses = 0;

  if (se_type == "classical") {
    Eigen::MatrixXd s2 = AtA(ei)/((double)n_use - (double)r_fe);
    Vcov_hat = Kr(s2, XtX_inv);
    res_var = s2.diagonal();

  } else {
    Eigen::MatrixXd temp_omega = ei.array().pow(2);

    res_var = temp_omega.colwise().sum()/((double)n_use - (double)r_fe);

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
          Rcpp::NumericVector fe_lev = Rcpp::as<Rcpp::NumericVector>(fe_leverage);
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
        //
        // The clamp and the count are deliberately different tests. The clamp
        // acts on the observations whose contribution has to be discarded. The
        // count decides whether R warns, and a strict `denom < 0` there would
        // put the warning at the mercy of one ulp: on an exactly saturated
        // design the solver used here returns a hat value of 1 + 2.2e-16 while
        // `qr()` and `stats::hatvalues()` return exactly 1, and the reported
        // standard error is identical in both cases. So the count uses a
        // tolerance, `sandwich::meatHC`'s `h > 1 - sqrt(eps)` on the same
        // quantity. A hat value is dimensionless and bounded by 1, so the
        // tolerance carries across packages in a way a tolerance on a column
        // norm or a condition number would not. It also reaches an observation
        // sitting just below 1, which is not dropped but whose contribution the
        // small divisor inflates by about 1e8; the warning covers both.
        const double lev_tol = 1.0 - std::sqrt(std::numeric_limits<double>::epsilon());
        n_leverage_near_one = (hii.array() > lev_tol).count();
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

        // One row per hypothesis over the kept coefficients, which are the
        // first columns of the design; any fixed-effect dummy columns after
        // them carry zero weight.
        if (has_hypotheses) {
          Rcpp::NumericMatrix hm(hypotheses);
          n_hypotheses = hm.nrow();
          C_hyp = Eigen::MatrixXd::Zero(n_hypotheses, meat_cols);
          for (int h = 0; h < n_hypotheses; ++h) {
            for (int c = 0; c < hm.ncol(); ++c) C_hyp(h, c) = hm(h, c);
          }
          H1c.resize(n_hypotheses, meat_cols*J);
          H2c.resize(n_hypotheses, meat_cols*J);
          H3c.resize(n_hypotheses, meat_cols*J);
          P_hyp.resize(n_hypotheses, J);
        }
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

              if (has_hypotheses) {
                const Eigen::MatrixXd MEc = C_hyp * ME;
                P_hyp.col(clust_num) = MEc.array().pow(2).rowwise().sum();
                const Eigen::MatrixXd MEUc = MEc * Xoriginal.block(start_pos, 0, len, meat_cols);
                H1c.block(0, p_pos, n_hypotheses, meat_cols) = MEUc * M_U_ct;
                H2c.block(0, p_pos, n_hypotheses, meat_cols) = MEc * X.block(start_pos, 0, len, meat_cols) * M_U_ct;
                H3c.block(0, p_pos, n_hypotheses, meat_cols) = MEUc * Omega_ct;
              }
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
      (double)n_use / ((double)n_use - (double)r_fe);

  } else if (se_type == "stata") {

    Vcov_hat =
      Vcov_hat *
      (((double)J * (n_use - 1)) / (((double)J - 1) * (n_use - r_fe)));
  }

  if (ci) {
    if ( !clustered ) {
      dof.fill(n_use - r_fe);
    } else if (!cr2) {
      dof.fill(J - 1);
    } else {
      for (int j = 0; j < r; j++) {
        if (which_covs[j]) {
          const double dof_j = cr2_satterthwaite(H1s, H2s, H3s, P_diags, j, meat_cols, J);
          for (int outcome_ix = 0; outcome_ix < ny; outcome_ix++) {
            dof(j + outcome_ix * r) = dof_j;
          }
        }
      }
      for (int h = 0; h < n_hypotheses; h++) {
        hypothesis_dof(h) = cr2_satterthwaite(H1c, H2c, H3c, P_hyp, h, meat_cols, J);
      }
    }
  }

  return List::create(_["Vcov_hat"]= Vcov_hat,
                      _["dof"]= dof,
                      _["res_var"]= res_var,
                      _["n_leverage_near_one"]= n_leverage_near_one,
                      _["hypothesis_dof"]= hypothesis_dof);
}

// ---------------------------------------------------------------------------
// demean_cpp: alternating-projections FE demeaning in C++
//
// mat        - N × P matrix to demean (modified in place, returned)
// fe_codes   - list of 1-indexed integer group vectors, one per FE variable
// weights    - length-N weight vector (pass numeric(0) for unweighted)
// eps        - convergence threshold, relative to max|mat|
// max_iter   - maximum number of full sweeps over all FE variables
//
// For one-way FE the algorithm converges in exactly 1 iteration.
// For multi-way FE it cycles through the FE variables until the maximum
// absolute change across all cells is at most eps * max|mat|.
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
  // The magnitude before any sweep, for the convergence floor below.
  const double orig_max = (n > 0 && p > 0) ? mat.cwiseAbs().maxCoeff() : 0.0;
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

    // Converged when the last sweep moved nothing by more than eps relative to
    // the matrix's own magnitude. The test was `eps * (1 + max_abs)`, and the 1
    // made it absolute whenever the matrix is small. The outcome is demeaned on
    // its own, so an outcome in units of 1e-6 stopped after half the sweeps a
    // unit-scale one takes, and on a two-way design whose coefficients were
    // right its residuals came back wrong by 4.4e-4, its fitted values by
    // 3.2e-3, and its standard errors by 7.1e-7. The floor, a millionth of the
    // magnitude before demeaning, is for a column the fixed effects absorb
    // entirely: it demeans to rounding error and would never satisfy a purely
    // relative test. `<=` lets an all-zero matrix converge at once.
    double scale = std::max(max_abs, 1e-6 * orig_max);
    if (max_delta <= eps * scale) { converged = true; break; }
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
