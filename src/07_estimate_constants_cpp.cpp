#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "00_internal.h"
using namespace Rcpp;
using namespace arma;

// The constants of the risk functions, on the split curves (core::Curves): each
// curve is smoothed once at the distinct points, the per-point vectors over the
// curves are built directly (instead of being looked up in a long table of
// smoothed values), and the lag products are taken by direct loops.

namespace {

// Mean over the finite entries; an empty set behaves as arma::mean() does.
double mean_finite(const double* x, arma::uword n) {
  double sum = 0;
  arma::uword k = 0;
  for (arma::uword i = 0; i < n; ++i) if (std::isfinite(x[i])) { sum += x[i]; ++k; }
  if (k == 0) return arma::mean(arma::vec());
  return sum / k;
}

// Mean over the finite entries, NaN and `empty` set when there is none (the
// caller raises arma::mean()'s error outside of any parallel region).
double mean_finite_flag(const double* x, arma::uword n, bool& empty) {
  double sum = 0;
  arma::uword k = 0;
  for (arma::uword i = 0; i < n; ++i) if (std::isfinite(x[i])) { sum += x[i]; ++k; }
  if (k == 0) { empty = true; return arma::datum::nan; }
  return sum / k;
}

// Mean over the finite entries, 0 when there is none (the NaN -> 0 convention
// of the empirical autocovariances).
double mean_finite_or_zero(const double* x, arma::uword n) {
  double sum = 0;
  arma::uword k = 0;
  for (arma::uword i = 0; i < n; ++i) if (std::isfinite(x[i])) { sum += x[i]; ++k; }
  return (k == 0) ? 0.0 : sum / k;
}

// Values of the curves at point t(j): curve by curve, and for each curve every
// position holding the same value as t(j) (as arma::find() on the long table).
arma::vec stack_by_curve(const arma::mat& X, const std::vector<arma::uword>& g) {
  arma::vec x(X.n_rows * g.size());
  arma::uword r = 0;
  for (arma::uword c = 0; c < X.n_rows; ++c)
    for (arma::uword jj : g) x(r++) = X(c, jj);
  return x;
}

} // anonymous namespace

arma::mat smooth_matrix(const core::Curves& C, const arma::vec& t, const arma::vec& h, core::Kern kern) {
  arma::mat X(C.N, t.n_elem);
  for (int c = 0; c < C.N; ++c) {
    double hc = h(C.ids(c) - 1);
    for (arma::uword j = 0; j < t.n_elem; ++j) X(c, j) = core::nw_at(C, c, t(j), hc, kern);
  }
  return X;
}

// The two observations closest to t (sort_index of |T - t| in the original
// order) give Z_n = (Y_i - Y_j)^2 / 2. They are found by binary search on the
// sorted points; when the second and third distances tie (the choice then
// depends on sort_index), the original computation is used.
arma::mat sigma_impl(const core::Curves& C, const arma::vec& t) {
  int n = t.size();
  arma::mat mat_res_sig(n, 2);
  for (int j = 0; j < n; ++j) {
    arma::vec zn_vec(C.N);
    for (int c = 0; c < C.N; ++c) {
      const std::vector<double>& ts = C.Ts[c];
      const arma::vec& Yn = C.Y[c];
      const arma::uword M = ts.size();
      bool done = false;
      if (M >= 3) {
        arma::uword pos = std::lower_bound(ts.begin(), ts.end(), t(j)) - ts.begin();
        arma::uword lo = (pos >= 3) ? pos - 3 : 0, hi = std::min<arma::uword>(M, pos + 3);
        double d[3] = {arma::datum::inf, arma::datum::inf, arma::datum::inf};
        arma::uword q[3] = {0, 0, 0};
        for (arma::uword k = lo; k < hi; ++k) {
          double dk = std::abs(ts[k] - t(j));
          if (dk < d[2]) {
            int m = 2;
            while (m > 0 && dk < d[m - 1]) { d[m] = d[m - 1]; q[m] = q[m - 1]; --m; }
            d[m] = dk; q[m] = k;
          }
        }
        if (d[1] < d[2]) {
          double y0 = Yn(C.perm[c][q[0]]), y1 = Yn(C.perm[c][q[1]]);
          zn_vec(c) = 0.5 * (y0 - y1) * (y0 - y1);
          done = true;
        }
      }
      if (!done) {
        arma::vec abs_diff = arma::abs(C.T[c] - t(j));
        arma::uvec idx = arma::sort_index(abs_diff);
        zn_vec(c) = 0.5 * (Yn(idx(0)) - Yn(idx(1))) * (Yn(idx(0)) - Yn(idx(1)));
      }
    }
    double sig = sqrt(arma::mean(zn_vec));
    mat_res_sig.row(j) = {t(j), sig};
  }
  return mat_res_sig;
}

// [[Rcpp::export]]
arma::mat estimate_sigma_cpp(const Rcpp::DataFrame data, const arma::vec t){
  return sigma_impl(core::make_curves(data), t);
}

arma::mat mom_impl(const core::Curves& C, const arma::vec& t, const arma::vec& h,
                   const double mom_order, const double center, core::Kern kern) {
  int n = t.size();
  arma::mat X = smooth_matrix(C, t, h, kern);
  std::vector<std::vector<arma::uword>> groups = core::same_value_groups(t);
  arma::mat mat_mom_res(n, 3);
  for (int idx_t = 0; idx_t < n; ++idx_t) {
    arma::vec x = stack_by_curve(X, groups[idx_t]);

    // Estimate the mean and center if center = true
    if (center) {
      double muhat = arma::mean(x.elem(arma::find_finite(x)));
      x -= muhat;
    }
    arma::vec x_pow = arma::pow(x, mom_order);

    // Remove NaN values and compute mom_order th moment
    arma::uvec finiteIndices = arma::find_finite(x_pow);
    double Ex_mom_order = 0;
    if (finiteIndices.is_empty()) {
      x_pow.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
      Ex_mom_order = arma::mean(x_pow);
    } else {
      Ex_mom_order = arma::mean(x_pow.elem(finiteIndices));
    }
    mat_mom_res.row(idx_t) = {t(idx_t), mom_order, Ex_mom_order};
  }
  // Replace non-finite value by zero
  mat_mom_res.replace(datum::nan, 0);
  mat_mom_res.replace(datum::inf, 0);
  return mat_mom_res;
}

// [[Rcpp::export]]
arma::mat estimate_empirical_mom_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                     const arma::vec h, const double mom_order,
                                     const double center,
                                     const std::string kernel_name = "epanechnikov"){
  core::Kern kern = core::parse_kernel(kernel_name);
  return mom_impl(core::make_curves(data), t, h, mom_order, center, kern);
}

// Empirical autocovariance of the smoothed curves at each t and each lag, as an
// n x nlag matrix (NaN / Inf set to 0, as in the long table). The curves are
// smoothed and the mean is computed once per t, not once per (t, lag).
arma::mat emp_autocov_matrix(const core::Curves& C, const arma::vec& t, const arma::vec& h,
                             const arma::vec& lag, core::Kern kern) {
  int n = t.size();
  int nlag = lag.size();
  arma::mat X = smooth_matrix(C, t, h, kern);
  std::vector<std::vector<arma::uword>> groups = core::same_value_groups(t);
  arma::mat ac(n, nlag);
  std::vector<double> prod;
  for (int idx_t = 0; idx_t < n; ++idx_t) {
    arma::vec x = stack_by_curve(X, groups[idx_t]);
    const arma::uword N = x.n_elem;
    double muhat = mean_finite(x.memptr(), N);
    for (int l = 0; l < nlag; ++l) {
      const arma::uword L = lag(l);
      if (L > N - 1) Rcpp::stop("'lag' must be lower than the number of curves.");
      const arma::uword m = N - L;
      prod.resize(m);
      for (arma::uword k = 0; k < m; ++k) prod[k] = (x(k) - muhat) * (x(k + L) - muhat);
      double v = mean_finite_or_zero(prod.data(), m);
      ac(idx_t, l) = std::isfinite(v) ? v : 0.0;
    }
  }
  return ac;
}

// [[Rcpp::export]]
arma::mat estimate_empirical_autocov_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                         const arma::vec h, const arma::vec lag,
                                         const std::string kernel_name = "epanechnikov"){
  int n = t.size();
  int nlag = lag.size();
  core::Kern kern = core::parse_kernel(kernel_name);
  arma::mat ac = emp_autocov_matrix(core::make_curves(data), t, h, lag, kern);
  arma::mat mat_autocov_res(nlag * n, 3);
  for (int l = 0; l < nlag; ++l)
    for (int idx_t = 0; idx_t < n; ++idx_t)
      mat_autocov_res.row(l * n + idx_t) = {t(idx_t), lag(l), ac(idx_t, l)};
  mat_autocov_res.replace(datum::nan, 0);
  mat_autocov_res.replace(datum::inf, 0);
  return mat_autocov_res;
}

// Empirical autocovariance of X_n(s) X_{n+cross_lag}(t) for each pair and each
// lag. For pair j the products P_k = X_k(s) X_{k+cross_lag}(t) are formed once
// (curve by curve, pairs equal to (s_j, t_j) pooled as the long-table look-up
// did); gamma = mean of the finite P_k, and the autocovariance at each lag is
// the mean of the finite (P_k - gamma)(P_{k+lag} - gamma), 0 if none. Output:
// pairs x lags matrices of gamma and autocovariance (non-finite set to 0).
void xsxt_matrices(const core::Curves& C, const arma::vec& s, const arma::vec& t,
                   const int cross_lag, const arma::vec& lag, const arma::vec& h,
                   core::Kern kern, const bool center, arma::vec& gamma, arma::mat& acov) {
  const int n = s.size();
  const int nlag = lag.size();
  // Smooth at the distinct values of s and of t only
  arma::vec us, ut;
  std::vector<arma::uword> fs = core::distinct(s, us), ft = core::distinct(t, ut);
  arma::mat Xs = smooth_matrix(C, us, h, kern);
  arma::mat Xt = smooth_matrix(C, ut, h, kern);
  std::vector<std::vector<arma::uword>> groups = core::same_pair_groups(s, t);
  gamma.set_size(n);
  acov.set_size(n, nlag);
  // Size checks before the parallel loop (no exception may leave it)
  for (int j = 0; j < n; ++j) {
    const arma::uword N = C.N * groups[j].size();
    if ((arma::uword) cross_lag > N - 1) Rcpp::stop("'cross_lag' must be lower than the number of curves.");
    for (int l = 0; l < nlag; ++l)
      if ((arma::uword) lag(l) > N - cross_lag - 1)
        Rcpp::stop("'lag' must be lower than the number of curves minus 'cross_lag'.");
  }
  bool empty = false;
#pragma omp parallel for schedule(dynamic)
  for (int j = 0; j < n; ++j) {
    bool empty_j = false;
    const std::vector<arma::uword>& g = groups[j];
    const arma::uword N = C.N * g.size();
    std::vector<double> a(N), b(N);
    arma::uword r = 0;
    for (int c = 0; c < C.N; ++c)
      for (arma::uword jj : g) { a[r] = Xs(c, fs[jj]); b[r] = Xt(c, ft[jj]); ++r; }

    // Estimate the mean, center if necessary
    double muhat_s = mean_finite_flag(a.data(), N, empty_j);
    double muhat_t = mean_finite_flag(b.data(), N, empty_j);
    if (center) for (arma::uword k = 0; k < N; ++k) { a[k] -= muhat_s; b[k] -= muhat_t; }

    const arma::uword m = N - cross_lag;
    std::vector<double> P(m), v(m);
    for (arma::uword k = 0; k < m; ++k) P[k] = a[k] * b[k + cross_lag];
    double g0 = mean_finite_flag(P.data(), m, empty_j);
    gamma(j) = std::isfinite(g0) ? g0 : 0.0;

    for (int l = 0; l < nlag; ++l) {
      const arma::uword L = lag(l);
      for (arma::uword k = 0; k < m - L; ++k) v[k] = (P[k] - g0) * (P[k + L] - g0);
      double val = mean_finite_or_zero(v.data(), m - L);
      acov(j, l) = std::isfinite(val) ? val : 0.0;
    }
    if (empty_j) {
#pragma omp critical
      empty = true;
    }
  }
  if (empty) {
    double err = arma::mean(arma::vec()); // same error as the former computation
    (void) err;
  }
}

// Per-pair quantities of the dependence term of the (auto)covariance risks: the
// lag-0 autocovariance of X(s)X(t) and sum |2 autocov| over the non-zero lags,
// counted once per pair equal to (s_j, t_j), as the long-table look-up did.
void xsxt_dependence(const core::Curves& C, const arma::vec& s, const arma::vec& t,
                     const int cross_lag, const arma::vec& lag, const arma::vec& h,
                     core::Kern kern, const bool center, arma::vec& var0, arma::vec& lrsum) {
  arma::vec gamma;
  arma::mat acov;
  xsxt_matrices(C, s, t, cross_lag, lag, h, kern, center, gamma, acov);
  std::vector<std::vector<arma::uword>> groups = core::same_pair_groups(s, t);
  const int n = s.size();
  var0.set_size(n);
  lrsum.set_size(n);
  for (int j = 0; j < n; ++j) {
    double v0 = arma::datum::nan, sum = 0;
    bool found = false;
    for (arma::uword l = 0; l < lag.n_elem; ++l) {
      if (lag(l) == 0) { if (!found) { v0 = acov(j, l); found = true; } }
      else sum += std::abs(2 * acov(j, l));
    }
    if (!found) Rcpp::stop("The lag-0 autocovariance of X(s)X(t) is required.");
    var0(j) = v0;
    lrsum(j) = sum * groups[j].size();
  }
}

// Function to estimate empirical autocovariance for curves X(s)X(t)
// [[Rcpp::export]]
arma::mat estimate_empirical_XsXt_autocov_cpp(const Rcpp::DataFrame data,
                                              const arma::vec s,
                                              const arma::vec t,
                                              const int cross_lag,
                                              const arma::vec lag,
                                              const arma::vec h,
                                              const std::string kernel_name,
                                              const bool center) {
  int n = s.size();
  int nlag = lag.size();
  double clag = cross_lag;
  core::Kern kern = core::parse_kernel(kernel_name);
  arma::vec gamma;
  arma::mat acov;
  xsxt_matrices(core::make_curves(data), s, t, cross_lag, lag, h, kern, center, gamma, acov);

  arma::mat mat_autocov_XsXt_res(nlag * n, 6);
  for (int j = 0; j < n; ++j)
    for (int l = 0; l < nlag; ++l)
      mat_autocov_XsXt_res.row(j * nlag + l) = {s(j), t(j), clag, lag(l), gamma(j), acov(j, l)};

  // Replace non-finite value by zero and sort
  mat_autocov_XsXt_res.replace(datum::nan, 0);
  mat_autocov_XsXt_res.replace(datum::inf, 0);
  mat_autocov_XsXt_res = mat_autocov_XsXt_res.rows(arma::sort_index(mat_autocov_XsXt_res.col(0)));
  return mat_autocov_XsXt_res ;
}

// :::::::::::::: Estimate \mathbb{D}(t, h_t) numerator ::::::::::::
// Estimate \mathbb{D}(t) for a single t
// [[Rcpp::export]]
double estimate_numerator_DD_single_t(const arma::vec& Znvec_raw, const int max_lag) {
  arma::vec Znvec = Znvec_raw;
  Znvec.replace(datum::nan, 0.0);  // Replace NaN with zero

  int N = Znvec.n_elem;
  double D_val_1 = 0.0;
  double D_val_2 = 0.0;

  // Quartic term
  for (int l1 = 1; l1 <= max_lag - 2; ++l1) {
    for (int l2 = l1 + 1; l2 <= max_lag - 1; ++l2) {
      for (int l3 = l2 + 1; l3 <= max_lag; ++l3) {
        int max_k = N - l3;
        if (max_k <= 0) continue;

        double sum_k = 0.0;
        for (int k = 0; k < max_k; ++k) {
          sum_k += Znvec(k) * Znvec(k + l1) * Znvec(k + l2) * Znvec(k + l3);
        }
        D_val_1 += 24.0 * std::abs(sum_k / max_k);
      }
    }
  }

  // Cubic-square term
  for (int l1 = 1; l1 <= max_lag - 1; ++l1) {
    for (int l2 = l1 + 1; l2 <= max_lag; ++l2) {
      int max_k = N - l2;
      if (max_k <= 0) continue;

      double sum_k = 0.0;
      for (int k = 0; k < max_k; ++k) {
        double z1 = Znvec(k);
        double z2 = Znvec(k + l1);
        double z3 = Znvec(k + l2);
        sum_k += z1 * z1 * z2 * z3 + z1 * z2 * z2 * z3 + z1 * z2 * z3 * z3;
      }
      D_val_2 += 6.0 * std::abs(sum_k / max_k);
    }
  }
  double D_val = D_val_1 + D_val_2;
  return D_val;
}

// Numerator of \mathbb{D}(t; h) and P_N(t; h) for each point and bandwidth, as
// n x n_bw matrices. The centred smoothed curves at t do not depend on the
// bandwidth: they are computed once per point, and only pi_n(t; h) per (t, h).
void dd_matrices(const core::Curves& C, const arma::vec& t, const arma::vec& bw_vec,
                 const arma::vec& h, const int max_lag, core::Kern kern, const bool center,
                 arma::mat& num, arma::mat& PN) {
  const int n = t.size();
  const int n_bw = bw_vec.size();
  arma::mat X = smooth_matrix(C, t, h, kern);
  std::vector<std::vector<arma::uword>> groups = core::same_value_groups(t);
  num.set_size(n, n_bw);
  PN.set_size(n, n_bw);
  for (int j = 0; j < n; ++j) {
    const std::vector<arma::uword>& g = groups[j];
    arma::vec Xhat_t = stack_by_curve(X, g);
    Xhat_t.replace(arma::datum::nan, 0.0);
    Xhat_t.replace(arma::datum::inf, 0.0);
    double muhat_t = arma::mean(Xhat_t.elem(arma::find_finite(Xhat_t)));
    if (center) Xhat_t -= muhat_t;
    for (int b = 0; b < n_bw; ++b) {
      arma::vec pn_vec(Xhat_t.n_elem);
      arma::uword r = 0;
      for (int c = 0; c < C.N; ++c) {
        double pn = core::pn_scaled(C, c, t(j), bw_vec(b));
        for (arma::uword q = 0; q < g.size(); ++q) pn_vec(r++) = pn;
      }
      arma::vec Zn_vec_pn = pn_vec % Xhat_t;
      PN(j, b) = arma::accu(pn_vec);
      num(j, b) = estimate_numerator_DD_single_t(Zn_vec_pn, max_lag);
    }
  }
}

// Function to estimate empirical autocovariance for curves X(s)X(t)
// [[Rcpp::export]]
arma::mat estimate_numerator_dependence_term_DD_cpp(const Rcpp::DataFrame data,
                                                    const arma::vec& t,
                                                    const arma::vec& bw_vec,
                                                    const arma::vec& h,
                                                    const int max_lag,
                                                    const std::string kernel_name,
                                                    const bool center) {
  int n = t.size();
  int n_bw = bw_vec.size();
  core::Kern kern = core::parse_kernel(kernel_name);
  arma::mat num, PN;
  dd_matrices(core::make_curves(data), t, bw_vec, h, max_lag, kern, center, num, PN);

  arma::mat mat_DD_res(n * n_bw, 4);
  for (int b = 0; b < n_bw; ++b) {
    for (int j = 0; j < n; ++j) {
      int idx_res = b * n + j;
      mat_DD_res(idx_res, 0) = t(j);
      mat_DD_res(idx_res, 1) = bw_vec(b);
      mat_DD_res(idx_res, 2) = num(j, b);
      mat_DD_res(idx_res, 3) = PN(j, b);
    }
  }
  return mat_DD_res;
}
