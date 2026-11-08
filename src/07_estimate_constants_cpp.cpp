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

// Every finite entry of x is within 1e150 of c (so the products of centred
// entries cannot overflow); the non-finite entries are NaN, not infinite.
bool bounded_or_nan(const double* x, arma::uword n, double c) {
  for (arma::uword i = 0; i < n; ++i) {
    if (std::isfinite(x[i])) { if (!(std::abs(x[i] - c) < 1e150)) return false; }
    else if (!std::isnan(x[i])) return false;
  }
  return true;
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

// Lag sums, four lags per pass over the data (one load of z[k] for four
// independent accumulators); each lag keeps its own sequential order in k.
//   out[l] = sum_{k < m - L_l} z[k] z[k + L_l]           (cnt == nullptr)
//   and cnt[l] = sum_{k < m - L_l} f[k] f[k + L_l]       (indicators f)
template <bool WITH_CNT>
void lag_sums_t(const double* z, const double* f, arma::uword m, const arma::vec& lag,
                double* out, double* cnt) {
  const arma::uword nlag = lag.n_elem;
  arma::uword l = 0;
  for (; l + 4 <= nlag; l += 4) {
    const arma::uword L0 = lag(l), L1 = lag(l + 1), L2 = lag(l + 2), L3 = lag(l + 3);
    const arma::uword Lmax = std::max(std::max(L0, L1), std::max(L2, L3));
    const arma::uword common = m - Lmax;
    double s0 = 0, s1 = 0, s2 = 0, s3 = 0, c0 = 0, c1 = 0, c2 = 0, c3 = 0;
    for (arma::uword k = 0; k < common; ++k) {
      const double zk = z[k];
      s0 += zk * z[k + L0]; s1 += zk * z[k + L1]; s2 += zk * z[k + L2]; s3 += zk * z[k + L3];
      if (WITH_CNT) { const double fk = f[k]; c0 += fk * f[k + L0]; c1 += fk * f[k + L1]; c2 += fk * f[k + L2]; c3 += fk * f[k + L3]; }
    }
    double* sv[4] = {&s0, &s1, &s2, &s3};
    double* cv[4] = {&c0, &c1, &c2, &c3};
    const arma::uword Ls[4] = {L0, L1, L2, L3};
    for (int j = 0; j < 4; ++j) {
      for (arma::uword k = common; k < m - Ls[j]; ++k) {
        *sv[j] += z[k] * z[k + Ls[j]];
        if (WITH_CNT) *cv[j] += f[k] * f[k + Ls[j]];
      }
      out[l + j] = *sv[j];
      if (WITH_CNT) cnt[l + j] = *cv[j];
    }
  }
  for (; l < nlag; ++l) {
    const arma::uword L = lag(l);
    double sum = 0, c = 0;
    for (arma::uword k = 0; k < m - L; ++k) { sum += z[k] * z[k + L]; if (WITH_CNT) c += f[k] * f[k + L]; }
    out[l] = sum;
    if (WITH_CNT) cnt[l] = c;
  }
}

void lag_sums(const double* z, const double* f, arma::uword m, const arma::vec& lag,
              double* out, double* cnt) {
  if (f) lag_sums_t<true>(z, f, m, lag, out, cnt);
  else lag_sums_t<false>(z, f, m, lag, out, cnt);
}

} // anonymous namespace

arma::mat smooth_matrix(const core::Curves& C, const arma::vec& t, const arma::vec& h, core::Kern kern) {
  arma::mat X(C.N, t.n_elem);
#pragma omp parallel for schedule(static)
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
    bool all_finite = std::isfinite(muhat);
    for (arma::uword k = 0; k < N && all_finite; ++k)
      all_finite = std::isfinite(x(k)) && std::abs(x(k) - muhat) < 1e150;
    for (int l = 0; l < nlag; ++l)
      if ((arma::uword) lag(l) > N - 1) Rcpp::stop("'lag' must be lower than the number of curves.");
    if (all_finite) {
      // Every product finite: plain (vectorised) means
      std::vector<double> z(N), sums(nlag);
      for (arma::uword k = 0; k < N; ++k) z[k] = x(k) - muhat;
      lag_sums(z.data(), nullptr, N, lag, sums.data(), nullptr);
      for (int l = 0; l < nlag; ++l) ac(idx_t, l) = sums[l] / (N - (arma::uword) lag(l));
      continue;
    }
    for (int l = 0; l < nlag; ++l) {
      const arma::uword L = lag(l);
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

    // Fast path when every product is finite (and far from overflow): the
    // lag products (P_k - g)(P_{k+L} - g) are all finite, so the mean of the
    // finite ones is their plain mean, vectorised.
    bool all_finite = std::isfinite(g0);
    for (arma::uword k = 0; k < m && all_finite; ++k)
      all_finite = std::isfinite(P[k]) && std::abs(P[k] - g0) < 1e150;
    if (all_finite) {
      for (arma::uword k = 0; k < m; ++k) v[k] = P[k] - g0;
      std::vector<double> sums(nlag);
      lag_sums(v.data(), nullptr, m, lag, sums.data(), nullptr);
      for (int l = 0; l < nlag; ++l) acov(j, l) = sums[l] / (m - (arma::uword) lag(l));
    } else if (std::isfinite(g0) && bounded_or_nan(P.data(), m, g0)) {
      // Some products are NaN (empty smoothing windows): the mean of the finite
      // lag products is sum(z_k z_{k+L}) / sum(f_k f_{k+L}) with z = P - g set to
      // 0 and f = 0 where P is not finite (two vectorised dot products).
      std::vector<double> f(m);
      for (arma::uword k = 0; k < m; ++k) {
        const bool ok = std::isfinite(P[k]);
        v[k] = ok ? P[k] - g0 : 0.0;
        f[k] = ok ? 1.0 : 0.0;
      }
      std::vector<double> sums(nlag), cnts(nlag);
      lag_sums(v.data(), f.data(), m, lag, sums.data(), cnts.data());
      for (int l = 0; l < nlag; ++l) acov(j, l) = (cnts[l] == 0) ? 0.0 : sums[l] / cnts[l];
    } else {
      for (int l = 0; l < nlag; ++l) {
        const arma::uword L = lag(l);
        for (arma::uword k = 0; k < m - L; ++k) v[k] = (P[k] - g0) * (P[k + L] - g0);
        double val = mean_finite_or_zero(v.data(), m - L);
        acov(j, l) = std::isfinite(val) ? val : 0.0;
      }
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
// ---------------------------------------------------------------------------
// X(s)X(t) dependence quantities on a grid of pairs, by matrix products.
//
// For the pair (s, t), P_k = a_k(s) b_{k+l}(t), k < m (l = cross lag, a and b
// the centred smoothed curves), and the lag-L autocovariance is the mean, over
// the k < m - L where P_k and P_{k+L} are finite, of (P_k - g)(P_{k+L} - g).
// With P0 = P set to 0 where not finite and Z the (few) k where it is not:
//   sum (P_k - g)(P_{k+L} - g) = T1 - g (T2 + T3) + g^2 n_f,
//   T1 = sum_k [a0_k a0_{k+L}] [b0_{k+l} b0_{k+l+L}]   (a0, b0: 0 where not finite)
//   T2 = sum_{k < m-L, k+L not in Z} P0_k,   T3 = sum_{k < m-L, k not in Z} P0_{k+L},
//   n_f = #{k < m - L : k, k + L not in Z}.
// T1 is, for fixed L, a product of an s-only and a t-only factor: for every
// (s, t) of the grid it is one matrix product per lag (4 x 4 blocks of dot
// products). T2, T3 and n_f come from prefix sums of P0, corrected on Z.
namespace {

inline void dot_block(const double* const* A, int ni, const double* const* Bm, int nj,
                      arma::uword len, double out[4][4]) {
  if (ni == 4 && nj == 4) {
    const double *a0 = A[0], *a1 = A[1], *a2 = A[2], *a3 = A[3];
    const double *b0 = Bm[0], *b1 = Bm[1], *b2 = Bm[2], *b3 = Bm[3];
    double c00 = 0, c01 = 0, c02 = 0, c03 = 0, c10 = 0, c11 = 0, c12 = 0, c13 = 0;
    double c20 = 0, c21 = 0, c22 = 0, c23 = 0, c30 = 0, c31 = 0, c32 = 0, c33 = 0;
#pragma omp simd reduction(+:c00,c01,c02,c03,c10,c11,c12,c13,c20,c21,c22,c23,c30,c31,c32,c33)
    for (arma::uword k = 0; k < len; ++k) {
      const double x0 = a0[k], x1 = a1[k], x2 = a2[k], x3 = a3[k];
      const double y0 = b0[k], y1 = b1[k], y2 = b2[k], y3 = b3[k];
      c00 += x0 * y0; c01 += x0 * y1; c02 += x0 * y2; c03 += x0 * y3;
      c10 += x1 * y0; c11 += x1 * y1; c12 += x1 * y2; c13 += x1 * y3;
      c20 += x2 * y0; c21 += x2 * y1; c22 += x2 * y2; c23 += x2 * y3;
      c30 += x3 * y0; c31 += x3 * y1; c32 += x3 * y2; c33 += x3 * y3;
    }
    out[0][0] = c00; out[0][1] = c01; out[0][2] = c02; out[0][3] = c03;
    out[1][0] = c10; out[1][1] = c11; out[1][2] = c12; out[1][3] = c13;
    out[2][0] = c20; out[2][1] = c21; out[2][2] = c22; out[2][3] = c23;
    out[3][0] = c30; out[3][1] = c31; out[3][2] = c32; out[3][3] = c33;
    return;
  }
  for (int i = 0; i < ni; ++i)
    for (int j = 0; j < nj; ++j) {
      const double *x = A[i], *y = Bm[j];
      double c = 0;
#pragma omp simd reduction(+:c)
      for (arma::uword k = 0; k < len; ++k) c += x[k] * y[k];
      out[i][j] = c;
    }
}

} // anonymous namespace

// Returns false (and leaves the outputs untouched) when the grid method does
// not apply; the caller then uses the direct computation.
bool xsxt_dependence_grid(const core::Curves& C, const arma::vec& s, const arma::vec& t,
                          const int cl, const arma::vec& lag, const arma::vec& h,
                          core::Kern kern, const bool center, arma::vec& var0, arma::vec& lrsum) {
  const int n = s.n_elem, N = C.N;
  if (cl < 0 || cl > N - 1) return false;
  const arma::uword m = N - cl;
  if (lag.n_elem != m) return false;
  for (arma::uword l = 0; l < m; ++l) if (lag(l) != (double) l) return false;
  // Distinct points, sorted (an upper-triangular set of pairs then leaves the
  // blocks below the diagonal empty)
  arma::vec vs = arma::unique(s), vt = arma::unique(t);
  const int Us = vs.n_elem, Ut = vt.n_elem;
  if ((double) n < 0.25 * Us * Ut || Us < 4 || Ut < 4) return false;
  std::vector<int> pos((std::size_t) Us * Ut, -1);
  std::vector<arma::uword> is(n), it(n);
  for (int k = 0; k < n; ++k) {
    is[k] = std::lower_bound(vs.begin(), vs.end(), s(k)) - vs.begin();
    it[k] = std::lower_bound(vt.begin(), vt.end(), t(k)) - vt.begin();
    std::size_t idx = (std::size_t) is[k] * Ut + it[k];
    if (pos[idx] >= 0) return false;   // duplicated pair: direct path
    pos[idx] = k;
  }

  // Centred values (0 where not finite), curve-contiguous rows; not-finite
  // positions per row
  arma::mat Xs = smooth_matrix(C, vs, h, kern), Xt = smooth_matrix(C, vt, h, kern);
  auto prep = [&](const arma::mat& X, int U, std::vector<double>& v, std::vector<std::vector<int>>& nanpos) {
    v.assign((std::size_t) U * N, 0.0);
    nanpos.assign(U, std::vector<int>());
    for (int u = 0; u < U; ++u) {
      double sum = 0; int cnt = 0;
      for (int c = 0; c < N; ++c) if (std::isfinite(X(c, u))) { sum += X(c, u); ++cnt; }
      if (cnt == 0) return false;      // the direct path raises the error
      const double mu = sum / cnt;
      for (int c = 0; c < N; ++c) {
        const double x = X(c, u);
        if (!std::isfinite(x)) { nanpos[u].push_back(c); continue; }
        const double y = center ? x - mu : x;
        if (!(std::abs(y) < 1e75)) return false;
        v[(std::size_t) u * N + c] = y;
      }
    }
    return true;
  };
  std::vector<double> a, b;
  std::vector<std::vector<int>> nan_s, nan_t;
  if (!prep(Xs, Us, a, nan_s) || !prep(Xt, Ut, b, nan_t)) return false;

  // Z for a pair: the k < m where a_k or b_{k+l} is not finite
  auto pair_Z = [&](int i, int j, std::vector<int>& Z) {
    Z.clear();
    for (int c : nan_s[i]) if (c < (int) m) Z.push_back(c);
    for (int c : nan_t[j]) if (c - cl >= 0) Z.push_back(c - cl);
    std::sort(Z.begin(), Z.end());
    Z.erase(std::unique(Z.begin(), Z.end()), Z.end());
  };
  // Pairs with many missing values: the corrections would cost more than the
  // direct sums, so those pairs are computed directly (below)
  std::vector<int> direct_pairs;
  {
    std::vector<int> Z;
    for (int k = 0; k < n; ++k) {
      pair_Z(is[k], it[k], Z);
      if ((arma::uword) Z.size() >= m) return false;    // no finite product: error in the direct path
      if (2 * Z.size() > m) {
        direct_pairs.push_back(k);
        pos[(std::size_t) is[k] * Ut + it[k]] = -2;      // not in the grid computation
      }
    }
  }

  std::vector<char> block_needed(((Us + 3) / 4) * ((Ut + 3) / 4), 0);
  const int SB = (Us + 3) / 4, TB = (Ut + 3) / 4;
  for (int k = 0; k < n; ++k)
    if (pos[(std::size_t) is[k] * Ut + it[k]] >= 0) block_needed[(is[k] / 4) * TB + it[k] / 4] = 1;

  std::vector<double> V0((std::size_t) Us * Ut, 0.0), LR((std::size_t) Us * Ut, 0.0);
#pragma omp parallel for schedule(dynamic)
  for (int sb = 0; sb < SB; ++sb) {
    const int i0 = sb * 4, ni = std::min(4, Us - i0);
    bool any = false;
    for (int tb = 0; tb < TB; ++tb) any = any || block_needed[sb * TB + tb];
    if (!any) continue;
    // Per pair of the row block (requested pairs only): g, and for every lag L
    // the sums T2(L), T3(L) and the count n_f(L) (prefix sums of P0 corrected
    // on Z, all the lags at once): cT2, cT3, cNF, each m long per pair.
    std::vector<double> cT2((std::size_t) ni * Ut * m), cT3((std::size_t) ni * Ut * m),
      cNF((std::size_t) ni * Ut * m), gv((std::size_t) ni * Ut, 0.0), pr(m + 1), P0(m);
    std::vector<char> req((std::size_t) ni * Ut, 0);
    std::vector<int> Z;
    for (int i = 0; i < ni; ++i)
      for (int j = 0; j < Ut; ++j) {
        const std::size_t q = (std::size_t) i * Ut + j;
        if (pos[(std::size_t) (i0 + i) * Ut + j] < 0) continue;
        req[q] = 1;
        pair_Z(i0 + i, j, Z);
        const double *ai = &a[(std::size_t) (i0 + i) * N], *bj = &b[(std::size_t) j * N + cl];
        pr[0] = 0;
        for (arma::uword k = 0; k < m; ++k) { P0[k] = ai[k] * bj[k]; pr[k + 1] = pr[k] + P0[k]; }
        // g: the mean of the finite products (the zeros at Z add nothing)
        gv[q] = pr[m] / (double) (m - Z.size());
        double *t2 = &cT2[q * m], *t3 = &cT3[q * m], *nf = &cNF[q * m];
        for (arma::uword L = 0; L < m; ++L) {
          t2[L] = pr[m - L];               // sum_{k < m-L} P0_k
          t3[L] = pr[m] - pr[L];           // sum_{k < m-L} P0_{k+L}
          nf[L] = (double) (m - L);
        }
        for (int z : Z) {
          // k = z - L (k + L in Z), L <= z: drop P0_k from T2, and that k
          for (int L = 0; L <= z; ++L) { t2[L] -= P0[z - L]; nf[L] -= 1; }
          // k = z (k in Z), L < m - z: drop P0_{k+L} from T3, and that k
          for (int L = 0; L < (int) m - z; ++L) { t3[L] -= P0[z + L]; nf[L] -= 1; }
        }
        // k and k + L both in Z were dropped twice
        nf[0] += (double) Z.size();
        for (std::size_t u = 0; u < Z.size(); ++u)
          for (std::size_t v = u + 1; v < Z.size(); ++v) nf[Z[v] - Z[u]] += 1;
      }
    std::vector<double> As(4 * m), Bt(4 * m);
    for (arma::uword L = 0; L < m; ++L) {
      const arma::uword len = m - L;
      for (int i = 0; i < ni; ++i) {
        const double* ai = &a[(std::size_t) (i0 + i) * N];
        double* A1 = &As[i * m];
        for (arma::uword k = 0; k < len; ++k) A1[k] = ai[k] * ai[k + L];
      }
      for (int tb = 0; tb < TB; ++tb) {
        if (!block_needed[sb * TB + tb]) continue;
        const int j0 = tb * 4, nj = std::min(4, Ut - j0);
        for (int j = 0; j < nj; ++j) {
          const double* bj = &b[(std::size_t) (j0 + j) * N + cl];
          double* B1 = &Bt[j * m];
          for (arma::uword k = 0; k < len; ++k) B1[k] = bj[k] * bj[k + L];
        }
        double T[4][4];
        const double* Ap[4] = {&As[0], &As[m], &As[2 * m], &As[3 * m]};
        const double* Bp[4] = {&Bt[0], &Bt[m], &Bt[2 * m], &Bt[3 * m]};
        dot_block(Ap, ni, Bp, nj, len, T);
        for (int i = 0; i < ni; ++i)
          for (int j = 0; j < nj; ++j) {
            const std::size_t q = (std::size_t) i * Ut + (j0 + j);
            if (!req[q]) continue;
            const double g = gv[q];
            const double T2 = cT2[q * m + L], T3 = cT3[q * m + L], nf = cNF[q * m + L];
            double acov = (nf == 0) ? 0.0 : (T[i][j] - g * (T2 + T3) + g * g * nf) / nf;
            if (!std::isfinite(acov)) acov = 0.0;
            const std::size_t idx = (std::size_t) (i0 + i) * Ut + (j0 + j);
            if (L == 0) V0[idx] = acov; else LR[idx] += std::abs(2 * acov);
          }
      }
    }
  }

  // The pairs with many missing values: the direct sums (NaN-aware lag sums)
#pragma omp parallel for schedule(dynamic)
  for (int q = 0; q < (int) direct_pairs.size(); ++q) {
    const int k = direct_pairs[q], i = is[k], j = it[k];
    const double *ai = &a[(std::size_t) i * N], *bj = &b[(std::size_t) j * N + cl];
    std::vector<double> f(m, 1.0), z(m);
    for (int c : nan_s[i]) if (c < (int) m) f[c] = 0.0;
    for (int c : nan_t[j]) if (c - cl >= 0) f[c - cl] = 0.0;
    double sum = 0, cnt = 0;
    for (arma::uword kk = 0; kk < m; ++kk) { z[kk] = ai[kk] * bj[kk]; sum += z[kk]; cnt += f[kk]; }
    const double g = sum / cnt;
    for (arma::uword kk = 0; kk < m; ++kk) z[kk] = (f[kk] != 0) ? z[kk] - g : 0.0;
    std::vector<double> sums(m), cnts(m);
    lag_sums(z.data(), f.data(), m, lag, sums.data(), cnts.data());
    double lr = 0, v0 = 0;
    for (arma::uword L = 0; L < m; ++L) {
      double acov = (cnts[L] == 0) ? 0.0 : sums[L] / cnts[L];
      if (!std::isfinite(acov)) acov = 0.0;
      if (L == 0) v0 = acov; else lr += std::abs(2 * acov);
    }
    const std::size_t idx = (std::size_t) i * Ut + j;
    V0[idx] = v0;
    LR[idx] = lr;
  }

  var0.set_size(n);
  lrsum.set_size(n);
  for (int k = 0; k < n; ++k) {
    const std::size_t idx = (std::size_t) is[k] * Ut + it[k];
    var0(k) = V0[idx];
    lrsum(k) = LR[idx];
  }
  return true;
}

void xsxt_dependence(const core::Curves& C, const arma::vec& s, const arma::vec& t,
                     const int cross_lag, const arma::vec& lag, const arma::vec& h,
                     core::Kern kern, const bool center, arma::vec& var0, arma::vec& lrsum) {
  // Dense grids of pairs: matrix products (see xsxt_dependence_grid)
  if (xsxt_dependence_grid(C, s, t, cross_lag, lag, h, kern, center, var0, lrsum)) return;
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
#pragma omp parallel for schedule(dynamic)
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
