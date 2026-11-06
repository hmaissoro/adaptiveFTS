#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "00_internal.h"
using namespace Rcpp;
using namespace arma;

// C++ cores of the Rubin-Panaretos (2020) estimators (R/04, R/05). The R
// functions build the observation vectors exactly as before and call these for
// the sums over the observations, which only visit the kernel windows.

//' Kernel sums of the Rubin-Panaretos mean estimator (C++ core)
//'
//' For each t, the local-linear mean
//' \eqn{(Q_0 S_2 - Q_1 S_1) / (S_0 S_2 - S_1^2)} from the pooled observations,
//' with \eqn{Q_p = \sum (T - t)^p Y K((T - t)/h) / (h N)} and
//' \eqn{S_p = \sum (T - t)^p K((T - t)/h) / (h N)}.
//'
//' @param Tn,Yn Pooled observation times (sorted) and values.
//' @param t Evaluation points.
//' @param bw Bandwidth.
//' @param N Number of curves.
//' @param kernel_name Kernel name.
//' @return The mean estimates at `t`.
//' @keywords internal
// [[Rcpp::export]]
arma::vec rp_mean_cpp(const arma::vec Tn, const arma::vec Yn, const arma::vec t,
                      const double bw, const double N, const std::string kernel_name) {
  core::Kern kern = core::parse_kernel(kernel_name);
  std::vector<double> ts(Tn.begin(), Tn.end());
  arma::vec out(t.n_elem);
  for (arma::uword i = 0; i < t.n_elem; ++i) {
    arma::uword lo, hi;
    core::window(ts, t(i), bw, lo, hi);
    double Q0 = 0, Q1 = 0, S0 = 0, S1 = 0, S2 = 0;
    for (arma::uword j = lo; j < hi; ++j) {
      double d = ts[j] - t(i);
      double w = (1 / bw) * core::kval(kern, d / bw);
      Q0 += Yn(j) * w; Q1 += d * Yn(j) * w;
      S0 += w; S1 += d * w; S2 += d * d * w;
    }
    Q0 /= N; Q1 /= N; S0 /= N; S1 /= N; S2 /= N;
    out(i) = (Q0 * S2 - Q1 * S1) / (S0 * S2 - S1 * S1);
  }
  return out;
}

//' Cross sums of the Rubin-Panaretos autocovariance estimator (C++ core)
//'
//' For each pair (s, t), the sum over every couple of observations (k, j) of
//' \eqn{w_j w_k ((x_j - t)/h)^p ((x_k - s)/h)^q K((x_j - t)/h) K((x_k - s)/h) / h^2},
//' the couples k = j excluded when `drop_diag` (lag 0). The double sum
//' factorises into the product of the two one-dimensional sums, each over the
//' kernel window only.
//'
//' @param xk,xj Sorted observation times of the s side and of the t side.
//' @param wk,wj Weights aligned with `xk` and `xj`.
//' @param s,t The pairs.
//' @param bw Bandwidth.
//' @param p,q Powers of the t side and of the s side.
//' @param drop_diag Whether to exclude the couples k = j (requires xk = xj).
//' @param kernel_name Kernel name.
//' @return The sums, one per pair.
//' @keywords internal
// [[Rcpp::export]]
arma::vec rp_pq_cpp(const arma::vec xk, const arma::vec xj, const arma::vec wk, const arma::vec wj,
                    const arma::vec s, const arma::vec t, const double bw, const double p,
                    const double q, const bool drop_diag, const std::string kernel_name) {
  core::Kern kern = core::parse_kernel(kernel_name);
  std::vector<double> vk(xk.begin(), xk.end()), vj(xj.begin(), xj.end());
  auto side = [&](const std::vector<double>& x, const arma::vec& w, double c, double pw) {
    arma::uword lo, hi;
    core::window(x, c, bw, lo, hi);
    double sum = 0;
    for (arma::uword i = lo; i < hi; ++i) {
      double u = (x[i] - c) / bw;
      sum += w(i) * std::pow(u, pw) * core::kval(kern, u);
    }
    return sum;
  };
  std::map<double, double> cache_s, cache_t;
  arma::vec out(s.n_elem);
  for (arma::uword m = 0; m < s.n_elem; ++m) {
    auto is = cache_s.find(s(m));
    double Ss = (is != cache_s.end()) ? is->second : (cache_s[s(m)] = side(vk, wk, s(m), q));
    auto it = cache_t.find(t(m));
    double St = (it != cache_t.end()) ? it->second : (cache_t[t(m)] = side(vj, wj, t(m), p));
    double tot = St * Ss;
    if (drop_diag) {
      arma::uword lo, hi;
      core::window(vk, s(m), bw, lo, hi);
      double d = 0;
      for (arma::uword i = lo; i < hi; ++i) {
        double u = (vk[i] - s(m)) / bw, v = (vj[i] - t(m)) / bw;
        d += wj(i) * wk(i) * std::pow(v, p) * core::kval(kern, v) * std::pow(u, q) * core::kval(kern, u);
      }
      tot -= d;
    }
    out(m) = tot / (bw * bw);
  }
  return out;
}

//' Cross sums of the Rubin-Panaretos estimator with recycled mean vectors (C++ core)
//'
//' Exact emulation, couple by couple, of the former data.table computation when
//' the mean vectors do not have the length of the observation vectors (then the
//' tables of observations and of means are indexed or recycled separately, as R
//' does). Used only in that case; it visits every couple of observations.
//'
//' @param xk,xj,Yk,Yj,Mk,Mj Sorted observation times, values and means of the
//'   s side and of the t side.
//' @param s,t The pairs.
//' @param bw Bandwidth.
//' @param p,q Powers of the t side and of the s side.
//' @param drop_diag Whether to exclude the couples k = j.
//' @param kernel_name Kernel name.
//' @return The sums, one per pair.
//' @keywords internal
// [[Rcpp::export]]
arma::vec rp_q_recycled_cpp(const arma::vec xk, const arma::vec xj, const arma::vec Yk, const arma::vec Yj,
                            const arma::vec Mk, const arma::vec Mj, const arma::vec s, const arma::vec t,
                            const double bw, const double p, const double q, const bool drop_diag,
                            const std::string kernel_name) {
  core::Kern kern = core::parse_kernel(kernel_name);
  const arma::uword nk = xk.n_elem, nj = xj.n_elem, mk = Mk.n_elem, mj = Mj.n_elem;
  const arma::uword Ly = nk * nj, Lm = mk * mj;
  arma::vec out(s.n_elem);
  for (arma::uword m = 0; m < s.n_elem; ++m) {
    double sum = 0;
    auto F = [&](arma::uword k, arma::uword j) {
      double u = (xk(k) - s(m)) / bw, v = (xj(j) - t(m)) / bw;
      return std::pow(v, p) * std::pow(u, q) * (1 / (bw * bw)) * core::kval(kern, v) * core::kval(kern, u);
    };
    if (drop_diag) {
      // One index set (the couples k != j) applied to every table; past the
      // end of the table of means, the means are NA.
      for (arma::uword r = 0; r < Ly; ++r) {
        arma::uword k = r / nj, j = r % nj;
        if (k == j) continue;
        double mu_k = NA_REAL, mu_j = NA_REAL;
        if (r < Lm) { mu_k = Mk(r / mj); mu_j = Mj(r % mj); }
        sum += (Yj(j) - mu_j) * (Yk(k) - mu_k) * F(k, j);
      }
    } else {
      // Element-wise operations on vectors of lengths Ly and Lm, recycled.
      const arma::uword L = std::max(Ly, Lm);
      for (arma::uword r = 0; r < L; ++r) {
        arma::uword ry = r % Ly, rm = r % Lm;
        arma::uword k = ry / nj, j = ry % nj;
        double g = (Yj(j) - Mj(rm % mj)) * (Yk(k) - Mk(rm / mj));
        sum += g * F(k, j);
      }
    }
    out(m) = sum;
  }
  return out;
}
