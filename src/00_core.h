#ifndef ADAPTIVEFTS_CORE_H
#define ADAPTIVEFTS_CORE_H

// =============================================================================
// Shared computational core of the C++ estimators (internal, header only).
//
// Every estimator used to convert the DataFrame to a matrix, look the rows of
// each curve up with arma::find() inside its loops, and evaluate the kernel on
// whole vectors through a std::function that allocates at each call. The
// pieces below replace that once for all:
//   - Curves: the data split by curve once (ids sorted as arma::unique, rows of
//     a curve in their original order, and their observation points sorted for
//     the compact-support windows);
//   - scalar kernels (same formulas as the vector kernels of 02_smoothing);
//   - window(): the observation points of a curve that can lie in the support
//     [x - b, x + b] of a kernel, by binary search on the sorted points. The
//     window is a superset (small safety margin): callers keep evaluating the
//     original predicates (|(T - x) / b| <= 1, |T - x| <= b) and kernel values,
//     which vanish outside the support, so the 0/1 indicators are unchanged.
// =============================================================================

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace core {

enum class Kern { epanechnikov, biweight, triweight, tricube, triangular, uniform };

inline Kern parse_kernel(const std::string& name) {
  if (name == "epanechnikov") return Kern::epanechnikov;
  if (name == "biweight") return Kern::biweight;
  if (name == "triweight") return Kern::triweight;
  if (name == "tricube") return Kern::tricube;
  if (name == "triangular") return Kern::triangular;
  if (name == "uniform") return Kern::uniform;
  Rcpp::stop("Unsupported kernel name. Choose from: epanechnikov, biweight, triweight, tricube, triangular, uniform.");
  return Kern::epanechnikov;
}

// Scalar kernel, 0 outside [-1, 1] (a NaN argument gives NaN, as the vector
// kernels do).
inline double kval(Kern k, double u) {
  if (std::abs(u) > 1) return 0.0;
  switch (k) {
  case Kern::epanechnikov: return 0.75 * (1 - u * u);
  case Kern::biweight:     return 15.0 / 16.0 * std::pow(1 - u * u, 2);
  case Kern::triweight:    return 35.0 / 32.0 * std::pow(1 - u * u, 3);
  case Kern::tricube:      return 70.0 / 81.0 * std::pow(1 - std::pow(std::abs(u), 3), 3);
  case Kern::triangular:   return 1 - std::abs(u);
  case Kern::uniform:      return 0.5;
  }
  return 0.0;
}

// Sorted copy of a design and the permutation: ts[j] = t[perm[j]] (ties keep
// their original order).
inline void sort_design(const arma::vec& t, std::vector<double>& ts, std::vector<arma::uword>& perm) {
  const arma::uword M = t.n_elem;
  perm.resize(M);
  for (arma::uword j = 0; j < M; ++j) perm[j] = j;
  const double* Tc = t.memptr();
  std::stable_sort(perm.begin(), perm.end(), [Tc](arma::uword a, arma::uword b) { return Tc[a] < Tc[b]; });
  ts.resize(M);
  for (arma::uword j = 0; j < M; ++j) ts[j] = Tc[perm[j]];
}

// The data split by curve once.
struct Curves {
  int N = 0;
  arma::vec ids;                    // sorted unique curve ids
  arma::vec id_col;                 // the id column (for the lambda estimate)
  std::vector<arma::vec> T, Y;      // per curve, rows in their original order
  std::vector<std::vector<double>> Ts;        // per curve, sorted points
  std::vector<std::vector<arma::uword>> perm; // Ts[j] = T[perm[j]]

  // Mean number of points per curve, computed as the estimators always did.
  double lambdahat() const { return arma::mean(arma::hist(id_col, ids)); }
};

inline Curves make_curves(const Rcpp::DataFrame& data) {
  Curves C;
  C.id_col = Rcpp::as<arma::vec>(data["id_curve"]);
  arma::vec tobs = Rcpp::as<arma::vec>(data["tobs"]);
  arma::vec X = Rcpp::as<arma::vec>(data["X"]);
  C.ids = arma::unique(C.id_col);
  C.N = C.ids.n_elem;
  std::map<double, int> pos;
  for (int c = 0; c < C.N; ++c) pos.emplace(C.ids(c), c);
  std::vector<std::vector<arma::uword>> rows(C.N);
  for (arma::uword r = 0; r < C.id_col.n_elem; ++r) rows[pos[C.id_col(r)]].push_back(r);
  C.T.resize(C.N); C.Y.resize(C.N); C.Ts.resize(C.N); C.perm.resize(C.N);
  for (int c = 0; c < C.N; ++c) {
    arma::uvec ri(rows[c]);
    C.T[c] = tobs.elem(ri);
    C.Y[c] = X.elem(ri);
    sort_design(C.T[c], C.Ts[c], C.perm[c]);
  }
  return C;
}

// Positions [lo, hi) in sorted points ts that can lie within b of x.
// A non-finite or non-positive b, or a non-finite x, keeps every point.
inline void window(const std::vector<double>& ts, double x, double b, arma::uword& lo, arma::uword& hi) {
  if (!(b > 0) || !std::isfinite(b) || !std::isfinite(x)) { lo = 0; hi = ts.size(); return; }
  const double margin = b * 1e-9 + 1e-12;
  lo = std::lower_bound(ts.begin(), ts.end(), x - b - margin) - ts.begin();
  hi = std::upper_bound(ts.begin(), ts.end(), x + b + margin) - ts.begin();
}

inline void window(const Curves& C, int c, double x, double b, arma::uword& lo, arma::uword& hi) {
  window(C.Ts[c], x, b, lo, hi);
}

// Nadaraya-Watson estimate at x with bandwidth b on a sorted design, as
// estimate_nw_cpp: sum_j K((T_j - x) / b) Y_j / sum_j K((T_j - x) / b), NaN when
// the sum is 0. Optionally returns the kernel sum.
inline double nw_sorted(const std::vector<double>& ts, const std::vector<arma::uword>& p,
                        const double* Y, double x, double b, Kern k, double* sumK = nullptr) {
  arma::uword lo, hi;
  window(ts, x, b, lo, hi);
  double S = 0, num = 0;
  for (arma::uword j = lo; j < hi; ++j) {
    double w = kval(k, (ts[j] - x) / b);
    S += w;
    num += w * Y[p[j]];
  }
  if (sumK) *sumK = S;
  return num / S;
}

inline double nw_at(const Curves& C, int c, double x, double b, Kern k) {
  return nw_sorted(C.Ts[c], C.perm[c], C.Y[c].memptr(), x, b, k);
}

// One window pass on curve c: the Nadaraya-Watson estimate at x (NaN when the
// kernel sum is 0), the kernel sum, and the indicator |T - x| <= b.
inline double nw_pn_abs(const Curves& C, int c, double x, double b, Kern k, double& sumK, double& pn) {
  arma::uword lo, hi;
  window(C, c, x, b, lo, hi);
  const std::vector<double>& ts = C.Ts[c];
  const std::vector<arma::uword>& p = C.perm[c];
  const double* Y = C.Y[c].memptr();
  double S = 0, num = 0;
  pn = 0;
  for (arma::uword j = lo; j < hi; ++j) {
    double w = kval(k, (ts[j] - x) / b);
    S += w;
    num += w * Y[p[j]];
    if (std::abs(ts[j] - x) <= b) pn = 1.0;
  }
  sumK = S;
  return num / S;
}

// Indicator that curve c has a point with |(T - x) / b| <= 1.
inline double pn_scaled(const Curves& C, int c, double x, double b) {
  arma::uword lo, hi;
  window(C, c, x, b, lo, hi);
  const std::vector<double>& ts = C.Ts[c];
  for (arma::uword j = lo; j < hi; ++j) if (std::abs((ts[j] - x) / b) <= 1) return 1.0;
  return 0.0;
}

// Indicator that curve c has a point with |T - x| <= b.
inline double pn_abs(const Curves& C, int c, double x, double b) {
  arma::uword lo, hi;
  window(C, c, x, b, lo, hi);
  const std::vector<double>& ts = C.Ts[c];
  for (arma::uword j = lo; j < hi; ++j) if (std::abs(ts[j] - x) <= b) return 1.0;
  return 0.0;
}

// Per-curve quantities of the risk functions at (x, b): pi_n (|(T-x)/b| <= 1),
// the max normalised weight, and b_n(2H) = sum |(T-x)/b|^{2H} |W|, with the
// weights W = K / sum K set to 0 when sum K = 0 (the NaN -> 0 convention).
struct RiskStats { double pn, wmax, sbn; };

inline RiskStats risk_stats(const Curves& C, int c, double x, double b, double H, Kern k) {
  arma::uword lo, hi;
  window(C, c, x, b, lo, hi);
  const std::vector<double>& ts = C.Ts[c];
  double S = 0, kmax = 0, s2h = 0, pn = 0;
  for (arma::uword j = lo; j < hi; ++j) {
    double u = (ts[j] - x) / b;
    double w = kval(k, u);
    if (!std::isfinite(w)) w = 0;
    if (std::abs(u) <= 1) pn = 1.0;
    S += w;
    if (w > kmax) kmax = w;
    s2h += std::pow(std::abs(u), 2 * H) * std::abs(w);
  }
  RiskStats r{pn, 0.0, 0.0};
  if (S != 0 && std::isfinite(S)) { r.wmax = kmax / S; r.sbn = s2h / S; }
  return r;
}

// risk_stats() for every bandwidth of `bws` at once. |(T - x)/h|^{2H} =
// |T - x|^{2H} h^{-2H}: the power of each observation is computed once (on the
// widest window) and reused for every bandwidth. Results for bandwidth b go to
// pn[b * stride], wmax[b * stride], sbn[b * stride].
inline void risk_stats_multi(const Curves& C, int c, double x, const double* bws, int B, double H, Kern k,
                             double* pn, double* wmax, double* sbn, std::size_t stride) {
  const std::vector<double>& ts = C.Ts[c];
  double bmax = 0;
  bool ok = std::isfinite(x);
  for (int b = 0; b < B; ++b) { if (!(bws[b] > 0) || !std::isfinite(bws[b])) ok = false; bmax = std::max(bmax, bws[b]); }
  if (!ok) {   // degenerate inputs: the direct computation
    for (int b = 0; b < B; ++b) {
      RiskStats r = risk_stats(C, c, x, bws[b], H, k);
      pn[b * stride] = r.pn; wmax[b * stride] = r.wmax; sbn[b * stride] = r.sbn;
    }
    return;
  }
  arma::uword lo0, hi0;
  window(ts, x, bmax, lo0, hi0);
  double apow_buf[256];
  std::vector<double> apow_vec;
  double* apow = apow_buf;
  if (hi0 - lo0 > 256) { apow_vec.resize(hi0 - lo0); apow = apow_vec.data(); }
  for (arma::uword j = lo0; j < hi0; ++j) apow[j - lo0] = std::pow(std::abs(ts[j] - x), 2 * H);
  for (int b = 0; b < B; ++b) {
    const double bw = bws[b];
    arma::uword lo, hi;
    window(ts, x, bw, lo, hi);
    double S = 0, kmax = 0, s2h = 0, p = 0;
    for (arma::uword j = lo; j < hi; ++j) {
      double u = (ts[j] - x) / bw;
      double w = kval(k, u);
      if (!std::isfinite(w)) w = 0;
      if (std::abs(u) <= 1) p = 1.0;
      S += w;
      if (w > kmax) kmax = w;
      s2h += apow[j - lo0] * std::abs(w);
    }
    pn[b * stride] = p;
    if (S != 0 && std::isfinite(S)) {
      wmax[b * stride] = kmax / S;
      sbn[b * stride] = s2h * std::pow(bw, -2 * H) / S;
    } else {
      wmax[b * stride] = 0; sbn[b * stride] = 0;
    }
  }
}

// First position of each distinct value of v, and for each element its group.
inline std::vector<int> first_index(const arma::vec& v) {
  std::map<double, int> first;
  std::vector<int> out(v.n_elem);
  for (arma::uword i = 0; i < v.n_elem; ++i) {
    auto it = first.find(v(i));
    if (it == first.end()) { first.emplace(v(i), (int) i); out[i] = (int) i; }
    else out[i] = it->second;
  }
  return out;
}

// Distinct values of v in order of first appearance (vals) and, for each
// element, the position of its value in vals.
inline std::vector<arma::uword> distinct(const arma::vec& v, arma::vec& vals) {
  std::map<double, arma::uword> pos;
  std::vector<double> u;
  std::vector<arma::uword> out(v.n_elem);
  for (arma::uword i = 0; i < v.n_elem; ++i) {
    auto it = pos.find(v(i));
    if (it == pos.end()) { pos.emplace(v(i), u.size()); out[i] = u.size(); u.push_back(v(i)); }
    else out[i] = it->second;
  }
  vals = arma::vec(u);
  return out;
}

// For each element, every position holding the same value, in ascending order
// (what arma::find(v == v(i)) returns).
inline std::vector<std::vector<arma::uword>> same_value_groups(const arma::vec& v) {
  std::map<double, std::vector<arma::uword>> g;
  for (arma::uword i = 0; i < v.n_elem; ++i) g[v(i)].push_back(i);
  std::vector<std::vector<arma::uword>> out(v.n_elem);
  for (arma::uword i = 0; i < v.n_elem; ++i) out[i] = g[v(i)];
  return out;
}

inline std::vector<std::vector<arma::uword>> same_pair_groups(const arma::vec& a, const arma::vec& b) {
  std::map<std::pair<double, double>, std::vector<arma::uword>> g;
  for (arma::uword i = 0; i < a.n_elem; ++i) g[std::make_pair(a(i), b(i))].push_back(i);
  std::vector<std::vector<arma::uword>> out(a.n_elem);
  for (arma::uword i = 0; i < a.n_elem; ++i) out[i] = g[std::make_pair(a(i), b(i))];
  return out;
}

} // namespace core

#endif
