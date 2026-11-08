#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <map>
#include <set>
#include <utility>
#include <cstdint>

// [[Rcpp::depends(RcppArmadillo)]]

#include "00_internal.h"
using namespace Rcpp;
using namespace arma;

// Default bandwidth grid of the (auto)covariance risk.
arma::vec default_autocov_bw_grid(const core::Curves& C) {
  double n_curve = C.N;
  double lambdahat = C.lambdahat();
  double b0 = 5 * std::max(std::pow(n_curve * lambdahat, -1 / (2 * 0.1 + 1)),
                           std::pow(n_curve * lambdahat * lambdahat, -1 / (2 * 0.1 + 1)));
  double bK = 0.2;
  return arma::logspace(log10(b0), log10(bK), 15);
}

// Risk of the (auto)covariance estimator (Myriam Vimond's phase B restructuring,
// in C++).
//   - The local regularity, noise level, second moments and dependence
//     numerators depend on a single point: they are computed once on the
//     distinct points of s and t together.
//   - The per-curve quantities pi_n, b_n(2H) and max W depend on (point, h, curve)
//     only: they are computed once, from the kernel window of each curve, and laid
//     out contiguously over the curves.
//   - Every sum over the curves is a dot product of an s-side and a t-side
//     vector (see below), run contiguously over the curves (vectorised).
//   - The point constants can be computed once for several risks on the same
//     points (AutocovPoints, used by the BLUP for C0 and C1).
// With `full`, the risk matrix of estimate_autocov_risk_cpp is filled; with
// `best`, the first minimiser of the finite risks is kept per pair, in the row
// order of that matrix, without materialising it.
AutocovPoints autocov_points(const core::Curves& C, const arma::vec& pts, const arma::vec& bw_grid_to_use,
                             const bool center, core::Kern kern, const LocregCfg& cfg) {
  AutocovPoints P;
  const int nc = C.N, B = bw_grid_to_use.size();
  P.U = arma::unique(pts);
  P.B = B;
  P.center = center;
  P.bw_grid = bw_grid_to_use;
  const int nU = P.U.n_elem;
  P.locU = locreg_impl(C, P.U, true, kern, cfg);
  P.h0 = P.locU(0, 1);
  arma::vec h(nc, fill::value(P.h0)); // the presmoothing bandwidth
  P.sigU = sigma_impl(C, P.U);
  P.momU = mom_impl(C, P.U, h, 2, center, kern);
  dd_matrices(C, P.U, bw_grid_to_use, h, 3, kern, center, P.ddnum, P.ddPN);
  const std::size_t cube = (std::size_t) nU * B * nc;
  P.c_pn.resize(cube); P.c_sbn.resize(cube); P.c_max.resize(cube);
  // Every bandwidth at once per (point, curve) (shared powers |T - x|^{2H})
#pragma omp parallel for schedule(dynamic)
  for (int uc = 0; uc < nU * nc; ++uc) {
    const int u = uc / nc, c = uc % nc;
    const double H = P.locU(u, 4), x = P.U(u);
    const std::size_t off = (std::size_t) u * B * nc + c;
    core::risk_stats_multi(C, c, x, bw_grid_to_use.memptr(), B, H, kern,
                           &P.c_pn[off], &P.c_max[off], &P.c_sbn[off], nc);
  }
  return P;
}

arma::mat autocov_risk_impl(const core::Curves& C, const arma::vec& s, const arma::vec& t,
                            const int lag, const arma::vec& bw_grid_to_use, const bool common_bw,
                            const bool center, core::Kern kern, const LocregCfg& cfg,
                            const bool full, AutocovBest* best, const AutocovPoints* shared) {
  const int n = t.size();
  const int nc = C.N;
  const int B = bw_grid_to_use.size();

  // ---- Point constants, once on the distinct points of s and t (or shared)
  AutocovPoints own;
  if (!shared) own = autocov_points(C, arma::join_cols(s, t), bw_grid_to_use, center, kern, cfg);
  const AutocovPoints& PC = shared ? *shared : own;
  const arma::vec& U = PC.U;
  const arma::mat &locU = PC.locU, &sigU = PC.sigU, &momU = PC.momU, &ddnum = PC.ddnum, &ddPN = PC.ddPN;
  const double h0 = PC.h0;
  arma::vec h(nc, fill::value(h0)); // the presmoothing bandwidth
  const std::vector<double> &c_pn = PC.c_pn, &c_sbn = PC.c_sbn, &c_max = PC.c_max;
  arma::vec XsXtvar_k, sumabslr_k;
  xsxt_dependence(C, s, t, lag, arma::regspace(0, nc - lag - 1), h, kern, center, XsXtvar_k, sumabslr_k);

  auto posU = [&U](double x) { return (int) (std::lower_bound(U.begin(), U.end(), x) - U.begin()); };
  std::vector<int> ks(n), kt(n);
  for (int k = 0; k < n; ++k) { ks[k] = posU(s(k)); kt[k] = posU(t(k)); }

  const int nU = U.n_elem;
  // Every sum of the risk over the curves is a dot product between an s-side
  // vector (pi_n(s), b_n(s), max W_n(s), curve n) and a t-side vector (same
  // quantities for the curve n + lag): pi_n = 0 implies b_n = max W_n = 0, so
  // the inactive curves contribute exact zeros and the products can run over
  // every curve, contiguously (vectorised). The only exception, sum of
  // max W_{n+lag}(t) of the two-bandwidth variance term, does not depend on s
  // and is computed once per (point, bandwidth).
  const int nn = std::max(nc - lag, 0);
  std::vector<double> tmax_sum(nU * B, 0.0);
  for (int ub = 0; ub < nU * B; ++ub) {
    const std::size_t off = (std::size_t) ub * nc;
    for (int i = 0; i < nn; ++i) if (c_pn[off + i + lag] != 0) tmax_sum[ub] += c_max[off + i + lag];
  }

  // Nested supports. pi_n(x; h) = 1{some |T_n - x| / h <= 1} can only switch
  // from 0 to 1 as h grows (IEEE division is monotone), so along the bandwidths
  // sorted increasingly each curve enters the support once and stays. Sums of
  // s-side quantities over the curves active on the t side (P_N, B_1, C_1),
  // and of t-side quantities over the curves active on the s side (B_2), are
  // then cumulative sums along the entry order: one pass over the curves gives
  // every bandwidth, instead of one dot product per bandwidth. The nesting is
  // checked; a point where it fails uses the dense products.
  std::vector<int> ord(B);
  for (int b = 0; b < B; ++b) ord[b] = b;
  std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return bw_grid_to_use(a) < bw_grid_to_use(b); });
  // Curves grouped by entry rank, per point and side: entry_*[u] lists the
  // curves (index n of the s-side curve n, t-side curve n + lag) by increasing
  // rank, start_*[u][p] .. start_*[u][p + 1] those entering at rank p.
  std::vector<std::vector<int>> entry_s(nU), entry_t(nU), start_s(nU), start_t(nU);
  std::vector<char> nested_s(nU, 1), nested_t(nU, 1);
  auto build_entry = [&](int u, int shift, std::vector<int>& entry, std::vector<int>& start, char& nested) {
    std::vector<int> rank(nn, B);
    for (int i = 0; i < nn; ++i) {
      for (int pr = 0; pr < B; ++pr) {
        const double pn = c_pn[((std::size_t) u * B + ord[pr]) * nc + i + shift];
        if (rank[i] == B) { if (pn != 0) rank[i] = pr; }
        else if (pn == 0) nested = 0;
      }
    }
    start.assign(B + 1, 0);
    for (int i = 0; i < nn; ++i) if (rank[i] < B) ++start[rank[i] + 1];
    for (int pr = 0; pr < B; ++pr) start[pr + 1] += start[pr];
    entry.assign(start[B], 0);
    std::vector<int> pos(start.begin(), start.end() - 1);
    for (int i = 0; i < nn; ++i) if (rank[i] < B) entry[pos[rank[i]]++] = i;
  };
  if (!common_bw) {
#pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < nU; ++u) {
      build_entry(u, 0, entry_s[u], start_s[u], nested_s[u]);
      build_entry(u, lag, entry_t[u], start_t[u], nested_t[u]);
    }
  }

  // D(u; h) = num_D / P_N^3 once per (point, bandwidth)
  arma::mat DU(nU, B);
  for (int u = 0; u < nU; ++u)
    for (int b = 0; b < B; ++b) DU(u, b) = ddnum(u, b) / std::pow(ddPN(u, b), 3);

  // ---- Risk
  const int n_rows_res = common_bw ? B * n : B * B * n;
  arma::mat mat_res_risk;
  if (full) mat_res_risk.set_size(n_rows_res, 14);
  if (best) {
    best->bw_s.set_size(n); best->bw_t.set_size(n);
    best->Hs.set_size(n); best->Ls.set_size(n); best->Ht.set_size(n); best->Lt.set_size(n);
  }
  bool any_missing = false;

#pragma omp parallel for schedule(dynamic)
  for (int k = 0; k < n; ++k) {
    const int us = ks[k], ut = kt[k];
    const double Hs = locU(us, 4), Ls = locU(us, 5), Ht = locU(ut, 4), Lt = locU(ut, 5);
    const double mom_s_square = momU(us, 2), mom_t_square = momU(ut, 2);
    const double sig_s_square = std::pow(sigU(us, 1), 2), sig_t_square = std::pow(sigU(ut, 1), 2);
    const double XsXtvar = XsXtvar_k(k), sumabslr = sumabslr_k(k);
    // As the former selection: the position of the first minimum among the
    // finite risks, applied to the list of all the bandwidth combinations.
    double best_risk = 0;
    arma::uword best_rank = 0, n_finite = 0;

    auto finish = [&](double bw_s, double bw_t, int idx_bw_s, int idx_bw_t, double PN_lag,
                      double bias_term_num, double variance_term_num, arma::uword row) {
      double bias_term = 4 * bias_term_num / PN_lag;
      double variance_term = 4 * variance_term_num / (PN_lag * PN_lag);
      double num_Ds = ddnum(us, idx_bw_s), num_Dt = ddnum(ut, idx_bw_t);
      double Ds = DU(us, idx_bw_s);
      double Dt = DU(ut, idx_bw_t);
      // P_N is a count: its cube is exact (= std::pow(P_N, 3))
      const double PN3 = PN_lag * PN_lag * PN_lag;
      double Ds_t = num_Ds / PN3;
      double Dt_s = num_Dt / PN3;
      double first_dependence_term = 15 * (sqrt(Ds_t * Dt) + sqrt(Dt_s * Ds) + 3 * sqrt(Ds * Dt)) / PN_lag;
      double second_dependence_term = (XsXtvar + sumabslr / PN_lag) / PN_lag;
      double dependence_term = first_dependence_term + second_dependence_term;
      double autocov_risk = bias_term + variance_term + dependence_term;
      if (full) {
        double* M = mat_res_risk.memptr();
        const arma::uword nr = n_rows_res;
        const double vals[14] = {s(k), t(k), bw_s, bw_t, PN_lag, h0, Hs, Ls, Ht, Lt,
                                 bias_term, variance_term, dependence_term, autocov_risk};
        for (int col = 0; col < 14; ++col) M[row + col * nr] = vals[col];
      }
      if (std::isfinite(autocov_risk)) {
        if (n_finite == 0 || autocov_risk < best_risk) { best_risk = autocov_risk; best_rank = n_finite; }
        ++n_finite;
      }
    };

    if (common_bw) {
      for (int b = 0; b < B; ++b) {
        const double bw = bw_grid_to_use(b);
        const double bw_s_pow = std::pow(bw, 2 * Hs), bw_t_pow = std::pow(bw, 2 * Ht);
        const std::size_t os = ((std::size_t) us * B + b) * nc, ot = ((std::size_t) ut * B + b) * nc;
        const double *sp = &c_pn[os], *ss = &c_sbn[os], *sm = &c_max[os];
        const double *tp = &c_pn[ot + lag], *tsb = &c_sbn[ot + lag], *tm = &c_max[ot + lag];
        double PN_lag = 0, B1 = 0, B2 = 0, V1 = 0, V2 = 0, V3 = 0;
#pragma omp simd reduction(+:PN_lag,B1,B2,V1,V2,V3)
        for (int i = 0; i < nn; ++i) {
          PN_lag += sp[i] * tp[i];
          B1 += tp[i] * ss[i];
          B2 += sp[i] * tsb[i];
          V1 += tp[i] * sm[i];
          V2 += sp[i] * tm[i];
          V3 += sm[i] * tm[i];
        }
        double bias_term_num = mom_t_square * Ls * bw_s_pow * B1 + mom_s_square * Lt * bw_t_pow * B2;
        double variance_term_num = sig_s_square * mom_t_square * V1 + sig_t_square * mom_s_square * V2 +
          sig_s_square * sig_t_square * V3;
        finish(bw, bw, b, b, PN_lag, bias_term_num, variance_term_num, (arma::uword) b * n + k);
      }
    } else if (nested_s[us] && nested_t[ut]) {
      // P_N, B_1, C_1 by cumulative sums along the t-side entry order, B_2 along
      // the s-side entry order; C_3 as a dense product.
      std::vector<double> PNm(B * B), B1m(B * B), C1m(B * B), B2m(B * B);
      for (int bs = 0; bs < B; ++bs) {
        const std::size_t os = ((std::size_t) us * B + bs) * nc;
        const double *sp = &c_pn[os], *ss = &c_sbn[os], *sm = &c_max[os];
        double a_pn = 0, a_b1 = 0, a_c1 = 0;
        for (int pr = 0; pr < B; ++pr) {
          for (int q = start_t[ut][pr]; q < start_t[ut][pr + 1]; ++q) {
            const int i = entry_t[ut][q];
            a_pn += sp[i]; a_b1 += ss[i]; a_c1 += sm[i];
          }
          const int bt = ord[pr];
          PNm[bs * B + bt] = a_pn; B1m[bs * B + bt] = a_b1; C1m[bs * B + bt] = a_c1;
        }
      }
      for (int bt = 0; bt < B; ++bt) {
        const double* tsb = &c_sbn[((std::size_t) ut * B + bt) * nc + lag];
        double a_b2 = 0;
        for (int pr = 0; pr < B; ++pr) {
          for (int q = start_s[us][pr]; q < start_s[us][pr + 1]; ++q) a_b2 += tsb[entry_s[us][q]];
          B2m[ord[pr] * B + bt] = a_b2;
        }
      }
      std::vector<double> pow_t(B);
      for (int bt = 0; bt < B; ++bt) pow_t[bt] = std::pow(bw_grid_to_use(bt), 2 * Ht);
      for (int bs = 0; bs < B; ++bs) {
        const double bw_s = bw_grid_to_use(bs);
        const double bw_s_pow = std::pow(bw_s, 2 * Hs);
        const double* sm = &c_max[((std::size_t) us * B + bs) * nc];
        for (int bt = 0; bt < B; ++bt) {
          const double bw_t = bw_grid_to_use(bt);
          const double bw_t_pow = pow_t[bt];
          const double* tm = &c_max[((std::size_t) ut * B + bt) * nc + lag];
          // NB: this branch historically squares pn_t (not pn_s * pn_t) in the
          // variance term; preserved: C2 runs over every curve with pn_t = 1.
          const double C2 = tmax_sum[ut * B + bt];
          double C3 = 0;
#pragma omp simd reduction(+:C3)
          for (int i = 0; i < nn; ++i) C3 += sm[i] * tm[i];
          const int ix = bs * B + bt;
          double bias_term_num = mom_t_square * Ls * bw_s_pow * B1m[ix] + mom_s_square * Lt * bw_t_pow * B2m[ix];
          double variance_term_num = sig_s_square * mom_t_square * C1m[ix] + sig_t_square * mom_s_square * C2 +
            sig_s_square * sig_t_square * C3;
          finish(bw_s, bw_t, bs, bt, PNm[ix], bias_term_num, variance_term_num,
                 (arma::uword) bs * B * n + (arma::uword) bt * n + k);
        }
      }
    } else {
      for (int bs = 0; bs < B; ++bs) {
        const double bw_s = bw_grid_to_use(bs);
        const double bw_s_pow = std::pow(bw_s, 2 * Hs);
        const std::size_t os = ((std::size_t) us * B + bs) * nc;
        for (int bt = 0; bt < B; ++bt) {
          const double bw_t = bw_grid_to_use(bt);
          const double bw_t_pow = std::pow(bw_t, 2 * Ht);
          const std::size_t ot = ((std::size_t) ut * B + bt) * nc;
          // NB: this branch historically squares pn_t (not pn_s * pn_t) in the
          // variance term; preserved: C2 runs over every curve with pn_t = 1.
          const double C2 = tmax_sum[ut * B + bt];
          const double *sp = &c_pn[os], *ss = &c_sbn[os], *sm = &c_max[os];
          const double *tp = &c_pn[ot + lag], *tsb = &c_sbn[ot + lag], *tm = &c_max[ot + lag];
          double PN_lag = 0, B1 = 0, B2 = 0, C1 = 0, C3 = 0;
#pragma omp simd reduction(+:PN_lag,B1,B2,C1,C3)
          for (int i = 0; i < nn; ++i) {
            PN_lag += sp[i] * tp[i];
            B1 += tp[i] * ss[i];
            B2 += sp[i] * tsb[i];
            C1 += tp[i] * sm[i];
            C3 += sm[i] * tm[i];
          }
          double bias_term_num = mom_t_square * Ls * bw_s_pow * B1 + mom_s_square * Lt * bw_t_pow * B2;
          double variance_term_num = sig_s_square * mom_t_square * C1 + sig_t_square * mom_s_square * C2 +
            sig_s_square * sig_t_square * C3;
          finish(bw_s, bw_t, bs, bt, PN_lag, bias_term_num, variance_term_num,
                 (arma::uword) bs * B * n + (arma::uword) bt * n + k);
        }
      }
    }
    if (best) {
      if (n_finite == 0) any_missing = true;
      if (common_bw) {
        best->bw_s(k) = bw_grid_to_use(best_rank); best->bw_t(k) = bw_grid_to_use(best_rank);
      } else {
        best->bw_s(k) = bw_grid_to_use(best_rank / B); best->bw_t(k) = bw_grid_to_use(best_rank % B);
      }
      best->Hs(k) = Hs; best->Ls(k) = Ls; best->Ht(k) = Ht; best->Lt(k) = Lt;
    }
  }
  if (best && any_missing) Rcpp::stop("index_min(): object has no elements");
  return mat_res_risk;
}

//' Estimate the risk of the covariance or autocovariance function
 //'
 //' Estimate the risk function of the lag-\eqn{\ell}, \eqn{\ell} = 0, 1,..., autocovariance function estimator of
 //' \insertCite{maissoro2024adaptive;textual}{adaptiveFTS} and \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
 //'
 //' @param data A DataFrame containing the columns "id_curve", "tobs", and "X". Typically, the output of the function \link{format_data}.
 //' @param s \code{vector (numeric)}. First argument of the autocovariance function.
 //' It corresponds to the observation points \code{s} in the pair (\code{s}, \code{t}).
 //' It has to be of the same length as the \code{t}
 //' @param t \code{vector (numeric)}. Second argument of the autocovariance function.
 //' It corresponds to the observation points \code{t} in the pair (\code{s}, \code{t}).
 //' It has to be of the same length as the \code{s}.
 //' @param lag \code{integer (positive integer)}. Lag of the autocovariance.
 //' @param bw_grid \code{vector (numeric)}. The bandwidth grid in which the best smoothing parameter is selected for each pair (\code{s}, \code{t}).
 //' It can be \code{NULL} and that way it will be defined as an exponential grid of \eqn{N\times\lambda}.
 //' @param common_bw A logical value indicating if the same bandwidth should be used for \code{s} and \code{t}. Default is \code{false}.
 //' @param center A logical value indicating if the data should be centered before estimation. Default is \code{true}.
 //' @param kernel_name \code{string}. Specifies the kernel function for estimation; default is "epanechnikov".
 //' @param presmooth_bw Numeric (positive vector or scalar). Bandwidth used to presmooth
 //' each curve in the local regularity step, see `estimate_locreg_cpp`. Default
 //' \code{NULL} selects it by cross-validation.
 //' @param Delta Numeric (positive). Length of the neighborhood of each point used in the
 //' local regularity step. Default \code{NULL} estimates it from the data.
 //' @param presmooth_bw_grid Numeric vector. Candidate bandwidths of the cross-validation
 //' that selects \code{presmooth_bw}. Default \code{NULL} uses the default grid.
 //' @param presmooth_nsubset Integer (positive). Number of curves used by that
 //' cross-validation. Default \code{NULL} uses min(70, floor(N / 2)) curves.
 //' Supported kernels include: "epanechnikov", "biweight", "triweight", "tricube", "triangular", and "uniform".
 //'
 //' @return A \code{matrix} containing the following fourteen columns in order:
 //'          \itemize{
 //'            \item{s : The first argument of the autocovariance function.}
 //'            \item{t : The second argument of the autocovariance function.}
 //'            \item{hs : The candidate bandwidth for the first argument of the autocovariance function. If \code{common_bw = TRUE}, the same bandwidth candidate is used for \code{s} and for \code{t}, so the 3rd and 4th columns contain the same values.}
 //'            \item{ht : The candidate bandwidth for the second argument of the autocovariance function.}
 //'            \item{PNl : The number of curves used in the estimation the autocovariance at (s,t). It corresponds to \eqn{P_{N,\ell}(s,t;h_s, h_t)}.}
 //'            \item{locreg_bw : The bandwidth used to estimate the local regularity parameters.}
 //'            \item{Hs : The estimates of the local exponent for each t. It corresponds to \eqn{H_s}.}
 //'            \item{Ls2 : The estimates of the Hölder constant for each t. It corresponds to \eqn{L_s^2}.}
 //'            \item{Ht : The estimates of the local exponent for each t. It corresponds to \eqn{H_t}.}
 //'            \item{Lt2 : The estimates of the Hölder constant for each t. It corresponds to \eqn{L_t^2}.}
 //'            \item{bias_term : The bias term of the risk function.}
 //'            \item{variance_term : The variance term of the risk function.}
 //'            \item{dependence_term : The dependence term of the risk function.}
 //'            \item{autocov_risk : The estimates of the risk function of the covariance/autocovariance function.}
 //'         }
 //'
 //'
 //' @import Rdpack
 //'
 //' @references
 //' \insertAllCited{}
 //'
 //' @examples
 //' \dontrun{
 //' # Load required libraries
 //' library(Rcpp)
 //' library(data.table)
 //'
 //' # Sample data
 //' set.seed(123)
 //' n <- 100
 //' data <- data.table(
 //'   id_curve = rep(1:10, each = n),
 //'   tobs = rep(seq(0, 1, length.out = n), times = 10),
 //'   X = rnorm(1000)
 //' )
 //'
 //' # Observation points
 //' s <- seq(0, 1, length.out = 50)
 //' t <- seq(0, 1, length.out = 50)
 //'
 //' # Estimate autocovariance risk
 //' autocov_risk <- estimate_autocov_risk_cpp(data, s, t, 1)
 //'
 //' # Print the result
 //' print(autocov_risk)
 //' }
 //'
 //'
 // [[Rcpp::export]]
 arma::mat estimate_autocov_risk_cpp(const Rcpp::DataFrame data,
                                     const arma::vec s,
                                     const arma::vec t,
                                     const int lag,
                                     const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                     const bool common_bw = false,
                                     const bool center = true,
                                     const std::string kernel_name = "epanechnikov",
                                     const Rcpp::Nullable<arma::vec> presmooth_bw = R_NilValue,
                                     const Rcpp::Nullable<double> Delta = R_NilValue,
                                     const Rcpp::Nullable<arma::vec> presmooth_bw_grid = R_NilValue,
                                     const Rcpp::Nullable<int> presmooth_nsubset = R_NilValue){
   // NB : if common_bw = TRUE, then the same bandwidth is used for s and t
   //      otherwise, a different bandwidth is use for s and for t.
   if (s.size() != t.size()) {
     stop("Arguments 's' and 't' must be of equal length.");
   } else if (s.size() == 0) {
     stop("Arguments 's' and 't' must be a numeric vectors or scalar value(s) between 0 and 1.");
   }
   core::Kern kern = core::parse_kernel(kernel_name);
   core::Curves C = core::make_curves(data);
   arma::vec bw_grid_to_use = bw_grid.isNull() ? default_autocov_bw_grid(C) : Rcpp::as<arma::vec>(bw_grid);
   LocregCfg cfg = resolve_locreg(C, presmooth_bw, Delta, presmooth_bw_grid, presmooth_nsubset, kern);
   return autocov_risk_impl(C, s, t, lag, bw_grid_to_use, common_bw, center, kern, cfg, true, nullptr, nullptr);
 }

 //' Get Unique Pairs of Elements in Upper Triangular Form
 //'
 //' This function takes two vectors `s` and `t` of the same length and generates
 //' unique pairs (s(k), t(k)) in an upper triangular form. This means for each
 //' pair, the smaller element is the first component and the larger element is
 //' the second component. Duplicate pairs are removed, and the result is sorted
 //' by the first and then the second column.
 //'
 //' @param s A numeric vector.
 //' @param t A numeric vector of the same length as `s`.
 //'
 //' @return A matrix where each row contains a unique pair (s(k), t(k)) in
 //' upper triangular form, sorted by the first and then the second column.
 //'
 //' @examples
 //' \dontrun{
 //'   library(RcppArmadillo)
 //'   s <- c(0.2, 0.4, 0.6, 0.4)
 //'   t <- c(0.4, 0.2, 0.6, 0.6)
 //'   get_upper_tri_couple(s, t)
 //' }
 //'
 // [[Rcpp::export]]
 arma::mat get_upper_tri_couple(const arma::vec& s, const arma::vec& t) {
   // Ensure vectors s and t have the same length
   if (s.size() != t.size()) {
     Rcpp::stop("Vectors s and t must have the same length.");
   }

   int n = s.size();
   arma::mat mat_st(n, 2);

   // Construct the matrix of pairs (s(k), t(k))
   for (int k = 0; k < n; ++k) {
     mat_st.row(k) = {std::min(s(k), t(k)), std::max(s(k), t(k))};
   }

   // Remove duplicated rows, keying on the EXACT (min, max) doubles so the
   // couple set stays consistent with the exact `==` look-ups downstream (the
   // lower-triangle fill in estimate_autocov_cpp). NB: do not key on a rounded
   // string here -- it can collapse distinct points and drop a couple, causing
   // an out-of-bounds look-up later.
   std::set<std::pair<double, double>> seen;
   std::vector<arma::uword> unique_indices;

   for (arma::uword i = 0; i < mat_st.n_rows; ++i) {
     std::pair<double, double> key(mat_st(i, 0), mat_st(i, 1));
     if (seen.insert(key).second) {
       unique_indices.push_back(i);
     }
   }

   arma::mat unique_mat(unique_indices.size(), mat_st.n_cols);
   for (arma::uword i = 0; i < unique_indices.size(); ++i) {
     unique_mat.row(i) = mat_st.row(unique_indices[i]);
   }

   // Order by first and second column
   // Step 1: Create a combined key for sorting
   arma::vec combined_key = unique_mat.col(0) * 1e6 + unique_mat.col(1);

   // Step 2: Get the sort indices based on the combined key
   arma::uvec sort_indices = arma::sort_index(combined_key);

   // Step 3: Reorder the rows of the matrix using the sort indices
   arma::mat sorted_mat = unique_mat.rows(sort_indices);

   return sorted_mat;
 }

 //' Sort Matrix by Specified Columns
 //'
 //' This function sorts a matrix by two specified columns in sequence.
 //'
 //' @param mat A numeric matrix to be sorted.
 //' @param first_col_idx An integer indicating the first column index to sort by.
 //' @param second_col_idx An integer indicating the second column index to sort by.
 //' @return A matrix sorted first by the column indicated by \code{first_col_idx} and then by the column indicated by \code{second_col_idx}.
 //' @export
 //'
 //' @examples
 //' \dontrun{
 //' mat <- matrix(c(3, 2, 5, 1, 3, 4, 2, 2, 3, 1, 1, 2), ncol = 3, byrow = TRUE)
 //' sorted_mat <- sort_by_columns(mat, 1, 2)
 //' print(sorted_mat)
 //' }
 // [[Rcpp::export]]
 arma::mat sort_by_columns(const arma::mat& mat, arma::uword first_col_idx, arma::uword second_col_idx) {
   // Order by first and second column
   // Step 1: Create a combined key for sorting
   arma::vec combined_key = mat.col(first_col_idx) * 1e6 + mat.col(second_col_idx);

   // Step 2: Get the sort indices based on the combined key
   arma::uvec sort_indices = arma::sort_index(combined_key);

   // Step 3: Reorder the rows of the matrix using the sort indices
   arma::mat sorted_mat = mat.rows(sort_indices);

   return sorted_mat;
 }

// The (auto)covariance estimator. The bandwidths come from the risk (fused
// argmin) when `has_bw` is false. Each curve is smoothed once on the distinct
// (point, bandwidth) tuples; the diagonal correction sums W_s W_t over the
// intersection of the two kernel windows.
arma::mat autocov_impl(const core::Curves& C, const arma::vec& s, const arma::vec& t, const int lag,
                       const bool has_bw, const arma::vec& bw_s_in, const arma::vec& bw_t_in,
                       const arma::vec& bw_grid, const bool common_bw, const bool center,
                       const bool correct_diagonal, core::Kern kern, const LocregCfg& cfg,
                       const AutocovPoints* shared) {
  const int n_curve = C.N;

  // Vector to be used in the estimation procedure
  arma::vec svec, tvec;
  if (lag == 0) {
    arma::mat mat_st_upper = get_upper_tri_couple(s, t);
    svec = mat_st_upper.col(0);
    tvec = mat_st_upper.col(1);
  } else {
    svec = s;
    tvec = t;
  }
  int n = tvec.n_elem;

  arma::vec optbw_s_to_use(n), optbw_t_to_use(n);
  arma::mat mat_locreg(n, 6);

  if (!has_bw) {
    AutocovBest best;
    autocov_risk_impl(C, svec, tvec, lag, bw_grid, common_bw, center, kern, cfg, false, &best, shared);
    optbw_s_to_use = best.bw_s;
    optbw_t_to_use = best.bw_t;
    mat_locreg.col(0) = svec;
    mat_locreg.col(1) = tvec;
    mat_locreg.col(2) = best.Hs;
    mat_locreg.col(3) = best.Ls;
    mat_locreg.col(4) = best.Ht;
    mat_locreg.col(5) = best.Lt;
  } else {
    int optbw_s_to_use_temp_size = bw_s_in.size();
    int optbw_t_to_use_temp_size = bw_t_in.size();

    if ((optbw_s_to_use_temp_size == n) && (optbw_t_to_use_temp_size == n)) {
      optbw_s_to_use = bw_s_in;
      optbw_t_to_use = bw_t_in;
    } else if ((optbw_s_to_use_temp_size == 1) && (optbw_t_to_use_temp_size == 1)){
      optbw_s_to_use = arma::ones(n) * (bw_s_in(0));
      optbw_t_to_use = arma::ones(n) * (bw_t_in(0));
    } else if (lag == 0) {
      // Bandwidth of the first pair holding the value in s (else in t)
      std::map<double, arma::uword> first_s, first_t;
      for (arma::uword i = 0; i < s.n_elem; ++i) first_s.emplace(s(i), i);
      for (arma::uword i = 0; i < t.n_elem; ++i) first_t.emplace(t(i), i);
      for (int j = 0; j < n; ++j) {
        auto is = first_s.find(svec(j));
        arma::uword i_s = (is != first_s.end()) ? is->second : first_t.at(svec(j));
        auto it = first_t.find(tvec(j));
        arma::uword i_t = (it != first_t.end()) ? it->second : first_s.at(tvec(j));
        optbw_s_to_use(j) = bw_s_in(i_s);
        optbw_t_to_use(j) = bw_t_in(i_t);
      }
    } else {
      stop("If 'bw_s' and 'bw_t' are not NULL, they must be the same length as 's' and as 't' or of length 1.");
    }
    mat_locreg.col(0) = svec;
    mat_locreg.col(1) = tvec;
    mat_locreg.cols(2, 5).zeros();
  }

  // Deduplicate (point, bandwidth) tuples: each estimate below depends on a
  // single (s, h_s) or (t, h_t) tuple. Compute each distinct tuple once.
  std::map<std::pair<double, double>, int> umap_s, umap_t;
  arma::uvec map_s(n), map_t(n);
  std::vector<double> svec_u_v, obs_u_v, tvec_u_v, obt_u_v;
  for (int k = 0; k < n; ++k) {
    std::pair<double, double> kys(svec(k), optbw_s_to_use(k));
    auto it = umap_s.find(kys);
    if (it == umap_s.end()) {
      int id = (int) svec_u_v.size();
      umap_s.emplace(kys, id); svec_u_v.push_back(svec(k)); obs_u_v.push_back(optbw_s_to_use(k));
      map_s(k) = id;
    } else { map_s(k) = it->second; }
    std::pair<double, double> kyt(tvec(k), optbw_t_to_use(k));
    auto jt = umap_t.find(kyt);
    if (jt == umap_t.end()) {
      int id = (int) tvec_u_v.size();
      umap_t.emplace(kyt, id); tvec_u_v.push_back(tvec(k)); obt_u_v.push_back(optbw_t_to_use(k));
      map_t(k) = id;
    } else { map_t(k) = jt->second; }
  }
  arma::vec svec_u(svec_u_v), obs_u(obs_u_v), tvec_u(tvec_u_v), obt_u(obt_u_v);
  int Us = svec_u.n_elem, Ut = tvec_u.n_elem;

  // Each curve smoothed once on each distinct tuple: estimate, support
  // indicator (|T - x| <= b) and kernel sum, stored curve-contiguous per tuple.
  // The same values give the means (all curves) and the cross products (curves
  // n and n + lag), as estimate_mean_cpp and the former per-curve loop did.
  for (int i = 0; i < n_curve; ++i)
    if (!C.Y[i].is_finite() || !C.T[i].is_finite())
      stop("The vectors y and t must be non-empty numeric vectors of the same length.");
  arma::mat XS(n_curve, Us), PS(n_curve, Us), SS(n_curve, Us);
  arma::mat XT(n_curve, Ut), PT(n_curve, Ut), ST(n_curve, Ut);
  auto smooth_tuples = [&](const arma::vec& pts, const arma::vec& bws, int nu, arma::mat& X, arma::mat& P, arma::mat& S) {
#pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < nu; ++j) {
      for (int c = 0; c < n_curve; ++c) {
        double sk, pn;
        double x = core::nw_pn_abs(C, c, pts(j), bws(j), kern, sk, pn);
        X(c, j) = std::isfinite(x) ? x : 0;
        P(c, j) = pn;
        S(c, j) = sk;
      }
    }
  };
  smooth_tuples(svec_u, obs_u, Us, XS, PS, SS);
  smooth_tuples(tvec_u, obt_u, Ut, XT, PT, ST);

  // Mean function at each tuple (P_N and mu_N over all the curves)
  auto tuple_means = [&](const arma::mat& X, const arma::mat& P, int nu, arma::vec& PNu, arma::vec& muu) {
    PNu.set_size(nu); muu.set_size(nu);
    for (int j = 0; j < nu; ++j) {
      double pn_sum = 0, num = 0;
      for (int c = 0; c < n_curve; ++c) { pn_sum += P(c, j); num += P(c, j) * X(c, j); }
      PNu(j) = pn_sum; muu(j) = num / pn_sum;
    }
  };
  arma::vec PN_s_u, mu_s_u, PN_t_u, mu_t_u;
  tuple_means(XS, PS, Us, PN_s_u, mu_s_u);
  tuple_means(XT, PT, Ut, PN_t_u, mu_t_u);
  arma::vec PN_s = PN_s_u.elem(map_s), muhat_s = mu_s_u.elem(map_s);
  arma::vec PN_t = PN_t_u.elem(map_t), muhat_t = mu_t_u.elem(map_t);

  // Error sd for the diagonal correction (once per distinct point)
  arma::vec vec_sig_s(n, arma::fill::zeros);
  arma::vec vec_sig_t(n, arma::fill::zeros);
  const bool diag = (lag == 0 && correct_diagonal);
  if (diag) {
    arma::vec sig_s_u = sigma_impl(C, svec_u).col(1);
    arma::vec sig_t_u = sigma_impl(C, tvec_u).col(1);
    vec_sig_s = sig_s_u.elem(map_s);
    vec_sig_t = sig_t_u.elem(map_t);
  }

  // Cross products, curve by curve (in order) for each pair. The diagonal
  // correction sums W_s W_t (W = K / sum K, 0 when the sum is 0) over the
  // intersection of the two kernel windows of the curve (lag 0: same curve).
  arma::vec PN_lag(n, arma::fill::zeros);
  arma::vec autocov_numerator(n, arma::fill::zeros);
  arma::vec diag_correct_numerator(n, arma::fill::zeros);
#pragma omp parallel for schedule(dynamic)
  for (int k = 0; k < n; ++k) {
    const arma::uword js = map_s(k), jt = map_t(k);
    const double xs0 = svec_u(js), bs0 = obs_u(js), xt0 = tvec_u(jt), bt0 = obt_u(jt);
    double pnl = 0, num = 0, dnum = 0;
    for (int i = 0; i < n_curve - lag; ++i) {
      const int cs = i, ct = i + lag;
      const double ps = PS(cs, js), pt = PT(ct, jt);
      const double xs = XS(cs, js), xt = XT(ct, jt);
      pnl += ps * pt;
      if (center) {
        num += ps * pt * (xs - muhat_s(k)) * (xt - muhat_t(k));
      } else {
        num += ps * pt * xs * xt;
      }
      if (diag && ps * pt != 0) {
        const double Ss = SS(cs, js), St = ST(ct, jt);
        double wp = 0;
        if (Ss != 0 && St != 0 && std::isfinite(Ss) && std::isfinite(St)) {
          arma::uword ls, hs, lt, ht;
          core::window(C, cs, xs0, bs0, ls, hs);
          core::window(C, ct, xt0, bt0, lt, ht);
          const std::vector<double>& ts = C.Ts[cs];
          for (arma::uword q = std::max(ls, lt); q < std::min(hs, ht); ++q)
            wp += (core::kval(kern, (ts[q] - xs0) / bs0) / Ss) * (core::kval(kern, (ts[q] - xt0) / bt0) / St);
        }
        dnum += vec_sig_s(k) * vec_sig_t(k) * ps * pt * wp;
      }
    }
    PN_lag(k) = pnl;
    autocov_numerator(k) = num;
    diag_correct_numerator(k) = dnum;
  }

  // Compute autocovariance estimate
  // // Note that if lag = 0 or correct_diagonal = FALSE, Then by definition diag_correct_numerator = 0
  arma::vec autocovhat;
  if (center) {
    autocovhat = autocov_numerator / PN_lag - diag_correct_numerator / PN_lag;
  } else {
    autocovhat = autocov_numerator / PN_lag - muhat_s % muhat_t - diag_correct_numerator / PN_lag;
  }

  // Return the result
  int n_couple = t.size();
  arma::mat mat_res_autocov(n_couple, 14);

  // // For autocovariance Or the upper part of the coavariance
  arma::mat res_autocov(n, 14);
  res_autocov.col(0) = svec;
  res_autocov.col(1) = tvec;
  res_autocov.col(2) = optbw_s_to_use;
  res_autocov.col(3) = optbw_t_to_use;
  res_autocov.col(4) = mat_locreg.col(2);
  res_autocov.col(5) = mat_locreg.col(3);
  res_autocov.col(6) = mat_locreg.col(4);
  res_autocov.col(7) = mat_locreg.col(5);
  res_autocov.col(8) = PN_s;
  res_autocov.col(9) = muhat_s;
  res_autocov.col(10) = PN_t;
  res_autocov.col(11) = muhat_t;
  res_autocov.col(12) = PN_lag;
  res_autocov.col(13) = autocovhat;

  if (lag == 0) {
    // Fill the lower part of the covariance matrix from the (unique) upper pairs
    std::map<std::pair<double, double>, arma::uword> row_of;
    for (int r = 0; r < n; ++r) row_of.emplace(std::make_pair(svec(r), tvec(r)), r);
    for (int k = 0; k < n_couple ; ++k) {
      auto up = row_of.find(std::make_pair(s(k), t(k)));
      if (up != row_of.end()) {
        mat_res_autocov.row(k) = res_autocov.row(up->second);
      } else {
        arma::uword r = row_of.at(std::make_pair(t(k), s(k)));
        mat_res_autocov.row(k) = {res_autocov(r, 1), res_autocov(r, 0),
                                  res_autocov(r, 3), res_autocov(r, 2),
                                  res_autocov(r, 6), res_autocov(r, 7),
                                  res_autocov(r, 4), res_autocov(r, 5),
                                  res_autocov(r, 10), res_autocov(r, 11),
                                  res_autocov(r, 8), res_autocov(r, 9),
                                  res_autocov(r, 12), res_autocov(r, 13)};
      }
    }
  } else {
    mat_res_autocov = res_autocov;
  }

  return sort_by_columns(mat_res_autocov, 0, 1);
}

 //' Estimate the autocovariance function
 //'
 //' This function estimates the lag-\eqn{\ell}, \eqn{\ell} = 0, 1,..., autocovariance function based on the methodology of
 //' \insertCite{maissoro2024adaptive;textual}{adaptiveFTS} and \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
 //'
 //' @param data A \code{DataFrame} containing the data with columns \code{"id_curve"}, \code{"tobs"}, and \code{"X"}.
 //' @param s A numeric vector specifying time points \code{s} for which to estimate autocovariance.
 //' @param t A numeric vector specifying time points \code{t} for which to estimate autocovariance.
 //' @param lag An integer specifying the lag value for autocovariance.
 //' @param bw_s Optional numeric vector specifying optimal bandwidths for \code{s}. Default is \code{NULL}.
 //' @param bw_t Optional numeric vector specifying optimal bandwidths for \code{t}. Default is \code{NULL}.
 //' @param bw_grid numeric vector of bandwidth grid values. Default is \code{NULL}.
 //' @param common_bw A logical value indicating if the same bandwidth should be used for \code{s} and \code{t}. Default is \code{false}.
 //' @param center A logical value indicating if the data should be centered before estimation. Default is \code{true}.
 //' @param correct_diagonal A logical value indicating whether the diagonal of the covariance should be corrected when \code{lag=0}.
 //' @param kernel_name A string specifying the kernel to use for estimation. Supported values are \code{"epanechnikov"}, \code{"biweight"}, \code{"triweight"}, \code{"tricube"}, \code{"triangular"}, \code{"uniform"}. Default is \code{"epanechnikov"}.
 //' @param presmooth_bw Numeric (positive vector or scalar). Bandwidth used to presmooth
 //' each curve in the local regularity step, see `estimate_locreg_cpp`. Default
 //' \code{NULL} selects it by cross-validation.
 //' @param Delta Numeric (positive). Length of the neighborhood of each point used in the
 //' local regularity step. Default \code{NULL} estimates it from the data.
 //' @param presmooth_bw_grid Numeric vector. Candidate bandwidths of the cross-validation
 //' that selects \code{presmooth_bw}. Default \code{NULL} uses the default grid.
 //' @param presmooth_nsubset Integer (positive). Number of curves used by that
 //' cross-validation. Default \code{NULL} uses min(70, floor(N / 2)) curves.
 //'
 //' @return A \code{matrix} containing the following fourteen columns in order:
 //'          \itemize{
 //'            \item{s : The first argument of the autocovariance function.}
 //'            \item{t : The second argument of the autocovariance function.}
 //'            \item{bw_s : The optimal bandwidth for the first argument of the autocovariance function. If \code{common_bw = TRUE}, the same bandwidth candidate is used for \code{s} and for \code{t}, so the 3rd and 4th columns contain the same values.}
 //'            \item{bw_t : The optimal bandwidth for the second argument of the autocovariance function.}
 //'            \item{Hs : The estimates of the local exponent for each \code{s}. It corresponds to \eqn{H_s}.}
 //'            \item{Ls2 : The estimates of the Hölder constant for each \code{s}. It corresponds to \eqn{L_s^2}.}
 //'            \item{Ht : The estimates of the local exponent for each \code{t}. It corresponds to \eqn{H_t}.}
 //'            \item{Lt2 : The estimates of the Hölder constant for each \code{t}. It corresponds to \eqn{L_t^2}.}
 //'            \item{PNs : The number of curves used to estimate the mean at \code{s}. It corresponds to \eqn{P_N(s;h)}.}
 //'            \item{muhat_s : The estimates of the mean at \code{s}. It corresponds to \eqn{\widehat{\mu}_N(s;h)}.}
 //'            \item{PNt : The number of curves used to estimate the mean at \code{t}. It corresponds to \eqn{P_N(t;h)}.}
 //'            \item{muhat_t : The estimates of the mean at \code{t}. It corresponds to \eqn{\widehat{\mu}_N(t;h)}.}
 //'            \item{PNl : The number of curves used to estimate the autocovariance at \code{(s,t)}. It corresponds to \eqn{P_{N,\ell}(s,t;h_s, h_t)}.}
 //'            \item{autocov : The estimates of the covariance/autocovariance.}
 //'         }
 //'
 //' @examples
 //' \dontrun{
 //' # Example usage
 //' data <- data.frame(id_curve = rep(1:10, each = 100),
 //'                    tobs = rep(seq(0, 1, length.out = 100), 10),
 //'                    X = rnorm(1000))
 //' s <- seq(0, 1, length.out = 10)
 //' t <- seq(0, 1, length.out = 10)
 //' result <- estimate_autocov_cpp(data, s, t, lag = 0)
 //' }
 //'
 //' @export
 //'
 //' @import Rdpack
 //'
 //' @references
 //' \insertAllCited{}
 //'
 // [[Rcpp::export]]
 arma::mat estimate_autocov_cpp(const Rcpp::DataFrame data,
                                const arma::vec s,
                                const arma::vec t,
                                const int lag,
                                const Rcpp::Nullable<arma::vec> bw_s = R_NilValue,
                                const Rcpp::Nullable<arma::vec> bw_t = R_NilValue,
                                const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                const bool common_bw = false,
                                const bool center = true,
                                const bool correct_diagonal = true,
                                const std::string kernel_name = "epanechnikov",
                                const Rcpp::Nullable<arma::vec> presmooth_bw = R_NilValue,
                                const Rcpp::Nullable<double> Delta = R_NilValue,
                                const Rcpp::Nullable<arma::vec> presmooth_bw_grid = R_NilValue,
                                const Rcpp::Nullable<int> presmooth_nsubset = R_NilValue) {
   if (s.size() != t.size()) {
     stop("Arguments 's' and 't' must be of equal length.");
   } else if (s.size() == 0) {
     stop("Arguments 's' and 't' must be numeric vectors or scalar value(s) between 0 and 1.");
   }
   core::Kern kern = core::parse_kernel(kernel_name);
   core::Curves C = core::make_curves(data);
   if (bw_s.isNull() || bw_t.isNull()) {
     arma::vec bw_grid_to_use = bw_grid.isNull() ? default_autocov_bw_grid(C) : Rcpp::as<arma::vec>(bw_grid);
     LocregCfg cfg = resolve_locreg(C, presmooth_bw, Delta, presmooth_bw_grid, presmooth_nsubset, kern);
     return autocov_impl(C, s, t, lag, false, arma::vec(), arma::vec(), bw_grid_to_use, common_bw,
                         center, correct_diagonal, kern, cfg, nullptr);
   }
   return autocov_impl(C, s, t, lag, true, Rcpp::as<arma::vec>(bw_s), Rcpp::as<arma::vec>(bw_t),
                       arma::vec(), common_bw, center, correct_diagonal, kern, LocregCfg(), nullptr);
 }
