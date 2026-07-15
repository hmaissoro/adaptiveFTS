#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "02_smoothing_rcpp.h"
#include "04_estimate_mean_cpp.h"
#include "05_estimate_autocov_cpp.h"
#include "07_estimate_constants_cpp.h"
using namespace Rcpp;
using namespace arma;

// =============================================================================
// C++ core of the design-weighted, Tikhonov-regularised adaptive functional
// BLUP (R/11_blup.R). `blup_fit_cpp` estimates every prediction-point-
// independent component (adaptive bandwidths, C0, mean, noise, regularised
// variance matrix V); `blup_predict_cpp` evaluates the one-step-ahead BLUP and
// loops for h-step-ahead prediction. Both reuse the compiled estimation
// primitives (estimate_mean(_risk)_cpp, estimate_autocov(_risk)_cpp,
// estimate_sigma_cpp) so the numerics match the R reference bit-for-bit.
//
// This translation unit is self-contained (it does not depend on
// 08_estimate_curve_cpp) so the old block-reconstruction predictor can be
// removed without touching the new one.
// =============================================================================

namespace {

// Index of the minimum among the finite entries of `v`, mapped back to the
// original position. Matches R's which.min() (first minimum, NAs ignored).
arma::uword argmin_finite(const arma::vec& v) {
  arma::uvec fin = arma::find_finite(v);
  return fin(arma::index_min(v.elem(fin)));
}

// Nearest-neighbour indices in one dimension (squared distance, first argmin).
arma::uvec nn1_1d(const arma::vec& ref, const arma::vec& query) {
  arma::uvec idx(query.n_elem);
  for (arma::uword i = 0; i < query.n_elem; ++i)
    idx(i) = arma::index_min(arma::square(ref - query(i)));
  return idx;
}

// Nearest-neighbour indices in two dimensions (Euclidean, first argmin).
arma::uvec nn1_2d(const arma::vec& ref_s, const arma::vec& ref_t,
                  const arma::vec& q_s, const arma::vec& q_t) {
  arma::uvec idx(q_s.n_elem);
  for (arma::uword i = 0; i < q_s.n_elem; ++i)
    idx(i) = arma::index_min(arma::square(ref_s - q_s(i)) + arma::square(ref_t - q_t(i)));
  return idx;
}

// Reshape a long (s, t, value) table to a [unique-s x unique-t] matrix, with
// both axes sorted ascending (as R's dcast / the package's reshape_matrix do).
arma::mat reshape_long(const arma::mat& A, arma::uword cs, arma::uword ct, arma::uword cv) {
  arma::vec svec = arma::sort(arma::unique(A.col(cs)), "ascend");
  arma::vec tvec = arma::sort(arma::unique(A.col(ct)), "ascend");
  arma::mat out(svec.n_elem, tvec.n_elem, arma::fill::zeros);
  for (arma::uword i = 0; i < svec.n_elem; ++i)
    for (arma::uword j = 0; j < tvec.n_elem; ++j) {
      arma::uvec k = arma::find((A.col(cs) == svec(i)) % (A.col(ct) == tvec(j)));
      if (!k.is_empty()) out(i, j) = A(k(0), cv);
    }
  return out;
}

// Per-t optimal mean bandwidth: argmin of the risk (col 9) over the candidate
// bandwidths (col 1), for each sub-grid point. Returns [t, optbw].
arma::mat select_mean_optbw(const arma::mat& risk, const arma::vec& tgrid) {
  arma::mat out(tgrid.n_elem, 2);
  out.col(0) = tgrid;
  for (arma::uword k = 0; k < tgrid.n_elem; ++k) {
    arma::uvec rows = arma::find(risk.col(0) == tgrid(k));
    arma::uword m = argmin_finite(risk(rows, arma::uvec({9})));
    out(k, 1) = risk(rows(m), 1);
  }
  return out;
}

// Per-(s,t) optimal (auto)covariance bandwidths: argmin of the risk (col 13)
// over the candidate bandwidths (optbw_s col 2, optbw_t col 3). Returns
// [s, t, optbw_s, optbw_t].
arma::mat select_autocov_optbw(const arma::mat& risk, const arma::vec& sgrid, const arma::vec& tgrid) {
  arma::mat out(sgrid.n_elem, 4);
  out.col(0) = sgrid;
  out.col(1) = tgrid;
  for (arma::uword k = 0; k < sgrid.n_elem; ++k) {
    arma::uvec rows = arma::find((risk.col(0) == sgrid(k)) % (risk.col(1) == tgrid(k)));
    arma::uword m = argmin_finite(risk(rows, arma::uvec({13})));
    out(k, 2) = risk(rows(m), 2);
    out(k, 3) = risk(rows(m), 3);
  }
  return out;
}

// Mean evaluated at `tq` (assumed sorted) reusing the cached bandwidths.
arma::vec mean_at(const DataFrame& data, const arma::mat& opt_mean,
                  const arma::vec& tq, const std::string& kernel) {
  arma::uvec idx = nn1_1d(opt_mean.col(0), tq);
  arma::vec optbw_all = opt_mean.col(1);
  arma::vec optbw = optbw_all.elem(idx);
  arma::mat m = estimate_mean_cpp(data, tq, Rcpp::wrap(optbw), R_NilValue, kernel);
  return m.col(5);
}

// (Auto)covariance block c_lag(s_i, t_j) reusing the cached bandwidths.
arma::mat autocov_at(const DataFrame& data, const arma::mat& opt_bw,
                     const arma::vec& s, const arma::vec& t, int lag,
                     bool correct_diagonal, const std::string& kernel) {
  arma::uword ns = s.n_elem, nt = t.n_elem;
  arma::vec gs(ns * nt), gt(ns * nt);
  for (arma::uword j = 0; j < nt; ++j)
    for (arma::uword i = 0; i < ns; ++i) {
      gs(j * ns + i) = s(i);
      gt(j * ns + i) = t(j);
    }
  arma::uvec idx = nn1_2d(opt_bw.col(0), opt_bw.col(1), gs, gt);
  arma::vec obs_all = opt_bw.col(2);
  arma::vec obt_all = opt_bw.col(3);
  arma::vec obs = obs_all.elem(idx);
  arma::vec obt = obt_all.elem(idx);
  arma::mat out = estimate_autocov_cpp(data, gs, gt, lag, Rcpp::wrap(obs), Rcpp::wrap(obt),
                                       R_NilValue, false, true, correct_diagonal, kernel);
  return reshape_long(out, 0, 1, 13);
}

// Design weights for a design `Td`: 1/M under the common design; the
// importance-sampling correction 1/(M * ghat) with a fixed-bandwidth
// leave-one-out Parzen-Rosenblatt estimate under the independent design.
arma::vec compute_rho(const arma::vec& Td, bool is_common, double density_bw,
                      std::function<arma::vec(const arma::vec)>& kfun) {
  arma::uword M = Td.n_elem;
  if (is_common) {
    arma::vec w(M);
    w.fill(1.0 / M);
    return w;
  }
  arma::vec ghat(M);
  arma::uvec all = arma::regspace<arma::uvec>(0, M - 1);
  for (arma::uword i = 0; i < M; ++i) {
    arma::uvec others = all.elem(arma::find(all != i));
    arma::vec u = (Td(i) - Td.elem(others)) / density_bw;
    ghat(i) = arma::sum(kfun(u)) / ((M - 1) * density_bw);
  }
  arma::vec rho = 1.0 / (M * arma::clamp(ghat, 1e-6, arma::datum::inf));
  return rho / arma::sum(rho);
}

// Median over the finite entries (matches R's median(x, na.rm = TRUE)).
double median_finite(const arma::vec& v) {
  arma::vec f = v.elem(arma::find_finite(v));
  return arma::median(f);
}

struct Cond {
  arma::vec muhat;    // mean at the design
  arma::mat c0;       // symmetrised covariance
  arma::vec sigma2;   // noise variance (length 1 if homoscedastic)
  arma::mat root_D;   // diag(sqrt(rho))
  arma::mat V;        // regularised variance matrix
  arma::vec resid;    // root_D (Yd - muhat)
};

// Conditioning quantities for an arbitrary design `Td` with values `Yd`.
Cond condition(const DataFrame& data, const arma::mat& opt_mean, const arma::mat& opt_cov,
               const arma::vec& Td, const arma::vec& Yd, const arma::vec& rho,
               bool homoscedastic, double tikhonov, const std::string& kernel) {
  Cond c;
  arma::uword M = Td.n_elem;
  c.root_D = arma::diagmat(arma::sqrt(rho));
  c.muhat = mean_at(data, opt_mean, Td, kernel);
  arma::mat c0 = autocov_at(data, opt_cov, Td, Td, 0, true, kernel);
  c.c0 = 0.5 * (c0 + c0.t());
  arma::vec sig2 = arma::square(estimate_sigma_cpp(data, Td).col(1));
  arma::vec noise;
  if (homoscedastic) {
    double s2 = median_finite(sig2);
    arma::vec s2v(1);
    s2v(0) = s2;
    c.sigma2 = s2v;
    noise = s2 * rho;
  } else {
    c.sigma2 = sig2;
    noise = sig2 % rho;
  }
  c.V = c.root_D * c.c0 * c.root_D + arma::diagmat(noise) + tikhonov * arma::eye(M, M);
  c.resid = c.root_D * (Yd - c.muhat);
  return c;
}

// One-step BLUP at `tpred` given conditioning quantities; also returns the mean
// at `tpred` through `muhat_out`.
arma::vec blup_one_step(const DataFrame& data, const arma::mat& opt_mean,
                        const arma::mat& opt_autocov, const std::string& kernel,
                        const arma::vec& tpred, const arma::vec& Td, const arma::mat& rootD,
                        const arma::mat& Vmat, const arma::vec& resid, arma::vec& muhat_out) {
  muhat_out = mean_at(data, opt_mean, tpred, kernel);
  arma::mat c1 = autocov_at(data, opt_autocov, Td, tpred, 1, false, kernel);
  return muhat_out + c1.t() * rootD * arma::solve(Vmat, resid);
}

} // anonymous namespace

//' Mean at new locations using cached bandwidths (C++ core)
//'
//' Reuses the adaptive mean bandwidths cached in a `blup_fit` object. Called by
//' `select_tikhonov_parameter()`; not intended to be used directly.
//'
//' @param data A DataFrame with columns \code{id_curve}, \code{tobs}, \code{X}.
//' @param opt_mean Cached mean adaptive-bandwidth matrix (`t`, `optbw`).
//' @param t Evaluation locations (assumed sorted).
//' @param kernel_name Kernel name.
//' @return The mean estimates at `t`.
//' @keywords internal
// [[Rcpp::export]]
arma::vec blup_mean_at_cpp(const Rcpp::DataFrame data, const arma::mat opt_mean,
                           const arma::vec t, const std::string kernel_name) {
  arma::vec muhat = mean_at(data, opt_mean, t, kernel_name);
  return muhat;
}

//' (Auto)covariance block at new locations using cached bandwidths (C++ core)
//'
//' Reuses the adaptive (auto)covariance bandwidths cached in a `blup_fit`
//' object. Called by `select_tikhonov_parameter()`; not intended to be used directly.
//'
//' @param data A DataFrame with columns \code{id_curve}, \code{tobs}, \code{X}.
//' @param opt_bw Cached (auto)covariance bandwidth matrix (`s`, `t`, `optbw_s`,
//'   `optbw_t`).
//' @param s,t Evaluation locations (rows indexed by `s`, columns by `t`).
//' @param lag 0 for the covariance, 1 for the lag-1 autocovariance.
//' @param correct_diagonal Whether to correct the covariance diagonal.
//' @param kernel_name Kernel name.
//' @return A `length(s)` by `length(t)` matrix of \eqn{\hat c_{lag}(s_i, t_j)}.
//' @keywords internal
// [[Rcpp::export]]
arma::mat blup_autocov_at_cpp(const Rcpp::DataFrame data, const arma::mat opt_bw,
                              const arma::vec s, const arma::vec t, const int lag,
                              const bool correct_diagonal, const std::string kernel_name) {
  arma::mat block = autocov_at(data, opt_bw, s, t, lag, correct_diagonal, kernel_name);
  return block;
}

//' Fit the adaptive functional BLUP (C++ core)
//'
//' Estimates every prediction-point-independent component of the adaptive
//' Best Linear Unbiased Predictor and caches the adaptive bandwidths. Called by
//' the R function \code{blup_fit()}; not intended to be used directly.
//'
//' @param data A DataFrame with columns \code{id_curve}, \code{tobs}, \code{X}.
//' @param id_lag Integer id of the conditioning curve.
//' @param bw_grid Bandwidth grid for the adaptive risk.
//' @param rho Design weights of the conditioning curve.
//' @param homoscedastic Whether to use a constant noise variance.
//' @param tikhonov Tikhonov regularisation parameter.
//' @param sub_grid_length Number of points per axis of the bandwidth sub-grid.
//' @param kernel_name Kernel name.
//'
//' @return A list with the cached bandwidths, the covariance operator, the
//'   mean, the noise level, the regularised variance matrix and the residual.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List blup_fit_cpp(const Rcpp::DataFrame data,
                        const int id_lag,
                        const arma::vec bw_grid,
                        const arma::vec rho,
                        const bool homoscedastic,
                        const double tikhonov,
                        const int sub_grid_length,
                        const std::string kernel_name) {
  arma::mat data_mat(data.nrows(), 3);
  data_mat.col(0) = as<arma::vec>(data["id_curve"]);
  data_mat.col(1) = as<arma::vec>(data["tobs"]);
  data_mat.col(2) = as<arma::vec>(data["X"]);

  arma::uvec idx = arma::find(data_mat.col(0) == id_lag);
  arma::vec Tn0 = arma::sort(data_mat(idx, arma::uvec({1})));
  // Values ordered by observation time.
  arma::vec Traw = data_mat(idx, arma::uvec({1}));
  arma::vec Yraw = data_mat(idx, arma::uvec({2}));
  arma::uvec ord = arma::sort_index(Traw);
  arma::vec Yn0 = Yraw.elem(ord);

  arma::vec sub_vec = arma::linspace(0.05, 0.95, sub_grid_length);
  arma::uword ng = sub_vec.n_elem;
  arma::vec gs(ng * ng), gt(ng * ng);
  for (arma::uword j = 0; j < ng; ++j)
    for (arma::uword i = 0; i < ng; ++i) {
      gs(j * ng + i) = sub_vec(i);
      gt(j * ng + i) = sub_vec(j);
    }
  Rcpp::Nullable<arma::vec> bwg = Rcpp::wrap(bw_grid);

  arma::mat mean_risk = estimate_mean_risk_cpp(data, sub_vec, bwg, kernel_name);
  arma::mat opt_mean = select_mean_optbw(mean_risk, sub_vec);

  arma::mat cov_risk = estimate_autocov_risk_cpp(data, gs, gt, 0, bwg, false, true, kernel_name);
  arma::mat opt_cov = select_autocov_optbw(cov_risk, gs, gt);

  arma::mat autocov_risk = estimate_autocov_risk_cpp(data, gs, gt, 1, bwg, false, true, kernel_name);
  arma::mat opt_autocov = select_autocov_optbw(autocov_risk, gs, gt);

  Cond c = condition(data, opt_mean, opt_cov, Tn0, Yn0, rho, homoscedastic, tikhonov, kernel_name);

  return Rcpp::List::create(
    Rcpp::Named("opt_mean") = opt_mean,
    Rcpp::Named("opt_cov") = opt_cov,
    Rcpp::Named("opt_autocov") = opt_autocov,
    Rcpp::Named("Tn0") = Tn0,
    Rcpp::Named("Yn0") = Yn0,
    Rcpp::Named("muhat_Tn0") = c.muhat,
    Rcpp::Named("c0hat") = c.c0,
    Rcpp::Named("sigma2") = c.sigma2,
    Rcpp::Named("V") = c.V,
    Rcpp::Named("resid") = c.resid);
}

//' Predict with the adaptive functional BLUP (C++ core)
//'
//' Evaluates the adaptive BLUP at the prediction points, looping for
//' h-step-ahead prediction. Called by \code{predict.blup_fit()}; not intended
//' to be used directly.
//'
//' @param data A DataFrame with columns \code{id_curve}, \code{tobs}, \code{X}.
//' @param opt_mean,opt_cov,opt_autocov Cached adaptive-bandwidth matrices.
//' @param Tn0 Conditioning-curve design points.
//' @param muhat_Tn0 Mean at the conditioning-curve design points.
//' @param V Regularised variance matrix of the fit.
//' @param root_D Square-root design-weight matrix of the fit.
//' @param Yn0 Conditioning-curve values (possibly overriding the fit's).
//' @param density_bw Fixed design-density bandwidth.
//' @param is_common Whether the design is common across curves.
//' @param homoscedastic Whether to use a constant noise variance.
//' @param tikhonov Tikhonov regularisation parameter.
//' @param t Prediction points.
//' @param horizon Prediction horizon (steps ahead).
//' @param kernel_name Kernel name.
//'
//' @return A matrix with columns \code{t}, \code{muhat}, \code{prediction}.
//' @keywords internal
// [[Rcpp::export]]
arma::mat blup_predict_cpp(const Rcpp::DataFrame data,
                           const arma::mat opt_mean,
                           const arma::mat opt_cov,
                           const arma::mat opt_autocov,
                           const arma::vec Tn0,
                           const arma::vec muhat_Tn0,
                           const arma::mat V,
                           const arma::mat root_D,
                           const arma::vec Yn0,
                           const double density_bw,
                           const bool is_common,
                           const bool homoscedastic,
                           const double tikhonov,
                           const arma::vec t,
                           const int horizon,
                           const std::string kernel_name) {
  std::function<arma::vec(const arma::vec)> kfun = select_kernel(kernel_name);
  arma::vec tpred = arma::sort(arma::unique(t));
  arma::uword nt = tpred.n_elem;

  // One row block per horizon step, stacked: [horizon, t, muhat, prediction].
  arma::mat out(horizon * nt, 4);

  arma::vec muhat_t;
  arma::vec resid = root_D * (Yn0 - muhat_Tn0);
  arma::vec pred = blup_one_step(data, opt_mean, opt_autocov, kernel_name, tpred,
                                 Tn0, root_D, V, resid, muhat_t);
  out.submat(0, 0, nt - 1, 0).fill(1.0);
  out.submat(0, 1, nt - 1, 1) = tpred;
  out.submat(0, 2, nt - 1, 2) = muhat_t;
  out.submat(0, 3, nt - 1, 3) = pred;

  // Feed each predicted curve (on `tpred`) back as the new conditioning curve;
  // the estimation data is left unchanged, only the conditioning values vary.
  for (int step = 2; step <= horizon; ++step) {
    arma::vec rho = compute_rho(tpred, is_common, density_bw, kfun);
    Cond c = condition(data, opt_mean, opt_cov, tpred, pred, rho,
                       homoscedastic, tikhonov, kernel_name);
    pred = blup_one_step(data, opt_mean, opt_autocov, kernel_name, tpred,
                         tpred, c.root_D, c.V, c.resid, muhat_t);
    arma::uword r0 = (step - 1) * nt;
    out.submat(r0, 0, r0 + nt - 1, 0).fill(static_cast<double>(step));
    out.submat(r0, 1, r0 + nt - 1, 1) = tpred;
    out.submat(r0, 2, r0 + nt - 1, 2) = muhat_t;
    out.submat(r0, 3, r0 + nt - 1, 3) = pred;
  }

  return out;
}
