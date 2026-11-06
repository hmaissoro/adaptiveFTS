#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "02_smoothing_rcpp.h"
#include "00_internal.h"
using namespace Rcpp;
using namespace arma;

// =============================================================================
// C++ core of the design-weighted, Tikhonov-regularised adaptive functional
// BLUP (R/11_blup.R).
//
// Building block: every plug-in estimate is computed directly at the points it
// is needed at, with its own risk-minimising bandwidth (as in the prediction
// scenario of M. Vimond): the mean at the design T_d and at the prediction
// points, the covariance C0 on T_d x T_d, the lag-1 autocovariance C1 on
// T_d x t_pred. No bandwidth sub-grid / nearest-neighbour look-up is involved.
// `blup_fit_cpp` estimates the prediction-point-independent part (C0, mean,
// noise, regularised variance matrix V); `blup_predict_cpp` evaluates the
// one-step-ahead BLUP and loops for h-step-ahead prediction; `blup_operators_cpp`
// returns the operators of one fold of the Tikhonov cross-validation.
// =============================================================================

namespace {

// Median over the finite entries (matches R's median(x, na.rm = TRUE)).
double median_finite(const arma::vec& v) {
  arma::vec f = v.elem(arma::find_finite(v));
  double med = arma::median(f);
  return med;
}

// Nearest PSD matrix: symmetrise and floor the eigenvalues at zero. The kernel
// covariance estimate C0 is not guaranteed PSD; without this the regularised
// system A0 + tikhonov*I stays indefinite for small tikhonov.
arma::mat psd_project(const arma::mat& M) {
  arma::vec eval;
  arma::mat evec;
  arma::eig_sym(eval, evec, 0.5 * (M + M.t()));
  eval = arma::clamp(eval, 0.0, arma::datum::inf);
  arma::mat out = evec * arma::diagmat(eval) * evec.t();
  out = 0.5 * (out + out.t());
  return out;
}

// Design weights for a design `Td`: 1/M under the common design; the
// importance-sampling correction 1/(M * ghat) with a fixed-bandwidth
// leave-one-out Parzen-Rosenblatt estimate under the independent design.
arma::vec compute_rho(const arma::vec& Td, bool is_common, double density_bw, core::Kern kern) {
  arma::uword M = Td.n_elem;
  if (is_common) {
    arma::vec w(M);
    w.fill(1.0 / M);
    return w;
  }
  arma::vec ghat(M);
  for (arma::uword i = 0; i < M; ++i) {
    double sum = 0;
    for (arma::uword j = 0; j < M; ++j) if (j != i) sum += core::kval(kern, (Td(i) - Td(j)) / density_bw);
    ghat(i) = sum / ((M - 1) * density_bw);
  }
  arma::vec rho = 1.0 / (M * arma::clamp(ghat, 1e-6, arma::datum::inf));
  rho /= arma::sum(rho);
  return rho;
}

// Adaptive mean at `pts` (bandwidth selected by the risk at each point).
void mean_direct(const core::Curves& C, const arma::vec& pts, const arma::vec& bwg,
                 core::Kern kern, const LocregCfg& cfg, arma::vec& muhat, arma::vec& bw) {
  arma::mat m = mean_impl(C, pts, arma::vec(), bwg, kern, cfg);
  muhat = m.col(5);
  bw = m.col(1);
}

// Adaptive (auto)covariance block c_lag(s_i, t_j), bandwidths (h_s, h_t)
// selected by the risk at each pair.
void cov_direct(const core::Curves& C, const arma::vec& s, const arma::vec& t, int lag,
                bool correct_diagonal, const arma::vec& bwg, core::Kern kern, const LocregCfg& cfg,
                arma::mat& block, arma::mat& bw_s, arma::mat& bw_t, const AutocovPoints* shared = nullptr) {
  const arma::uword ns = s.n_elem, nt = t.n_elem;
  arma::vec gs(ns * nt), gt(ns * nt);
  for (arma::uword j = 0; j < nt; ++j)
    for (arma::uword i = 0; i < ns; ++i) {
      gs(j * ns + i) = s(i);
      gt(j * ns + i) = t(j);
    }
  arma::mat out = autocov_impl(C, gs, gt, lag, false, arma::vec(), arma::vec(), bwg,
                               false, true, correct_diagonal, kern, cfg, shared);
  std::map<std::pair<double, double>, arma::uword> row_of;
  for (arma::uword r = 0; r < out.n_rows; ++r) row_of.emplace(std::make_pair(out(r, 0), out(r, 1)), r);
  block.set_size(ns, nt); bw_s.set_size(ns, nt); bw_t.set_size(ns, nt);
  for (arma::uword i = 0; i < ns; ++i)
    for (arma::uword j = 0; j < nt; ++j) {
      arma::uword r = row_of.at(std::make_pair(s(i), t(j)));
      block(i, j) = out(r, 13); bw_s(i, j) = out(r, 2); bw_t(i, j) = out(r, 3);
    }
}

// Noise variance at Td: pointwise, or the median when homoscedastic.
arma::vec noise_variance(const core::Curves& C, const arma::vec& Td, bool homoscedastic) {
  arma::vec sig2 = arma::square(sigma_impl(C, Td).col(1));
  if (!homoscedastic) return sig2;
  arma::vec s2v(1);
  s2v(0) = median_finite(sig2);
  return s2v;
}

// Regularised variance matrix V = D C0 D + diag(sigma^2 rho) + tikhonov I.
arma::mat variance_matrix(const arma::mat& c0, const arma::vec& sigma2, const arma::vec& rho,
                          double tikhonov) {
  const arma::uword M = rho.n_elem;
  arma::mat root_D = arma::diagmat(arma::sqrt(rho));
  arma::vec noise = (sigma2.n_elem == 1) ? arma::vec(sigma2(0) * rho) : arma::vec(sigma2 % rho);
  return root_D * c0 * root_D + arma::diagmat(noise) + tikhonov * arma::eye(M, M);
}

LocregCfg cached_cfg(const core::Curves& C, const arma::vec& presmooth_bw, double Delta, core::Kern kern) {
  // `Rcpp::Nullable` stores a bare SEXP without protecting it, so the wrapped
  // values are held in Rcpp types that keep them alive across the call.
  Rcpp::NumericVector h_r = Rcpp::wrap(presmooth_bw);
  Rcpp::NumericVector d_r = Rcpp::NumericVector::create(Delta);
  return resolve_locreg(C, Rcpp::Nullable<arma::vec>((SEXP) h_r), Rcpp::Nullable<double>((SEXP) d_r),
                        R_NilValue, R_NilValue, kern);
}

} // anonymous namespace

//' Project a matrix onto the PSD cone (C++ core)
//'
//' Symmetrises `M` and floors its eigenvalues at zero. Called by
//' `select_tikhonov_parameter()`; not intended to be used directly.
//'
//' @param M A numeric matrix.
//' @return The nearest positive-semidefinite matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::mat psd_project_cpp(const arma::mat M) {
  arma::mat out = psd_project(M);
  return out;
}

//' Plug-in operators of the adaptive BLUP at given points (C++ core)
//'
//' Estimates, each at its own risk-minimising bandwidth, the mean at the design
//' `Td` and at the prediction points `tpred`, the covariance on `Td x Td`, the
//' pointwise noise variance at `Td` and the lag-1 autocovariance on
//' `Td x tpred`. Called by `select_tikhonov_parameter()`; not intended to be
//' used directly.
//'
//' @param data A DataFrame with columns \code{id_curve}, \code{tobs}, \code{X}.
//' @param Td Design points (sorted).
//' @param tpred Prediction points (sorted); may be empty.
//' @param bw_grid Bandwidth grid of the adaptive risks.
//' @param kernel_name Kernel name.
//' @param presmooth_bw,Delta,presmooth_bw_grid,presmooth_nsubset Local
//'   regularity settings, see `estimate_locreg_cpp`.
//' @return A list with `muhat_d`, `c0` (not projected), `sigma2` (pointwise),
//'   `muhat_pred` and `c1`.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List blup_operators_cpp(const Rcpp::DataFrame data, const arma::vec Td, const arma::vec tpred,
                              const arma::vec bw_grid, const std::string kernel_name,
                              const Rcpp::Nullable<arma::vec> presmooth_bw = R_NilValue,
                              const Rcpp::Nullable<double> Delta = R_NilValue,
                              const Rcpp::Nullable<arma::vec> presmooth_bw_grid = R_NilValue,
                              const Rcpp::Nullable<int> presmooth_nsubset = R_NilValue) {
  core::Kern kern = core::parse_kernel(kernel_name);
  core::Curves C = core::make_curves(data);
  LocregCfg cfg = resolve_locreg(C, presmooth_bw, Delta, presmooth_bw_grid, presmooth_nsubset, kern);
  // One mean risk on the union of the points (each point is independent)
  arma::vec U = arma::unique(arma::join_cols(Td, tpred));
  arma::vec mu_U, bw_U;
  mean_direct(C, U, bw_grid, kern, cfg, mu_U, bw_U);
  auto at = [&U, &mu_U](const arma::vec& x) {
    arma::vec out(x.n_elem);
    for (arma::uword i = 0; i < x.n_elem; ++i)
      out(i) = mu_U(std::lower_bound(U.begin(), U.end(), x(i)) - U.begin());
    return out;
  };
  // C0 and C1 share their point constants (computed once on Td and tpred)
  AutocovPoints pts = autocov_points(C, U, bw_grid, true, kern, cfg);
  arma::mat c0, b0s, b0t;
  cov_direct(C, Td, Td, 0, true, bw_grid, kern, cfg, c0, b0s, b0t, &pts);
  arma::vec sig2 = arma::square(sigma_impl(C, Td).col(1));
  arma::mat c1, b1s, b1t;
  if (tpred.n_elem > 0) cov_direct(C, Td, tpred, 1, false, bw_grid, kern, cfg, c1, b1s, b1t, &pts);
  return Rcpp::List::create(
    Rcpp::Named("muhat_d") = at(Td),
    Rcpp::Named("c0") = c0,
    Rcpp::Named("sigma2") = sig2,
    Rcpp::Named("muhat_pred") = at(tpred),
    Rcpp::Named("c1") = c1);
}

//' Fit the adaptive functional BLUP (C++ core)
//'
//' Estimates every prediction-point-independent component of the adaptive
//' Best Linear Unbiased Predictor, with the bandwidths selected directly at the
//' conditioning-curve design points. Called by the R function \code{blup_fit()};
//' not intended to be used directly.
//'
//' @param data A DataFrame with columns \code{id_curve}, \code{tobs}, \code{X}.
//' @param id_conditioning_curve Integer id of the conditioning curve.
//' @param bw_grid Bandwidth grid for the adaptive risk.
//' @param rho Design weights of the conditioning curve.
//' @param homoscedastic Whether to use a constant noise variance.
//' @param tikhonov Tikhonov regularisation parameter.
//' @param kernel_name Kernel name.
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
//' @return A list with the selected bandwidths, the covariance operator, the
//'   mean, the noise level, the regularised variance matrix, the residual and
//'   the resolved local regularity settings.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List blup_fit_cpp(const Rcpp::DataFrame data,
                        const int id_conditioning_curve,
                        const arma::vec bw_grid,
                        const arma::vec rho,
                        const bool homoscedastic,
                        const double tikhonov,
                        const std::string kernel_name,
                        const Rcpp::Nullable<arma::vec> presmooth_bw = R_NilValue,
                        const Rcpp::Nullable<double> Delta = R_NilValue,
                        const Rcpp::Nullable<arma::vec> presmooth_bw_grid = R_NilValue,
                        const Rcpp::Nullable<int> presmooth_nsubset = R_NilValue) {
  core::Kern kern = core::parse_kernel(kernel_name);
  core::Curves C = core::make_curves(data);
  LocregCfg cfg = resolve_locreg(C, presmooth_bw, Delta, presmooth_bw_grid, presmooth_nsubset, kern);

  // Conditioning curve, ordered by observation time.
  arma::uword c0i = 0;
  while (c0i < (arma::uword) C.N && C.ids(c0i) != id_conditioning_curve) ++c0i;
  if (c0i == (arma::uword) C.N) stop("'id_conditioning_curve' is not a curve of 'data'.");
  arma::uvec ord = arma::sort_index(C.T[c0i]);
  arma::vec Tn0 = arma::sort(C.T[c0i]);
  arma::vec Yn0 = C.Y[c0i].elem(ord);

  arma::vec muhat, bw_mean;
  mean_direct(C, Tn0, bw_grid, kern, cfg, muhat, bw_mean);
  arma::mat c0raw, bw_s, bw_t;
  cov_direct(C, Tn0, Tn0, 0, true, bw_grid, kern, cfg, c0raw, bw_s, bw_t);
  arma::mat c0 = psd_project(c0raw);
  arma::vec sigma2 = noise_variance(C, Tn0, homoscedastic);
  arma::mat V = variance_matrix(c0, sigma2, rho, tikhonov);
  arma::vec resid = arma::diagmat(arma::sqrt(rho)) * (Yn0 - muhat);

  return Rcpp::List::create(
    Rcpp::Named("Tn0") = Tn0,
    Rcpp::Named("Yn0") = Yn0,
    Rcpp::Named("bw_mean") = bw_mean,
    Rcpp::Named("bw_cov_s") = bw_s,
    Rcpp::Named("bw_cov_t") = bw_t,
    Rcpp::Named("muhat_Tn0") = muhat,
    Rcpp::Named("c0hat") = c0,
    Rcpp::Named("sigma2") = sigma2,
    Rcpp::Named("V") = V,
    Rcpp::Named("resid") = resid,
    Rcpp::Named("presmooth_bw") = cfg.h,
    Rcpp::Named("Delta") = cfg.Delta);
}

//' Predict with the adaptive functional BLUP (C++ core)
//'
//' Evaluates the adaptive BLUP at the prediction points, looping for
//' h-step-ahead prediction. The mean at the prediction points and the lag-1
//' autocovariance are estimated directly at those points. Called by
//' \code{predict.blup_fit()}; not intended to be used directly.
//'
//' @param data A DataFrame with columns \code{id_curve}, \code{tobs}, \code{X}.
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
//' @param bw_grid Bandwidth grid of the adaptive risks.
//' @param presmooth_bw,Delta Local regularity settings resolved by the fit.
//'
//' @return A matrix with columns \code{horizon}, \code{t}, \code{muhat}, \code{prediction}.
//' @keywords internal
// [[Rcpp::export]]
arma::mat blup_predict_cpp(const Rcpp::DataFrame data,
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
                           const std::string kernel_name,
                           const arma::vec bw_grid,
                           const arma::vec presmooth_bw,
                           const double Delta) {
  core::Kern kern = core::parse_kernel(kernel_name);
  core::Curves C = core::make_curves(data);
  LocregCfg cfg = cached_cfg(C, presmooth_bw, Delta, kern);
  arma::vec tpred = arma::sort(arma::unique(t));
  arma::uword nt = tpred.n_elem;

  // One row block per horizon step, stacked: [horizon, t, muhat, prediction].
  arma::mat out(horizon * nt, 4);

  arma::vec muhat_t, bw_mu;
  mean_direct(C, tpred, bw_grid, kern, cfg, muhat_t, bw_mu);
  arma::mat c1, b1s, b1t;
  cov_direct(C, Tn0, tpred, 1, false, bw_grid, kern, cfg, c1, b1s, b1t);
  arma::vec resid = root_D * (Yn0 - muhat_Tn0);
  arma::vec pred = muhat_t + c1.t() * root_D * arma::solve(V, resid);
  out.submat(0, 0, nt - 1, 0).fill(1.0);
  out.submat(0, 1, nt - 1, 1) = tpred;
  out.submat(0, 2, nt - 1, 2) = muhat_t;
  out.submat(0, 3, nt - 1, 3) = pred;

  // Feed each predicted curve (on `tpred`) back as the new conditioning curve;
  // the estimation data is left unchanged, only the conditioning values vary,
  // so the operators on tpred are estimated once.
  if (horizon >= 2) {
    arma::vec rho = compute_rho(tpred, is_common, density_bw, kern);
    arma::mat rD = arma::diagmat(arma::sqrt(rho));
    AutocovPoints pts = autocov_points(C, tpred, bw_grid, true, kern, cfg);
    arma::mat c0raw, b0s, b0t;
    cov_direct(C, tpred, tpred, 0, true, bw_grid, kern, cfg, c0raw, b0s, b0t, &pts);
    arma::mat Vp = variance_matrix(psd_project(c0raw), noise_variance(C, tpred, homoscedastic), rho, tikhonov);
    arma::mat c1p, bps, bpt;
    cov_direct(C, tpred, tpred, 1, false, bw_grid, kern, cfg, c1p, bps, bpt, &pts);
    arma::mat gain = c1p.t() * rD;
    for (int step = 2; step <= horizon; ++step) {
      arma::vec r = rD * (pred - muhat_t);
      pred = muhat_t + gain * arma::solve(Vp, r);
      arma::uword r0 = (step - 1) * nt;
      out.submat(r0, 0, r0 + nt - 1, 0).fill(static_cast<double>(step));
      out.submat(r0, 1, r0 + nt - 1, 1) = tpred;
      out.submat(r0, 2, r0 + nt - 1, 2) = muhat_t;
      out.submat(r0, 3, r0 + nt - 1, 3) = pred;
    }
  }

  return out;
}
