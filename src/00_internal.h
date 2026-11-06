#ifndef ADAPTIVEFTS_INTERNAL_H
#define ADAPTIVEFTS_INTERNAL_H

// Internal (non-exported) entry points of the estimators, working on the curves
// split once (core::Curves). The exported *_cpp functions are thin wrappers
// around them, and the estimators call each other through them, so that the
// data are split, the kernel resolved and the presmoothing bandwidth selected
// once per call instead of once per sub-estimator.

#include <RcppArmadillo.h>
#include "00_core.h"

// 02_smoothing_rcpp.cpp
double nw_bw_impl(const arma::vec& y, const arma::vec& t, const arma::vec& bw_grid, core::Kern kern);
double nw_optimal_bw_impl(const core::Curves& C, const Rcpp::Nullable<arma::vec> bw_grid,
                          const Rcpp::Nullable<int> nsubset, core::Kern kern);

// 03_estimate_locreg_rcpp.cpp: presmoothing bandwidths (one per curve) and
// neighbourhood length, resolved once and shared by every sub-estimator.
struct LocregCfg {
  arma::vec h;
  double Delta = 0;
};
LocregCfg resolve_locreg(const core::Curves& C, const Rcpp::Nullable<arma::vec> h,
                         const Rcpp::Nullable<double> Delta,
                         const Rcpp::Nullable<arma::vec> presmooth_bw_grid,
                         const Rcpp::Nullable<int> presmooth_nsubset, core::Kern kern);
arma::mat locreg_impl(const core::Curves& C, const arma::vec& t, const bool center,
                      core::Kern kern, const LocregCfg& cfg);

// 07_estimate_constants_cpp.cpp
arma::mat smooth_matrix(const core::Curves& C, const arma::vec& t, const arma::vec& h, core::Kern kern);
arma::mat sigma_impl(const core::Curves& C, const arma::vec& t);
arma::mat mom_impl(const core::Curves& C, const arma::vec& t, const arma::vec& h,
                   const double mom_order, const double center, core::Kern kern);
arma::mat emp_autocov_matrix(const core::Curves& C, const arma::vec& t, const arma::vec& h,
                             const arma::vec& lag, core::Kern kern);
void xsxt_matrices(const core::Curves& C, const arma::vec& s, const arma::vec& t,
                   const int cross_lag, const arma::vec& lag, const arma::vec& h,
                   core::Kern kern, const bool center, arma::vec& gamma, arma::mat& acov);
void xsxt_dependence(const core::Curves& C, const arma::vec& s, const arma::vec& t,
                     const int cross_lag, const arma::vec& lag, const arma::vec& h,
                     core::Kern kern, const bool center, arma::vec& var0, arma::vec& lrsum);
void dd_matrices(const core::Curves& C, const arma::vec& t, const arma::vec& bw_vec,
                 const arma::vec& h, const int max_lag, core::Kern kern, const bool center,
                 arma::mat& num, arma::mat& PN);

// 04_estimate_mean_cpp.cpp
arma::vec default_mean_bw_grid(const core::Curves& C);
arma::mat mean_risk_impl(const core::Curves& C, const arma::vec& t, const arma::vec& bw_grid,
                         core::Kern kern, const LocregCfg& cfg);
arma::mat mean_impl(const core::Curves& C, const arma::vec& t, const arma::vec& bw,
                    const arma::vec& bw_grid, core::Kern kern, const LocregCfg& cfg);

// 05_estimate_autocov_cpp.cpp
struct AutocovBest {          // first minimiser of the risk, per pair
  arma::vec bw_s, bw_t, Hs, Ls, Ht, Lt;
};
// Point constants of the (auto)covariance risk on the distinct points U (local
// regularity, noise level, second moments, D numerators, per-curve cubes on
// (point, bandwidth, curve)); they do not depend on the pairs nor on the lag,
// so several risks on the same points can share them.
struct AutocovPoints {
  arma::vec U, bw_grid;
  arma::mat locU, sigU, momU, ddnum, ddPN;
  double h0 = 0;
  int B = 0;
  bool center = true;
  std::vector<double> c_pn, c_sbn, c_max;
};
AutocovPoints autocov_points(const core::Curves& C, const arma::vec& pts, const arma::vec& bw_grid,
                             const bool center, core::Kern kern, const LocregCfg& cfg);
arma::vec default_autocov_bw_grid(const core::Curves& C);
arma::mat autocov_risk_impl(const core::Curves& C, const arma::vec& s, const arma::vec& t,
                            const int lag, const arma::vec& bw_grid, const bool common_bw,
                            const bool center, core::Kern kern, const LocregCfg& cfg,
                            const bool full, AutocovBest* best, const AutocovPoints* shared = nullptr);
arma::mat autocov_impl(const core::Curves& C, const arma::vec& s, const arma::vec& t, const int lag,
                       const bool has_bw, const arma::vec& bw_s, const arma::vec& bw_t,
                       const arma::vec& bw_grid, const bool common_bw, const bool center,
                       const bool correct_diagonal, core::Kern kern, const LocregCfg& cfg,
                       const AutocovPoints* shared = nullptr);
arma::mat get_upper_tri_couple(const arma::vec& s, const arma::vec& t);
arma::mat sort_by_columns(const arma::mat& mat, arma::uword first_col_idx, arma::uword second_col_idx);

// 06_estimate_cov_segment_cpp.cpp
arma::mat cov_segment_risk_impl(const core::Curves& C, const arma::vec& t, const arma::vec& bw_grid,
                                const bool center, core::Kern kern, const LocregCfg& cfg);

#endif
