#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

#include "00_internal.h"
using namespace Rcpp;
using namespace arma;

// Default bandwidth grid of the mean risk.
arma::vec default_mean_bw_grid(const core::Curves& C) {
  double n_curve = C.N;
  double lambdahat = C.lambdahat();
  double b0 = std::pow(n_curve * lambdahat, -1 / (2 * 0.05 + 1)); // rate with minimum local exponent = 0.05
  double bK = 0.2; // rate with maximum local exponent = 1
  return arma::logspace(log10(b0), log10(bK), 20);
}

// Risk of the mean estimator on the bandwidth grid (Myriam Vimond's phase A
// restructuring, in C++). The local regularity, the noise level and the
// empirical autocovariances (t x lag matrix) are computed once; for each
// (h, t) the per-curve quantities pi_n, max W and b_n(2H) come from the kernel
// window of each curve, and the dependence term sum_l 2 gamma_l(t) rho_l(t;h)
// uses the counts c_l = #{n : pi_n pi_{n+l} = 1} over the active curves.
arma::mat mean_risk_impl(const core::Curves& C, const arma::vec& t, const arma::vec& bw_grid_to_use,
                         core::Kern kern, const LocregCfg& cfg) {
  const int n = t.size();
  const int n_curve = C.N;
  const int bw_size = bw_grid_to_use.size();

  // Estimate local regularity
  arma::mat mat_locreg = locreg_impl(C, t, true, kern, cfg);
  arma::vec h(n_curve, arma::fill::value(mat_locreg(0, 1))); // the presmoothing bandwidth

  // Estimate the error sd
  arma::mat mat_sig = sigma_impl(C, t);

  // Estimate the autocovariance, lags 0..N-1 (one row per t)
  arma::mat ac = emp_autocov_matrix(C, t, h, arma::regspace(0, n_curve - 1), kern);

  // Per-t constants (first occurrence of each value of t, as the look-ups did)
  std::vector<int> first = core::first_index(t);
  arma::vec Ht_k(n), Lt_k(n), sig2_k(n);
  for (int k = 0; k < n; ++k) {
    Ht_k(k) = mat_locreg(first[k], 4);
    Lt_k(k) = mat_locreg(first[k], 5);
    sig2_k(k) = mat_sig(first[k], 1) * mat_sig(first[k], 1);
  }

  arma::mat mat_res_risk(bw_size * n, 10);
  const int n_tasks = bw_size * n;
#pragma omp parallel for schedule(dynamic)
  for (int task = 0; task < n_tasks; ++task) {
    const int idx_bw = task / n, k = task % n;
    const double bw = bw_grid_to_use(idx_bw);
    const double Ht = Ht_k(k), Lt = Lt_k(k), sig_square = sig2_k(k);
    const double bw_pow = std::pow(bw, 2 * Ht);

    std::vector<int> active;
    active.reserve(n_curve);
    double bias_term_num = 0, variance_term_num = 0;
    for (int i = 0; i < n_curve; ++i) {
      core::RiskStats r = core::risk_stats(C, i, t(k), bw, Ht, kern);
      if (r.pn > 0) active.push_back(i);
      bias_term_num += Lt * bw_pow * r.pn * r.sbn;
      variance_term_num += sig_square * r.pn * r.wmax;
    }

    // P_N(t;h), bias and variance terms
    double PN = active.size();
    double bias_term = 2 * bias_term_num / PN;
    double variance_term = 2 * variance_term_num / (PN * PN);

    // Dependence term: lag 0, then 2 gamma_l rho_l with rho_l = c_l / P_N
    std::vector<int> cnt(n_curve, 0);
    const int na = active.size();
    for (int a = 0; a < na; ++a)
      for (int b = a + 1; b < na; ++b) ++cnt[active[b] - active[a]];
    const int k0 = first[k];
    double dep_sum = ac(k0, 0);
    for (int ell = 1; ell < n_curve; ++ell) {
      double rho = cnt[ell] / PN;
      dep_sum += 2 * ac(k0, ell) * rho;
    }
    double dependence_term = 2 * std::abs(dep_sum) / PN;
    double mean_risk = bias_term + variance_term + dependence_term;

    double* row = mat_res_risk.memptr();
    const arma::uword nr = mat_res_risk.n_rows, ri = idx_bw * n + k;
    const double vals[10] = {t(k), bw, PN, h(0), Ht, Lt, bias_term, variance_term, dependence_term, mean_risk};
    for (int col = 0; col < 10; ++col) row[ri + col * nr] = vals[col];
  }
  return mat_res_risk;
}

// Adaptive mean at t: the bandwidth minimising the risk at each t (first finite
// minimum) when `bw` is empty, else `bw`; then each curve is smoothed once at t.
arma::mat mean_impl(const core::Curves& C, const arma::vec& t, const arma::vec& bw,
                    const arma::vec& bw_grid, core::Kern kern, const LocregCfg& cfg) {
  const int n = t.size();
  arma::vec optbw_to_use(n);
  arma::vec Ht_used(n, arma::fill::zeros);
  arma::vec Lt_used(n, arma::fill::zeros);
  if (bw.n_elem == 0) {
    arma::mat mat_risk = mean_risk_impl(C, t, bw_grid, kern, cfg);
    const int bw_size = bw_grid.n_elem;
    std::vector<std::vector<arma::uword>> groups = core::same_value_groups(t);
    for (int k = 0; k < n; ++k) {
      // As the former selection: the position of the first minimum among the
      // finite risks, applied to the list of all the rows of t(k).
      std::vector<arma::uword> rows;
      double best = arma::datum::inf;
      arma::uword best_rank = 0, n_finite = 0;
      for (int b = 0; b < bw_size; ++b)
        for (arma::uword kk : groups[k]) {
          arma::uword row = b * n + kk;
          rows.push_back(row);
          double r = mat_risk(row, 9);
          if (std::isfinite(r)) {
            if (n_finite == 0 || r < best) { best = r; best_rank = n_finite; }
            ++n_finite;
          }
        }
      if (n_finite == 0) Rcpp::stop("index_min(): object has no elements");
      arma::uword best_row = rows[best_rank];
      optbw_to_use(k) = mat_risk(best_row, 1);
      Ht_used(k) = mat_risk(best_row, 4);
      Lt_used(k) = mat_risk(best_row, 5);
    }
  } else {
    if ((int) bw.n_elem != n) stop("If 'bw' is not NULL, it must be the same length as 't'.");
    optbw_to_use = bw;
  }

  // Compute the mean function
  arma::vec PN(n, fill::zeros);
  arma::vec mean_numerator(n, fill::zeros);
  for (int i = 0; i < C.N; ++i)
    if (!C.Y[i].is_finite() || !C.T[i].is_finite())
      stop("The vectors y and t must be non-empty numeric vectors of the same length.");
#pragma omp parallel for schedule(dynamic)
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < C.N; ++i) {
      double Xhat = core::nw_at(C, i, t(j), optbw_to_use(j), kern);
      if (!std::isfinite(Xhat)) Xhat = 0;
      double pn = core::pn_abs(C, i, t(j), optbw_to_use(j));
      PN(j) += pn;
      mean_numerator(j) += pn * Xhat;
    }
  }
  arma::vec muhat = mean_numerator / PN;

  arma::mat mat_res_mean(n, 6);
  mat_res_mean.col(0) = t;
  mat_res_mean.col(1) = optbw_to_use;
  mat_res_mean.col(2) = Ht_used;
  mat_res_mean.col(3) = Lt_used;
  mat_res_mean.col(4) = PN;
  mat_res_mean.col(5) = muhat;
  return mat_res_mean;
}

//' Estimate the risk of the mean function
 //'
 //' This function estimates the risk function \eqn{R_\mu(t;h)} of the mean function estimation proposed by \insertCite{maissoro2024adaptive;textual}{adaptiveFTS}.
 //'
 //' @param data A DataFrame containing the columns "id_curve", "tobs", and "X", typically the output of the function \link{format_data}.
 //' @param t A numeric vector of observation points at which we want to estimate the mean function of the underlying process.
 //' @param bw_grid A numeric vector representing the bandwidth grid in which the best smoothing parameter is selected for each t.
 //' It can be NULL, in which case it will be defined as an exponential grid of \eqn{N \times \lambda}.
 //' @param kernel_name A string specifying the kernel function of the Nadaraya-Watson estimator. Default is "epanechnikov".
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
 //' @return A matrix containing the following ten columns:
 //'          \itemize{
 //'            \item{t : The points at which the risk function is estimated.}
 //'            \item{h : The candidate bandwidth.}
 //'            \item{PN : The number of curves used to estimate the mean at t. It corresponds to \eqn{P_N(t;h)}.}
 //'            \item{locreg_bw : The bandwidth used to estimate the local regularity parameters.}
 //'            \item{Ht : The estimates of the local exponent for each t. It corresponds to \eqn{H_t}}
 //'            \item{Lt2 : The estimates of the Hölder constant for each t. It corresponds to \eqn{L_t^2}.}
 //'            \item{bias_term : The bias term of the risk function.}
 //'            \item{variance_term : The variance term of the risk function.}
 //'            \item{dependence_term : The dependence term of the risk function.}
 //'            \item{mean_risk : The estimates of the risk function of the mean.}
 //'         }
 //'
 //' @details
 //' The local regularity parameters are directly estimated inside the function using \link{estimate_locreg_cpp}.
 //'
 //' @examples
 //' \dontrun{
 //' library(Rcpp)
 //' library(data.table)
 //'
 //' # Sample data
 //'    set.seed(123)
 //'      n <- 100
 //'    data <- data.table(
 //'        id_curve = rep(1:10, each = n),
 //'        tobs = rep(seq(0, 1, length.out = n), times = 10),
 //'        X = rnorm(1000)
 //'    )
 //'
 //' # Observation points
 //'      t <- seq(0, 1, length.out = 50)
 //'
 //' # Estimate mean risk
 //'        mean_risk <- estimate_mean_risk_cpp(data, t)
 //'
 //' # Print the result
 //'          print(mean_risk)
 //' }
 //'
 //'
 //' @import Rdpack
 //'
 //' @references
 //' \insertAllCited{}
 //'
 // [[Rcpp::export]]
 arma::mat estimate_mean_risk_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                  const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                  const std::string kernel_name = "epanechnikov",
                                  const Rcpp::Nullable<arma::vec> presmooth_bw = R_NilValue,
                                  const Rcpp::Nullable<double> Delta = R_NilValue,
                                  const Rcpp::Nullable<arma::vec> presmooth_bw_grid = R_NilValue,
                                  const Rcpp::Nullable<int> presmooth_nsubset = R_NilValue) {
   if (t.size() == 0) {
     Rcpp::stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.");
   }
   core::Kern kern = core::parse_kernel(kernel_name);
   core::Curves C = core::make_curves(data);
   arma::vec bw_grid_to_use = bw_grid.isNull() ? default_mean_bw_grid(C) : Rcpp::as<arma::vec>(bw_grid);
   LocregCfg cfg = resolve_locreg(C, presmooth_bw, Delta, presmooth_bw_grid, presmooth_nsubset, kern);
   return mean_risk_impl(C, t, bw_grid_to_use, kern, cfg);
 }


 //' Estimate Mean Function for Functional Data
 //'
 //' This function estimates the mean function for functional data using the
 //' Nadaraya-Watson estimator with a specified kernel. The function takes a
 //' data frame of observations, a vector of evaluation points, and optional
 //' bandwidth parameters. See \insertCite{maissoro2024adaptive;textual}{adaptiveFTS}.
 //'
 //' @param data A data frame containing the functional data. It should have three columns:
 //' \code{id_curve} (curve identifiers), \code{tobs} (observation times), and \code{X} (observations).
 //' @param t A numeric vector of evaluation points between 0 and 1.
 //' @param bw An optional numeric vector of optimal bandwidths. If \code{NULL}, the function
 //' will estimate the optimal bandwidths.
 //' @param bw_grid An optional numeric vector of bandwidths to be used if \code{bw} is \code{NULL}.
 //' @param kernel_name A string specifying the kernel to be used. Supported kernels are:
 //' \code{"epanechnikov"}, \code{"biweight"}, \code{"triweight"}, \code{"tricube"},
 //' \code{"triangular"}, and \code{"uniform"}.
 //' @param presmooth_bw Numeric (positive vector or scalar). Bandwidth used to presmooth
 //' each curve in the local regularity step, see `estimate_locreg_cpp`. Default
 //' \code{NULL} selects it by cross-validation.
 //' @param Delta Numeric (positive). Length of the neighborhood of each point used in the
 //' local regularity step. Default \code{NULL} estimates it from the data.
 //' @param presmooth_bw_grid Numeric vector. Candidate bandwidths of the cross-validation
 //' that selects \code{presmooth_bw}. Default \code{NULL} uses the default grid.
 //' @param presmooth_nsubset Integer (positive). Number of curves used by that
 //' cross-validation. Default \code{NULL} uses min(70, floor(N / 2)) curves.
 //' @return A matrix where each row contains the following columns:
 //' \itemize{
 //'   \item \code{t} - The evaluation point.
 //'   \item \code{bw} - The optimal bandwidth used.
 //'   \item \code{Ht_used} - Intermediate result (Ht) used in the estimation.
 //'   \item \code{Lt_used} - Intermediate result (Lt) used in the estimation.
 //'   \item \code{PN} - The vector P_N(t;h).
 //'   \item \code{muhat} - The estimated mean function value at \code{t}.
 //' }
 //'
 //' @examples
 //' \dontrun{
 //' data <- data.frame(id_curve = rep(1:5, each = 10),
 //'                    tobs = rep(seq(0, 1, length.out = 10), 5),
 //'                    X = rnorm(50))
 //' t <- seq(0, 1, length.out = 100)
 //' estimate_mean_cpp(data, t)
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
 arma::mat estimate_mean_cpp(const Rcpp::DataFrame data, const arma::vec t,
                             const Rcpp::Nullable<arma::vec> bw = R_NilValue,
                             const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                             const std::string kernel_name = "epanechnikov",
                             const Rcpp::Nullable<arma::vec> presmooth_bw = R_NilValue,
                             const Rcpp::Nullable<double> Delta = R_NilValue,
                             const Rcpp::Nullable<arma::vec> presmooth_bw_grid = R_NilValue,
                             const Rcpp::Nullable<int> presmooth_nsubset = R_NilValue){
   if (t.size() == 0) {
     stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.");
   }
   core::Kern kern = core::parse_kernel(kernel_name);
   core::Curves C = core::make_curves(data);
   if (bw.isNull()) {
     arma::vec bw_grid_to_use = bw_grid.isNull() ? default_mean_bw_grid(C) : Rcpp::as<arma::vec>(bw_grid);
     LocregCfg cfg = resolve_locreg(C, presmooth_bw, Delta, presmooth_bw_grid, presmooth_nsubset, kern);
     return mean_impl(C, t, arma::vec(), bw_grid_to_use, kern, cfg);
   }
   arma::vec bw_vec = Rcpp::as<arma::vec>(bw);
   if ((int) bw_vec.n_elem != (int) t.n_elem)
     stop("If 'bw' is not NULL, it must be the same length as 't'.");
   return mean_impl(C, t, bw_vec, arma::vec(), kern, LocregCfg());
 }
