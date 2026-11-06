#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

#include "00_internal.h"
using namespace Rcpp;
using namespace arma;

// Risk of the covariance segment estimator, on the same structure as the mean
// risk: point constants once (moments, X(t)X(t) dependence, D numerators on the
// distinct points), per-curve quantities from the kernel windows.
arma::mat cov_segment_risk_impl(const core::Curves& C, const arma::vec& t, const arma::vec& bw_grid_to_use,
                                const bool center, core::Kern kern, const LocregCfg& cfg) {
  const int n = t.size();
  const int n_curve = C.N;
  const int bw_size = bw_grid_to_use.size();

  // Estimate local regularity
  arma::mat mat_locreg = locreg_impl(C, t, true, kern, cfg);
  arma::vec h(n_curve, arma::fill::value(mat_locreg(0, 1))); // the presmoothing bandwidth

  // Estimate the error sd
  arma::mat mat_sig = sigma_impl(C, t);

  // Point constants on the distinct points
  arma::vec U = arma::unique(t);
  arma::mat mat_mom_t = mom_impl(C, U, h, 2, center, kern);
  arma::vec var0, lrsum;
  xsxt_dependence(C, U, U, 0, arma::regspace(0, n_curve - 1), h, kern, center, var0, lrsum);
  arma::mat ddnum, ddPN;
  dd_matrices(C, U, bw_grid_to_use, h, 3, kern, center, ddnum, ddPN);

  std::vector<int> first = core::first_index(t);
  std::vector<int> pos(n);
  arma::vec Ht_k(n), Lt_k(n), sig2_k(n);
  for (int k = 0; k < n; ++k) {
    Ht_k(k) = mat_locreg(first[k], 4);
    Lt_k(k) = mat_locreg(first[k], 5);
    sig2_k(k) = mat_sig(first[k], 1) * mat_sig(first[k], 1);
    pos[k] = std::lower_bound(U.begin(), U.end(), t(k)) - U.begin();
  }

  arma::mat mat_res_risk(bw_size * n, 10);
#pragma omp parallel for schedule(dynamic)
  for (int k = 0; k < n; ++k) {
  const double Ht = Ht_k(k), Lt = Lt_k(k), sig_square = sig2_k(k);
  const double mom_t_square = mat_mom_t(pos[k], 2);
  std::vector<double> q_pn((std::size_t) bw_size * n_curve), q_wm((std::size_t) bw_size * n_curve),
    q_sb((std::size_t) bw_size * n_curve);
  for (int i = 0; i < n_curve; ++i)
    core::risk_stats_multi(C, i, t(k), bw_grid_to_use.memptr(), bw_size, Ht, kern,
                           &q_pn[i], &q_wm[i], &q_sb[i], n_curve);
  for (int idx_bw = 0; idx_bw < bw_size; ++idx_bw) {
    const double bw = bw_grid_to_use(idx_bw);
    const double bw_pow = std::pow(bw, 2 * Ht);
    const double *rp = &q_pn[(std::size_t) idx_bw * n_curve], *rw = &q_wm[(std::size_t) idx_bw * n_curve],
      *rs = &q_sb[(std::size_t) idx_bw * n_curve];

    double PN = 0, bias_term_num = 0, variance_term_num = 0;
    for (int i = 0; i < n_curve; ++i) {
      PN += rp[i];
      bias_term_num += Lt * bw_pow * rp[i] * rs[i];
      variance_term_num += sig_square * rp[i] * rw[i];
    }

    // Bias and variance terms
    double bias_term = 20 * mom_t_square * bias_term_num / PN;
    double variance_term = 20 * mom_t_square * variance_term_num / (PN * PN);

    // Dependence term
    double num_Dt = ddnum(pos[k], idx_bw);
    double PNt = ddPN(pos[k], idx_bw);
    double Dt = num_Dt / std::pow(PNt, 3);
    double first_dependence_term = 25 * Dt / PN;
    double second_dependence_term_num = var0(pos[k]) + lrsum(pos[k]) / PN;
    double second_dependence_term = second_dependence_term_num / PN;
    double dependence_term = first_dependence_term + second_dependence_term;

    double cov_segment_risk = bias_term + variance_term + dependence_term;

    double* M = mat_res_risk.memptr();
    const arma::uword nr = mat_res_risk.n_rows, ri = idx_bw * n + k;
    const double vals[10] = {t(k), bw, PN, h(0), Ht, Lt, bias_term, variance_term, dependence_term, cov_segment_risk};
    for (int col = 0; col < 10; ++col) M[ri + col * nr] = vals[col];
  }
  }
  return mat_res_risk;
}

arma::vec default_cov_segment_bw_grid(const core::Curves& C) {
  double n_curve = C.N;
  double lambdahat = C.lambdahat();
  double b0 = std::pow(n_curve * lambdahat, -1 / (2 * 0.05 + 1)); // rate with minimum local exponent = 0.05
  double bK = 0.2; // rate with maximum local exponent = 1
  return arma::logspace(log10(b0), log10(bK), 20);
}

//' Estimate the risk of the covariance segment function
 //'
 //' This function estimates the risk \eqn{R_{\Gamma_0}(t; h)} associated with the covariance segment line estimation proposed by \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
 //'
 //' @param data A DataFrame containing the columns "id_curve", "tobs", and "X", typically the output of the function \link{format_data}.
 //' @param t A numeric vector of observation points at which we want to estimate the mean function of the underlying process.
 //' @param bw_grid A numeric vector representing the bandwidth grid in which the best smoothing parameter is selected for each t.
 //' It can be NULL, in which case it will be defined as an exponential grid of \eqn{N \times \lambda}.
 //' @param center A logical value indicating if the data should be centered before estimation. Default is \code{true}.
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
 //'            \item{cov_segment_risk : The estimates of the risk function of the covariance segment function.}
 //'         }
 //'
 //' @details
 //' The local regularity parameters are directly estimated inside the function using \link{estimate_locreg_cpp}.
 //'
 //' The dependence term includes a term based on \eqn{\mathbb{D}(t; h_t)} derived from fourth-moment tensors, and another based on
 //' empirical autocovariance from \link{estimate_empirical_XsXt_autocov_cpp}.
 //'
 //' @examples
 //' \dontrun{
 //' library(Rcpp)
 //' library(data.table)
 //'
 //' # Simulated data
 //' set.seed(42)
 //' n <- 100
 //' data <- data.table(
 //'   id_curve = rep(1:10, each = n),
 //'   tobs = rep(seq(0, 1, length.out = n), times = 10),
 //'   X = rnorm(10 * n)
 //' )
 //'
 //' # Time points
 //' t <- seq(0.1, 0.9, length.out = 50)
 //'
 //' # Risk estimation
 //' res <- estimate_cov_segment_risk_cpp(data, t)
 //' head(res)
 //' }
 //'
 //' @import Rdpack
 //' @references
 //' \insertAllCited{}
 //'
 // [[Rcpp::export]]
 arma::mat estimate_cov_segment_risk_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                         const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                         const bool center = true,
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
   arma::vec bw_grid_to_use = bw_grid.isNull() ? default_cov_segment_bw_grid(C) : Rcpp::as<arma::vec>(bw_grid);
   LocregCfg cfg = resolve_locreg(C, presmooth_bw, Delta, presmooth_bw_grid, presmooth_nsubset, kern);
   return cov_segment_risk_impl(C, t, bw_grid_to_use, center, kern, cfg);
 }

 //' Estimate Covariance Segment Function for Functional Data
 //'
 //' Estimates the covariance segment function \eqn{\Gamma_{N,0}(t,t;h_t,h_t)} for functional data
 //' using the Nadaraya–Watson estimator with a specified kernel. This is part of the methodology
 //' described in \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
 //'
 //' @param data A data frame containing the functional data. It must have three columns:
 //' \code{id_curve} (curve identifiers), \code{tobs} (observation times), and \code{X} (observed values).
 //' @param t A numeric vector of evaluation points in the interval [0, 1].
 //' @param bw An optional numeric vector of bandwidths to be used at each evaluation point in \code{t}.
 //' If \code{NULL}, bandwidths are selected based on risk minimization.
 //' @param bw_grid An optional numeric vector of candidate bandwidths used for optimization when \code{bw} is \code{NULL}.
 //' @param center A logical value indicating if the data should be centered before estimation. Default is \code{true}.
 //' @param kernel_name A string indicating the kernel to use. Supported kernels are:
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
 //'
 //' @return A matrix with 8 columns and \code{length(t)} rows:
 //' \describe{
 //'   \item{\code{t}}{Evaluation point}
 //'   \item{\code{bw}}{Bandwidth used at \code{t}}
 //'   \item{\code{Ht_used}}{Intermediate quantity \eqn{H(t)} used in estimation}
 //'   \item{\code{Lt_used}}{Intermediate quantity \eqn{L(t)} used in estimation}
 //'   \item{\code{PN}}{Number of curves contributing to the estimation at \code{t}}
 //'   \item{\code{cov_segment_hat}}{Uncorrected covariance segment estimate}
 //'   \item{\code{covseg_correction}}{Correction term based on measurement error variance}
 //'   \item{\code{cov_segment_hat_correct}}{Final corrected covariance segment estimate}
 //' }
 //'
 //' @examples
 //' \dontrun{
 //' data <- data.frame(
 //'   id_curve = rep(1:5, each = 10),
 //'   tobs = rep(seq(0, 1, length.out = 10), 5),
 //'   X = rnorm(50)
 //' )
 //' t <- seq(0, 1, length.out = 100)
 //' estimate_cov_segment_cpp(data, t)
 //' }
 //'
 //' @import Rdpack
 //' @references
 //' \insertAllCited{}
 //'
 //'
 // [[Rcpp::export]]
 arma::mat estimate_cov_segment_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                    const Rcpp::Nullable<arma::vec> bw = R_NilValue,
                                    const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                    const bool center = true,
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
   int n = t.size();
   int n_curve = C.N;

   // Estimate the error sd
   arma::vec sig_square_vec = arma::square(sigma_impl(C, t).col(1));

   // Estimate / Verify optimal risk function
   arma::vec optbw_to_use(n);
   arma::vec Ht_used(n, arma::fill::zeros);
   arma::vec Lt_used(n, arma::fill::zeros);
   if (bw.isNull()) {
     arma::vec bw_grid_to_use = bw_grid.isNull() ? default_cov_segment_bw_grid(C) : Rcpp::as<arma::vec>(bw_grid);
     LocregCfg cfg = resolve_locreg(C, presmooth_bw, Delta, presmooth_bw_grid, presmooth_nsubset, kern);
     arma::mat mat_risk = cov_segment_risk_impl(C, t, bw_grid_to_use, center, kern, cfg);
     const int bw_size = bw_grid_to_use.n_elem;
     std::vector<std::vector<arma::uword>> groups = core::same_value_groups(t);
     for (int k = 0; k < n; ++k) {
       // As the former selection: the position of the first minimum among the
       // finite risks, applied to the list of all the rows of t(k).
       std::vector<arma::uword> rows;
       double best = 0;
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
       if (n_finite == 0) stop("index_min(): object has no elements");
       arma::uword best_row = rows[best_rank];
       optbw_to_use(k) = mat_risk(best_row, 1);
       Ht_used(k) = mat_risk(best_row, 4);
       Lt_used(k) = mat_risk(best_row, 5);
     }
   } else {
     arma::vec optbw_cur = as<arma::vec>(bw);
     int optbw_cur_size = optbw_cur.size();
     if (optbw_cur_size != n) {
       stop("If 'bw' is not NULL, it must be the same length as 't'.");
     } else {
       optbw_to_use = optbw_cur;
     }
   }

   // Compute the mean function
   arma::vec muhat = mean_impl(C, t, optbw_to_use, arma::vec(), kern, LocregCfg()).col(5);

   // Compute the covariance segment function
   arma::mat mat_res_covseg(n, 8);
   arma::vec PN(n, fill::zeros);
   arma::vec covseg_numerator(n, fill::zeros);
   arma::vec covseg_correction_numerator(n, fill::zeros);

   for (int i = 0 ; i < n_curve; ++i) {
     const std::vector<double>& ts = C.Ts[i];
     for (int k = 0; k < n; ++k) {
       const double x = t(k), b = optbw_to_use(k);
       // Smooth (NaN -> 0), support indicator, and sum of the squared weights
       double Xhat = core::nw_at(C, i, x, b, kern);
       if (!std::isfinite(Xhat)) Xhat = 0;
       double pn = core::pn_abs(C, i, x, b);
       arma::uword lo, hi;
       core::window(C, i, x, b, lo, hi);
       double S = 0, S2 = 0;
       for (arma::uword q = lo; q < hi; ++q) {
         double w = core::kval(kern, (ts[q] - x) / b);
         if (!std::isfinite(w)) w = 0;
         S += w;
         S2 += w * w;
       }
       double sum_weight_square = (S != 0 && std::isfinite(S)) ? S2 / (S * S) : 0.0;
       PN(k) += pn;
       covseg_numerator(k) += pn * (Xhat - muhat(k)) * (Xhat - muhat(k));
       covseg_correction_numerator(k) += pn * sum_weight_square;
     }
   }

   // Compute \widehat \Gamma_{N,0}(t,t;h_t, h_t)
   arma::vec cov_segment_hat = covseg_numerator / PN;
   arma::vec covseg_correction = (sig_square_vec % covseg_correction_numerator) / PN;
   arma::vec cov_segment_hat_correct = cov_segment_hat - covseg_correction;

   // return the result
   mat_res_covseg.col(0) = t;
   mat_res_covseg.col(1) = optbw_to_use;
   mat_res_covseg.col(2) = Ht_used;
   mat_res_covseg.col(3) = Lt_used;
   mat_res_covseg.col(4) = PN;
   mat_res_covseg.col(5) = cov_segment_hat;
   mat_res_covseg.col(6) = covseg_correction;
   mat_res_covseg.col(7) = cov_segment_hat_correct;
   return mat_res_covseg;
 }
