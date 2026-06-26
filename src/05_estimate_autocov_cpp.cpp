#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <map>
#include <utility>

// [[Rcpp::depends(RcppArmadillo)]]

#include "02_smoothing_rcpp.h"
#include "03_estimate_locreg_rcpp.h"
#include "07_estimate_constants_cpp.h"
#include "04_estimate_mean_cpp.h"
using namespace Rcpp;
using namespace arma;

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
 //' @param use_same_bw A logical value indicating if the same bandwidth should be used for \code{s} and \code{t}. Default is \code{false}.
 //' @param center A logical value indicating if the data should be centered before estimation. Default is \code{true}.
 //' @param kernel_name \code{string}. Specifies the kernel function for estimation; default is "epanechnikov".
 //' Supported kernels include: "epanechnikov", "biweight", "triweight", "tricube", "triangular", and "uniform".
 //'
 //' @return A \code{matrix} containing the following fourteen columns in order:
 //'          \itemize{
 //'            \item{s : The first argument of the autocovariance function.}
 //'            \item{t : The second argument of the autocovariance function.}
 //'            \item{hs : The candidate bandwidth for the first argument of the autocovariance function. If \code{use_same_bw = TRUE}, the same bandwidth candidate is used for \code{s} and for \code{t}, so the 3rd and 4th columns contain the same values.}
 //'            \item{ht : The candidate bandwidth for the second argument of the autocovariance function.}
 //'            \item{PNl : The number of curves used in the estimation the autocovariance at (s,t). It corresponds to \eqn{P_{N,\ell}(s,t;h_s, h_t)}.}
 //'            \item{locreg_bw : The bandwidth used to estimate the local regularity parameters.}
 //'            \item{Hs : The estimates of the local exponent for each t. It corresponds to \eqn{H_s}.}
 //'            \item{Ls : The estimates of the Hölder constant for each t. It corresponds to \eqn{L_s^2}.}
 //'            \item{Ht : The estimates of the local exponent for each t. It corresponds to \eqn{H_t}.}
 //'            \item{Lt : The estimates of the Hölder constant for each t. It corresponds to \eqn{L_t^2}.}
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
                                     const bool use_same_bw = false,
                                     const bool center = true,
                                     const std::string kernel_name = "epanechnikov"){
   // NB : We consider that no vector is null
   // NB : if use_same_bw = TRUE, then the same bandwidth is used for s and t
   //      otherwise, a different bandwidth is use for s and for t.

   if (s.size() != t.size()) {
     stop("Arguments 's' and 't' must be of equal length.");
   } else if (s.size() == 0) {
     stop("Arguments 's' and 't' must be a numeric vectors or scalar value(s) between 0 and 1.");
   }
   // Check if kernel_name is one of the supported kernels
   std::vector<std::string> supported_kernels = {"epanechnikov", "biweight", "triweight",
                                                 "tricube", "triangular", "uniform"};
   if (std::find(supported_kernels.begin(), supported_kernels.end(), kernel_name) == supported_kernels.end()) {
     stop("Unsupported kernel name. Choose from: epanechnikov, biweight, triweight, tricube, triangular, uniform.");
   }

   // Select the kernel function
   std::function<arma::vec(const arma::vec)> kernel_func = select_kernel(kernel_name);

   // Take the number of couples (s,t)
   int n = t.size();

   // Convert data to matrix
   arma::mat data_mat(data.nrows(), 3);
   data_mat.col(0) = as<arma::vec>(data["id_curve"]);
   data_mat.col(1) = as<arma::vec>(data["tobs"]);
   data_mat.col(2) = as<arma::vec>(data["X"]);
   arma::vec unique_id_curve = arma::sort(arma::unique(data_mat.col(0)), "ascend");
   double n_curve = unique_id_curve.n_elem;

   // Set the bandwidth grid
   arma::vec bw_grid_to_use;
   if (bw_grid.isNull()) {
     // Estimate lambda
     double lambdahat = arma::mean(hist(data_mat.col(0), unique_id_curve));
     double b0 = 5 * std::max(std::pow(n_curve * lambdahat, -1 / (2 * 0.1 + 1)),
                              std::pow(n_curve * lambdahat * lambdahat, -1 / (2 * 0.1 + 1))); // rate with minimum local exponent = 0.05
     // double bK = 4 * std::max(std::pow(n_curve * lambdahat, -1 / (2 * 0.8 + 1)),
     //                          std::pow(n_curve * lambdahat * lambdahat, -1 / (2 * 0.8 + 1))); // rate with maximum local exponent = 0.8
     double bK = 0.2;
     bw_grid_to_use = arma::logspace(log10(b0), log10(bK), 15);
   } else {
     bw_grid_to_use = as<arma::vec>(bw_grid);
   }
   int bw_size = bw_grid_to_use.size();

   // Precompute the row indices of each curve once. These index sets are
   // constant, but were previously recomputed with arma::find() inside the
   // innermost loop (once per curve, per (s,t) pair, per bandwidth pair),
   // which dominated the runtime. Caching them changes nothing numerically.
   std::vector<arma::uvec> curve_idx((arma::uword) n_curve);
   for (arma::uword c = 0; c < (arma::uword) n_curve; ++c)
     curve_idx[c] = arma::find(data_mat.col(0) == unique_id_curve(c));

   // Estimate local regularity
   arma::mat mat_locreg_s = estimate_locreg_cpp(data, arma::unique(s), true, kernel_name, R_NilValue, R_NilValue);
   arma::mat mat_locreg_t = estimate_locreg_cpp(data, arma::unique(t), true, kernel_name, R_NilValue, R_NilValue);
   arma::vec h(n_curve, fill::value(mat_locreg_s(0, 1))); // extract the presmoothing bandwidth

   // Estimate the error sd
   arma::mat mat_sig_s = estimate_sigma_cpp(data, arma::unique(s));
   arma::mat mat_sig_t = estimate_sigma_cpp(data, arma::unique(t));
   arma::vec sig2_vec_s = arma::pow(mat_sig_s.col(1), 2);
   arma::vec sig2_vec_t = arma::pow(mat_sig_t.col(1), 2);

   // Estimate moment
   arma::mat mat_mom_s = estimate_empirical_mom_cpp(data, arma::unique(s), h, 2, center, kernel_name);
   arma::mat mat_mom_t = estimate_empirical_mom_cpp(data, arma::unique(t), h, 2, center, kernel_name);
   arma::vec mom_vec_s = mat_mom_s.col(2);
   arma::vec mom_vec_t = mat_mom_t.col(2);

   // Estimate the autocovariance
   arma::mat mat_emp_autocov = estimate_empirical_XsXt_autocov_cpp(data, s, t, lag, arma::regspace(0, n_curve - lag - 1), h, kernel_name, center);

   // Compute the numerator of \mathbb{D}(t;h_t), \mathbb{D}(s;h_s), \mathbb{D}_\ell(t;h_t|s,h_s) \mathbb{D}_\ell(s;h_s|t,h_t)
   arma::mat mat_num_DD_s = estimate_numerator_dependence_term_DD_cpp(data, arma::unique(s), bw_grid_to_use, h, 3, kernel_name, center);
   arma::mat mat_num_DD_t = estimate_numerator_dependence_term_DD_cpp(data, arma::unique(t), bw_grid_to_use, h, 3, kernel_name, center);

   // Definie result matrix
   int n_rows_res = use_same_bw ? bw_size * n : bw_size * bw_size * n;
   arma::mat mat_res_risk(n_rows_res, 14);

   // ---- Precompute quantities shared across the bandwidth loops ----
   // Previously every constant below was recomputed inside the O(bw^2) x n loop.
   const int nc = (int) n_curve;

   // Per-pair (s(k), t(k)) constants: local regularity, second moments, error
   // variances, and the parts of the second dependence term independent of PN_lag.
   arma::vec Hs_k(n), Ls_k(n), Ht_k(n), Lt_k(n);
   arma::vec moms_k(n), momt_k(n), sigs_k(n), sigt_k(n);
   arma::vec XsXtvar_k(n), sumabslr_k(n);
   for (int k = 0; k < n; ++k) {
     arma::uvec is = arma::find(mat_locreg_s.col(0) == s(k));
     arma::uvec it = arma::find(mat_locreg_t.col(0) == t(k));
     Hs_k(k) = mat_locreg_s(is(0), 4); Ls_k(k) = mat_locreg_s(is(0), 5);
     Ht_k(k) = mat_locreg_t(it(0), 4); Lt_k(k) = mat_locreg_t(it(0), 5);

     arma::uvec ims = arma::find(mat_mom_s.col(0) == s(k));
     arma::uvec imt = arma::find(mat_mom_t.col(0) == t(k));
     moms_k(k) = mom_vec_s(ims(0)); momt_k(k) = mom_vec_t(imt(0));

     arma::uvec iss = arma::find(mat_sig_s.col(0) == s(k));
     arma::uvec ist = arma::find(mat_sig_t.col(0) == t(k));
     sigs_k(k) = sig2_vec_s(iss(0)); sigt_k(k) = sig2_vec_t(ist(0));

     arma::uvec idx_lag0 = arma::find( (mat_emp_autocov.col(0) == s(k)) % (mat_emp_autocov.col(1) == t(k)) % (mat_emp_autocov.col(3) == 0) );
     arma::uvec idx_lag  = arma::find( (mat_emp_autocov.col(0) == s(k)) % (mat_emp_autocov.col(1) == t(k)) % (mat_emp_autocov.col(3) != 0) );
     XsXtvar_k(k) = mat_emp_autocov(idx_lag0(0), 5);
     arma::mat XsXt_mat_lr_var = mat_emp_autocov.rows(idx_lag);
     sumabslr_k(k) = arma::accu(arma::abs(2 * XsXt_mat_lr_var.col(5)));
   }

   // Per-pair, per-bandwidth dependence-term numerators and counts.
   arma::mat numDs_kb(n, bw_size), PNs_kb(n, bw_size);
   arma::mat numDt_kb(n, bw_size), PNt_kb(n, bw_size);
   for (int k = 0; k < n; ++k) {
     for (int b = 0; b < bw_size; ++b) {
       double bw = bw_grid_to_use[b];
       arma::uvec ids = arma::find( (mat_num_DD_s.col(0) == s(k)) % (mat_num_DD_s.col(1) == bw) );
       arma::uvec idt = arma::find( (mat_num_DD_t.col(0) == t(k)) % (mat_num_DD_t.col(1) == bw) );
       numDs_kb(k, b) = mat_num_DD_s(ids(0), 2); PNs_kb(k, b) = mat_num_DD_s(ids(0), 3);
       numDt_kb(k, b) = mat_num_DD_t(idt(0), 2); PNt_kb(k, b) = mat_num_DD_t(idt(0), 3);
     }
   }

   // Per-curve kernel-weight aggregates. These depend only on the point value
   // and bandwidth, and s/t typically take few distinct values (an m x m grid
   // has only m distinct s and m distinct t), so compute each once per distinct
   // value and expand via ks_idx / kt_idx. Bit-identical; removes the previous
   // per-pair redundancy.
   //   *_pn  : 1 if the curve has an observation in the kernel support, else 0
   //   *_sbn : sum |(T - x)/bw|^{2H} * |normalised weight|   (bias numerator)
   //   *_max : max normalised weight                          (variance numerator)
   std::map<double, int> us_map, ut_map;
   arma::uvec ks_idx(n), kt_idx(n);
   std::vector<double> sval_u_v, tval_u_v;
   for (int k = 0; k < n; ++k) {
     auto is2 = us_map.find(s(k));
     if (is2 == us_map.end()) { int id = (int) sval_u_v.size(); us_map.emplace(s(k), id); sval_u_v.push_back(s(k)); ks_idx(k) = id; }
     else ks_idx(k) = is2->second;
     auto it2 = ut_map.find(t(k));
     if (it2 == ut_map.end()) { int id = (int) tval_u_v.size(); ut_map.emplace(t(k), id); tval_u_v.push_back(t(k)); kt_idx(k) = id; }
     else kt_idx(k) = it2->second;
   }
   arma::vec sval_u(sval_u_v), tval_u(tval_u_v);
   int Us = sval_u.n_elem, Ut = tval_u.n_elem;
   arma::vec Hs_u(Us), Ht_u(Ut);
   for (int j = 0; j < Us; ++j) { arma::uvec ii = arma::find(mat_locreg_s.col(0) == sval_u(j)); Hs_u(j) = mat_locreg_s(ii(0), 4); }
   for (int j = 0; j < Ut; ++j) { arma::uvec ii = arma::find(mat_locreg_t.col(0) == tval_u(j)); Ht_u(j) = mat_locreg_t(ii(0), 4); }

   arma::vec tobs_all = data_mat.col(1);
   // Cache each curve's observation points once (avoids re-gathering per bw).
   std::vector<arma::vec> curve_tobs((arma::uword) nc);
   for (int c = 0; c < nc; ++c) curve_tobs[c] = tobs_all(curve_idx[c]);

   arma::cube s_pn(nc, Us, bw_size), s_sbn(nc, Us, bw_size), s_max(nc, Us, bw_size);
   arma::cube t_pn(nc, Ut, bw_size), t_sbn(nc, Ut, bw_size), t_max(nc, Ut, bw_size);
#pragma omp parallel for
   for (int j = 0; j < Us; ++j) {
     double Hs = Hs_u(j);
     for (int b = 0; b < bw_size; ++b) {
       double bw = bw_grid_to_use[b];
       for (int c = 0; c < nc; ++c) {
         const arma::vec& Tc = curve_tobs[c];
         arma::vec ds = (Tc - sval_u(j)) / bw;
         arma::vec ws = kernel_func(ds);
         ws.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
         arma::vec wsn = ws / arma::accu(ws);
         wsn.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
         s_pn(c, j, b)  = arma::find(arma::abs(ds) <= 1).is_empty() ? 0.0 : 1.0;
         s_sbn(c, j, b) = arma::sum(arma::pow(arma::abs(ds), 2 * Hs) % arma::abs(wsn));
         s_max(c, j, b) = wsn.max();
       }
     }
   }
#pragma omp parallel for
   for (int j = 0; j < Ut; ++j) {
     double Ht = Ht_u(j);
     for (int b = 0; b < bw_size; ++b) {
       double bw = bw_grid_to_use[b];
       for (int c = 0; c < nc; ++c) {
         const arma::vec& Tc = curve_tobs[c];
         arma::vec dt = (Tc - tval_u(j)) / bw;
         arma::vec wt = kernel_func(dt);
         wt.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
         arma::vec wtn = wt / arma::accu(wt);
         wtn.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
         t_pn(c, j, b)  = arma::find(arma::abs(dt) <= 1).is_empty() ? 0.0 : 1.0;
         t_sbn(c, j, b) = arma::sum(arma::pow(arma::abs(dt), 2 * Ht) % arma::abs(wtn));
         t_max(c, j, b) = wtn.max();
       }
     }
   }

   // Compute the risk for each (s,t) and each bandwidth in bw_grid
   if (use_same_bw) {
#pragma omp parallel for
     for (int idx_bw = 0; idx_bw < bw_size; ++idx_bw) {
       double bw = bw_grid_to_use[idx_bw];
       for (int k = 0; k < n; ++k) {
         double Hs = Hs_k(k), Ls = Ls_k(k), Ht = Ht_k(k), Lt = Lt_k(k);
         double mom_s_square = moms_k(k), mom_t_square = momt_k(k);
         double sig_s_square = sigs_k(k), sig_t_square = sigt_k(k);
         double bw_s_pow = std::pow(bw, 2 * Hs);
         double bw_t_pow = std::pow(bw, 2 * Ht);

         const int ks = (int) ks_idx(k), kt = (int) kt_idx(k);
         double bias_term_num = 0, variance_term_num = 0, PN_lag = 0;
         for (int i = 0; i < n_curve - lag; ++i) {
           double pns = s_pn(i, ks, idx_bw);
           double pnt = t_pn(i + lag, kt, idx_bw);
           PN_lag += pns * pnt;
           bias_term_num += mom_t_square * Ls * bw_s_pow * pns * pnt * s_sbn(i, ks, idx_bw) +
             mom_s_square * Lt * bw_t_pow * pns * pnt * t_sbn(i + lag, kt, idx_bw);
           variance_term_num += sig_s_square * mom_t_square * pns * pnt * s_max(i, ks, idx_bw) +
             sig_t_square * mom_s_square * pns * pnt * t_max(i + lag, kt, idx_bw) +
             sig_s_square * sig_t_square * pns * pnt * s_max(i, ks, idx_bw) * t_max(i + lag, kt, idx_bw);
         }

         // Compute bias and variance terms
         double bias_term = 4 * bias_term_num / PN_lag;
         double variance_term = 4 * variance_term_num / (PN_lag * PN_lag);

         // Dependence term (numerators / counts precomputed)
         double num_Ds = numDs_kb(k, idx_bw), num_Dt = numDt_kb(k, idx_bw);
         double PNs = PNs_kb(k, idx_bw), PNt = PNt_kb(k, idx_bw);
         double Ds = num_Ds / std::pow(PNs, 3) ;
         double Dt = num_Dt / std::pow(PNt, 3) ;
         double Ds_t = num_Ds / std::pow(PN_lag, 3) ;
         double Dt_s = num_Dt / std::pow(PN_lag, 3) ;
         double first_dependence_term = 15 * ( sqrt(Ds_t * Dt) + sqrt(Dt_s * Ds) + 3 * sqrt(Ds * Dt) ) / PN_lag;
         double second_dependence_term_num = XsXtvar_k(k) + sumabslr_k(k) / PN_lag;
         double second_dependence_term = (second_dependence_term_num / PN_lag);
         double dependence_term = first_dependence_term + second_dependence_term;

         double autocov_risk = bias_term + variance_term + dependence_term;
         mat_res_risk.row(idx_bw * n + k) = {s(k), t(k), bw, bw, PN_lag, h(0), Hs, Ls, Ht, Lt, bias_term, variance_term, dependence_term, autocov_risk};
       }
     }

   } else {
     // If two bandwidth are used
#pragma omp parallel for
     for (int idx_bw_s = 0; idx_bw_s < bw_size; ++idx_bw_s) {
       double bw_s = bw_grid_to_use[idx_bw_s];
       for (int idx_bw_t = 0; idx_bw_t < bw_size; ++idx_bw_t) {
         double bw_t = bw_grid_to_use[idx_bw_t];
         for (int k = 0; k < n; ++k) {
           double Hs = Hs_k(k), Ls = Ls_k(k), Ht = Ht_k(k), Lt = Lt_k(k);
           double mom_s_square = moms_k(k), mom_t_square = momt_k(k);
           double sig_s_square = sigs_k(k), sig_t_square = sigt_k(k);
           double bw_s_pow = std::pow(bw_s, 2 * Hs);
           double bw_t_pow = std::pow(bw_t, 2 * Ht);

           const int ks = (int) ks_idx(k), kt = (int) kt_idx(k);
           double bias_term_num = 0, variance_term_num = 0, PN_lag = 0;
           for (int i = 0; i < n_curve - lag; ++i) {
             double pns = s_pn(i, ks, idx_bw_s);
             double pnt = t_pn(i + lag, kt, idx_bw_t);
             PN_lag += pns * pnt;
             bias_term_num += mom_t_square * Ls * bw_s_pow * pns * pnt * s_sbn(i, ks, idx_bw_s) +
               mom_s_square * Lt * bw_t_pow * pns * pnt * t_sbn(i + lag, kt, idx_bw_t);
             // NB: this branch historically squares pn_t (not pn_s * pn_t) in the
             // variance term; preserved exactly to keep results bit-identical.
             variance_term_num += sig_s_square * mom_t_square * pnt * pnt * s_max(i, ks, idx_bw_s) +
               sig_t_square * mom_s_square * pnt * pnt * t_max(i + lag, kt, idx_bw_t) +
               sig_s_square * sig_t_square * pnt * pnt * s_max(i, ks, idx_bw_s) * t_max(i + lag, kt, idx_bw_t);
           }

           double bias_term = 4 * bias_term_num / PN_lag;
           double variance_term = 4 * variance_term_num / (PN_lag * PN_lag);

           double num_Ds = numDs_kb(k, idx_bw_s), num_Dt = numDt_kb(k, idx_bw_t);
           double PNs = PNs_kb(k, idx_bw_s), PNt = PNt_kb(k, idx_bw_t);
           double Ds = num_Ds / std::pow(PNs, 3) ;
           double Dt = num_Dt / std::pow(PNt, 3) ;
           double Ds_t = num_Ds / std::pow(PN_lag, 3) ;
           double Dt_s = num_Dt / std::pow(PN_lag, 3) ;
           double first_dependence_term = 15 * ( sqrt(Ds_t * Dt) + sqrt(Dt_s * Ds) + 3 * sqrt(Ds * Dt) ) / PN_lag;
           double second_dependence_term_num = XsXtvar_k(k) + sumabslr_k(k) / PN_lag;
           double second_dependence_term = (second_dependence_term_num / PN_lag);
           double dependence_term = first_dependence_term + second_dependence_term;

           double autocov_risk = bias_term + variance_term + dependence_term;
           mat_res_risk.row(idx_bw_s * bw_size * n + idx_bw_t * n + k) = {s(k), t(k), bw_s, bw_t, PN_lag, h(0), Hs, Ls, Ht, Lt, bias_term, variance_term, dependence_term, autocov_risk};
         }
       }
     }
   }

   return mat_res_risk;
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

   // Remove duplicated rows
   std::unordered_set<std::string> seen;
   std::vector<arma::uword> unique_indices;

   for (arma::uword i = 0; i < mat_st.n_rows; ++i) {
     std::string row_str = std::to_string(mat_st(i, 0)) + "," + std::to_string(mat_st(i, 1));
     if (seen.find(row_str) == seen.end()) {
       seen.insert(row_str);
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

 //' Estimate the autocovariance function
 //'
 //' This function estimates the lag-\eqn{\ell}, \eqn{\ell} = 0, 1,..., autocovariance function based on the methodology of
 //' \insertCite{maissoro2024adaptive;textual}{adaptiveFTS} and \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
 //'
 //' @param data A \code{DataFrame} containing the data with columns \code{"id_curve"}, \code{"tobs"}, and \code{"X"}.
 //' @param s A numeric vector specifying time points \code{s} for which to estimate autocovariance.
 //' @param t A numeric vector specifying time points \code{t} for which to estimate autocovariance.
 //' @param lag An integer specifying the lag value for autocovariance.
 //' @param optbw_s Optional numeric vector specifying optimal bandwidths for \code{s}. Default is \code{NULL}.
 //' @param optbw_t Optional numeric vector specifying optimal bandwidths for \code{t}. Default is \code{NULL}.
 //' @param bw_grid numeric vector of bandwidth grid values. Default is \code{NULL}.
 //' @param use_same_bw A logical value indicating if the same bandwidth should be used for \code{s} and \code{t}. Default is \code{false}.
 //' @param center A logical value indicating if the data should be centered before estimation. Default is \code{true}.
 //' @param correct_diagonal A logical value indicating whether the diagonal of the covariance should be corrected when \code{lag=0}.
 //' @param kernel_name A string specifying the kernel to use for estimation. Supported values are \code{"epanechnikov"}, \code{"biweight"}, \code{"triweight"}, \code{"tricube"}, \code{"triangular"}, \code{"uniform"}. Default is \code{"epanechnikov"}.
 //'
 //' @return A \code{matrix} containing the following fourteen columns in order:
 //'          \itemize{
 //'            \item{s : The first argument of the autocovariance function.}
 //'            \item{t : The second argument of the autocovariance function.}
 //'            \item{optbw_s : The optimal bandwidth for the first argument of the autocovariance function. If \code{use_same_bw = TRUE}, the same bandwidth candidate is used for \code{s} and for \code{t}, so the 3rd and 4th columns contain the same values.}
 //'            \item{optbw_t : The optimal bandwidth for the second argument of the autocovariance function.}
 //'            \item{Hs : The estimates of the local exponent for each \code{s}. It corresponds to \eqn{H_s}.}
 //'            \item{Ls : The estimates of the Hölder constant for each \code{s}. It corresponds to \eqn{L_s^2}.}
 //'            \item{Ht : The estimates of the local exponent for each \code{t}. It corresponds to \eqn{H_t}.}
 //'            \item{Lt : The estimates of the Hölder constant for each \code{t}. It corresponds to \eqn{L_t^2}.}
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
                                const Rcpp::Nullable<arma::vec> optbw_s = R_NilValue,
                                const Rcpp::Nullable<arma::vec> optbw_t = R_NilValue,
                                const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                const bool use_same_bw = false,
                                const bool center = true,
                                const bool correct_diagonal = true,
                                const std::string kernel_name = "epanechnikov") {
   if (s.size() != t.size()) {
     stop("Arguments 's' and 't' must be of equal length.");
   } else if (s.size() == 0) {
     stop("Arguments 's' and 't' must be numeric vectors or scalar value(s) between 0 and 1.");
   }

   // Check if kernel_name is one of the supported kernels
   std::vector<std::string> supported_kernels = {"epanechnikov", "biweight", "triweight",
                                                 "tricube", "triangular", "uniform"};
   if (std::find(supported_kernels.begin(), supported_kernels.end(), kernel_name) == supported_kernels.end()) {
     stop("Unsupported kernel name. Choose from: epanechnikov, biweight, triweight, tricube, triangular, uniform.");
   }

   // Select the kernel function
   std::function<arma::vec(const arma::vec)> kernel_func = select_kernel(kernel_name);

   // Convert data to matrix
   arma::mat data_mat(data.nrows(), 3);
   data_mat.col(0) = as<arma::vec>(data["id_curve"]);
   data_mat.col(1) = as<arma::vec>(data["tobs"]);
   data_mat.col(2) = as<arma::vec>(data["X"]);
   arma::vec unique_id_curve = arma::sort(arma::unique(data_mat.col(0)), "ascend");
   double n_curve = unique_id_curve.n_elem;

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

   if (optbw_s.isNull() || optbw_t.isNull()) {
     arma::mat mat_risk = estimate_autocov_risk_cpp(data, svec, tvec, lag, bw_grid, use_same_bw, center, kernel_name);
     for (int k = 0; k < n; ++k) {
       arma::uvec idx_risk_cur = arma::find((mat_risk.col(0) == svec(k)) % (mat_risk.col(1) == tvec(k)));
       arma::vec risk = mat_risk(idx_risk_cur, arma::uvec({13}));
       arma::uword idx_min = arma::index_min(risk.elem(arma::find_finite(risk)));

       optbw_s_to_use(k) = mat_risk(idx_risk_cur(idx_min), 2);
       optbw_t_to_use(k) = mat_risk(idx_risk_cur(idx_min), 3);
       mat_locreg.row(k) = {mat_risk(idx_risk_cur(idx_min), 0), mat_risk(idx_risk_cur(idx_min), 1),
                            mat_risk(idx_risk_cur(idx_min), 6), mat_risk(idx_risk_cur(idx_min), 7),
                            mat_risk(idx_risk_cur(idx_min), 8), mat_risk(idx_risk_cur(idx_min), 9)};
     }
   } else {
     arma::vec optbw_s_to_use_temp = as<arma::vec>(optbw_s);
     arma::vec optbw_t_to_use_temp = as<arma::vec>(optbw_t);
     int optbw_s_to_use_temp_size = optbw_s_to_use_temp.size();
     int optbw_t_to_use_temp_size = optbw_t_to_use_temp.size();

     if ((optbw_s_to_use_temp_size == n) && (optbw_t_to_use_temp_size == n)) {
       optbw_s_to_use = optbw_s_to_use_temp;
       optbw_t_to_use = optbw_t_to_use_temp;
     } else if ((optbw_s_to_use_temp_size == 1) && (optbw_t_to_use_temp_size == 1)){
       optbw_s_to_use = arma::ones(n) * (optbw_s_to_use_temp(0));
       optbw_t_to_use = arma::ones(n) * (optbw_t_to_use_temp(0));
     } else if (lag == 0) {
       arma::vec optbw_svec(n);
       arma::vec optbw_tvec(n);
       for (int j = 0; j < n; ++j) {
         arma::uvec idx_optbw_svecj = arma::find(s == svec(j));
         if (idx_optbw_svecj.is_empty()) {
           // In case of svec(j) corresponds to a point in t
           idx_optbw_svecj = arma::find(t == svec(j));
         }
         arma::uvec idx_optbw_tvecj = arma::find(t == tvec(j));
         if (idx_optbw_tvecj.is_empty()) {
           // In case of tvec(j) corresponds to a point in s
           idx_optbw_tvecj = arma::find(s == tvec(j));
         }
         optbw_svec(j) = optbw_s_to_use_temp(idx_optbw_svecj(0));
         optbw_tvec(j) = optbw_t_to_use_temp(idx_optbw_tvecj(0));
       }
       optbw_s_to_use = optbw_svec;
       optbw_t_to_use = optbw_tvec;

     } else {
       stop("If 'optbw_s' and 'optbw_t' are not NULL, they must be the same length as 's' and as 't' or of length 1.");
     }

     // init log reg mat
     mat_locreg.col(0) = svec;
     mat_locreg.col(1) = tvec;
     mat_locreg.cols(2, 5).zeros();
   }

   // Deduplicate (point, bandwidth) tuples. The per-pair NW smoothing, support
   // indicators, diagonal weights and mean estimates depend only on a single
   // (s, h_s) or (t, h_t) tuple, yet svec/tvec typically repeat the same values
   // many times (e.g. an m x m covariance grid has only m distinct s and t).
   // Compute each distinct tuple once and expand back via map_s / map_t. This
   // is numerically identical: each estimate is independent of the others.
   std::map<std::pair<double, double>, int> umap_s, umap_t;
   arma::uvec map_s(n), map_t(n);
   std::vector<double> svec_u_v, obs_u_v, tvec_u_v, obt_u_v;
   for (int k = 0; k < n; ++k) {
     std::pair<double, double> ks(svec(k), optbw_s_to_use(k));
     auto it = umap_s.find(ks);
     if (it == umap_s.end()) {
       int id = (int) svec_u_v.size();
       umap_s.emplace(ks, id); svec_u_v.push_back(svec(k)); obs_u_v.push_back(optbw_s_to_use(k));
       map_s(k) = id;
     } else { map_s(k) = it->second; }
     std::pair<double, double> kt(tvec(k), optbw_t_to_use(k));
     auto jt = umap_t.find(kt);
     if (jt == umap_t.end()) {
       int id = (int) tvec_u_v.size();
       umap_t.emplace(kt, id); tvec_u_v.push_back(tvec(k)); obt_u_v.push_back(optbw_t_to_use(k));
       map_t(k) = id;
     } else { map_t(k) = jt->second; }
   }
   arma::vec svec_u(svec_u_v), obs_u(obs_u_v), tvec_u(tvec_u_v), obt_u(obt_u_v);
   int Us = svec_u.n_elem, Ut = tvec_u.n_elem;

   // Estimate mean function (once per distinct tuple, then expand by row)
   arma::mat mat_mean_s_u = estimate_mean_cpp(data, svec_u, Rcpp::wrap(obs_u), R_NilValue, kernel_name);
   arma::mat mat_mean_t_u = estimate_mean_cpp(data, tvec_u, Rcpp::wrap(obt_u), R_NilValue, kernel_name);
   arma::mat mat_mean_s = mat_mean_s_u.rows(map_s);
   arma::mat mat_mean_t = mat_mean_t_u.rows(map_t);
   arma::vec muhat_s = mat_mean_s.col(5);
   arma::vec muhat_t = mat_mean_t.col(5);

   // // Estimate the error sd for diagonal correction
   arma::vec vec_sig_s(n, arma::fill::zeros);
   arma::vec vec_sig_t(n, arma::fill::zeros);
   if (lag == 0 && correct_diagonal) {
     // The error sd depends only on the point; svec/tvec repeat few distinct
     // values (an m x m grid has only m distinct points), so estimate once per
     // distinct value and expand. Bit-identical (sigma at a point is independent
     // of the other points).
     arma::mat mat_sig_s = estimate_sigma_cpp(data, svec_u);
     arma::mat mat_sig_t = estimate_sigma_cpp(data, tvec_u);
     arma::vec sig_s_u = mat_sig_s.col(1);
     arma::vec sig_t_u = mat_sig_t.col(1);
     vec_sig_s = sig_s_u.elem(map_s);
     vec_sig_t = sig_t_u.elem(map_t);
   }

   // Compute the autocovariance function
   arma::vec PN_lag(n, arma::fill::zeros);
   arma::vec autocov_numerator(n, arma::fill::zeros);
   arma::vec diag_correct_numerator(n, arma::fill::zeros);

   for (int i = 0; i < n_curve - lag; ++i) {
     arma::uvec indices_cur_s = arma::find(data_mat.col(0) == unique_id_curve(i));
     arma::uvec indices_cur_t = arma::find(data_mat.col(0) == unique_id_curve(i + lag));
     arma::vec Tnvec_s = data_mat(indices_cur_s, arma::uvec({1}));
     arma::vec Ynvec_s = data_mat(indices_cur_s, arma::uvec({2}));
     arma::vec Tnvec_t = data_mat(indices_cur_t, arma::uvec({1}));
     arma::vec Ynvec_t = data_mat(indices_cur_t, arma::uvec({2}));

     // Smooth using Nadaraya-Watson estimator (once per distinct tuple)
     arma::vec Xhat_s_u = estimate_nw_cpp(Ynvec_s, Tnvec_s, svec_u, obs_u, kernel_name);
     arma::vec Xhat_t_u = estimate_nw_cpp(Ynvec_t, Tnvec_t, tvec_u, obt_u, kernel_name);
     Xhat_s_u.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
     Xhat_t_u.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
     arma::vec Xhat_s = Xhat_s_u.elem(map_s);
     arma::vec Xhat_t = Xhat_t_u.elem(map_t);

     // Compute the vector p_n(t;h) and update P_N(t;h) (once per distinct tuple)
     arma::vec pn_s_u(Us), pn_t_u(Ut);
     for (int j = 0; j < Us; ++j) pn_s_u(j) = arma::any(arma::abs(Tnvec_s - svec_u(j)) <= obs_u(j)) ? 1.0 : 0.0;
     for (int j = 0; j < Ut; ++j) pn_t_u(j) = arma::any(arma::abs(Tnvec_t - tvec_u(j)) <= obt_u(j)) ? 1.0 : 0.0;
     arma::vec pn_s = pn_s_u.elem(map_s);
     arma::vec pn_t = pn_t_u.elem(map_t);
     PN_lag += pn_s % pn_t;

     // Estimate the numerator function numerator
     if (center) {
       autocov_numerator += pn_s % pn_t % (Xhat_s - muhat_s) % (Xhat_t - muhat_t);
     } else {
       autocov_numerator += pn_s % pn_t % Xhat_s % Xhat_t;
     }

     // Diagonal correction
     if (lag == 0 && correct_diagonal) {
       // Normalized NW weight vectors, once per distinct tuple
       std::vector<arma::vec> Ws(Us), Wt(Ut);
       for (int j = 0; j < Us; ++j) {
         arma::vec w = kernel_func((Tnvec_s - svec_u(j)) / obs_u(j));
         w /= arma::accu(w);
         w.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
         Ws[j] = w;
       }
       for (int j = 0; j < Ut; ++j) {
         arma::vec w = kernel_func((Tnvec_t - tvec_u(j)) / obt_u(j));
         w /= arma::accu(w);
         w.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
         Wt[j] = w;
       }
       // Compute the weight product for each (s,t) pair from the cached vectors
       arma::vec weight_product(n);
       for (int k = 0; k < n; ++k)
         weight_product(k) = arma::accu(Ws[map_s(k)] % Wt[map_t(k)]);

       // Estimate the numerator of the diagonal part
       diag_correct_numerator += vec_sig_s % vec_sig_t % pn_s % pn_t % weight_product;
     }
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
   // // Init output
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
   res_autocov.col(8) = mat_mean_s.col(4);
   res_autocov.col(9) = mat_mean_s.col(5);
   res_autocov.col(10) = mat_mean_t.col(4);
   res_autocov.col(11) = mat_mean_t.col(5);
   res_autocov.col(12) = PN_lag;
   res_autocov.col(13) = autocovhat;

   if (lag == 0) {
     // If lag == 0, fill the lower part of the covariance matrix
     for (int k = 0; k < n_couple ; ++k) {
       arma::uvec idx_cov_st_upper = arma::find((res_autocov.col(0) == s(k)) % (res_autocov.col(1) == t(k)));
       if (! idx_cov_st_upper.is_empty()) {
         mat_res_autocov.row(k) = res_autocov.row(idx_cov_st_upper(0));
       } else {
         arma::uvec idx_cov_st_lower = arma::find((res_autocov.col(1) == s(k)) % (res_autocov.col(0) == t(k)));
         mat_res_autocov.row(k) = {res_autocov(idx_cov_st_lower(0), 1), res_autocov(idx_cov_st_lower(0), 0),
                                   res_autocov(idx_cov_st_lower(0), 3), res_autocov(idx_cov_st_lower(0), 2),
                                   res_autocov(idx_cov_st_lower(0), 6), res_autocov(idx_cov_st_lower(0), 7),
                                   res_autocov(idx_cov_st_lower(0), 4), res_autocov(idx_cov_st_lower(0), 5),
                                   res_autocov(idx_cov_st_lower(0), 10), res_autocov(idx_cov_st_lower(0), 11),
                                   res_autocov(idx_cov_st_lower(0), 8), res_autocov(idx_cov_st_lower(0), 9),
                                   res_autocov(idx_cov_st_lower(0), 12), res_autocov(idx_cov_st_lower(0), 13)};
       }

     }
   } else {
     mat_res_autocov = res_autocov;
   }

   return sort_by_columns(mat_res_autocov, 0, 1);

 }




