#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "02_smoothing_rcpp.h"
#include "03_estimate_locreg_rcpp.h"
#include "07_estimate_constants_cpp.h"
#include "04_estimate_mean_cpp.h"
using namespace Rcpp;
using namespace arma;

//' Estimate the risk of the covariance or autocovariance function
 //'
 //' Estimate the risk function of the lag-\eqn{\ell}, \eqn{\ell} = 0, 1,..., autocovariance function estimator of Maissoro et al. (2024).
 //'
 //' @param data A DataFrame containing the columns "id_curve", "tobs", and "X". Typically, the output of the function \link{.format_data}.
 //' @param s \code{vector (numeric)}. First argument of the autocovariance function.
 //' It corresponds to the observation points \code{s} in the pair (\code{s}, \code{t}).
 //' It has to be of the same length as the \code{t}
 //' @param t \code{vector (numeric)}. Second argument of the autocovariance function.
 //' It corresponds to the observation points \code{t} in the pair (\code{s}, \code{t}).
 //' It has to be of the same length as the \code{s}.
 //' @param lag \code{integer (positive integer)}. Lag of the autocovariance.
 //' @param bw_grid \code{vector (numeric)}. The bandwidth grid in which the best smoothing parameter is selected for each pair (\code{s}, \code{t}).
 //' It can be \code{NULL} and that way it will be defined as an exponential grid of \eqn{N\times\lambda}.
 //' @param smooth_ker \code{function}. The kernel function of the Nadaraya-Watson estimator.
 //' Default \code{smooth_ker = epanechnikov}.
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
 //' @export
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
   arma::vec unique_id_curve = arma::unique(data_mat.col(0));
   double n_curve = unique_id_curve.n_elem;

   // Set the bandwidth grid
   arma::vec bw_grid_to_use;
   if (bw_grid.isNull()) {
     // Estimate lambda
     double lambdahat = arma::mean(hist(data_mat.col(0), unique_id_curve));
     double b0 = 2 * std::max(std::pow(n_curve * lambdahat, -1 / (2 * 0.05 + 1)),
                              std::pow(n_curve * lambdahat * lambdahat, -1 / (2 * 0.05 + 1))); // rate with minimum local exponent = 0.05
     double bK = 4 * std::max(std::pow(n_curve * lambdahat, -1 / (2 * 0.8 + 1)),
                              std::pow(n_curve * lambdahat * lambdahat, -1 / (2 * 0.8 + 1))); // rate with maximum local exponent = 0.8
     bw_grid_to_use = arma::logspace(log10(b0), log10(bK), 15);
   } else {
     bw_grid_to_use = as<arma::vec>(bw_grid);
   }

   int bw_size = bw_grid_to_use.size();


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

   // Definie result matrix
   int n_rows_res = 0;
   if (use_same_bw) {
     n_rows_res = bw_size * n;
   } else {
     n_rows_res = bw_size * bw_size * n;
   }
   arma::mat mat_res_risk(n_rows_res, 14);

   // Compute the risk for each (s,t) and each bandwidth in bw_grid
   if (use_same_bw) {
     for (int idx_bw = 0; idx_bw < bw_size; ++idx_bw) {
       // extract the bandwidth
       double bw = bw_grid_to_use[idx_bw];

       // Do the computation for each s, t
       for (int k = 0; k < n; ++k) {
         // Extract local regularity parameters
         arma::uvec idx_locreg_cur_s = arma::find(mat_locreg_s.col(0) == s(k));
         arma::uvec idx_locreg_cur_t = arma::find(mat_locreg_t.col(0) == t(k));
         double Hs = mat_locreg_s(idx_locreg_cur_s(0), 4);
         double Ls = mat_locreg_s(idx_locreg_cur_s(0), 5);
         double Ht = mat_locreg_t(idx_locreg_cur_t(0), 4);
         double Lt = mat_locreg_t(idx_locreg_cur_t(0), 5);

         // Init. output
         arma::vec pn_s_vec = arma::zeros(n_curve - lag);
         arma::vec pn_lag_t_vec = arma::zeros(n_curve - lag);
         double bias_term_num = 0;
         double variance_term_num = 0;

         for(int i = 0 ; i < n_curve - lag; ++i){
           // Extrat the current curve index data
           arma::mat mat_cur_s = data_mat.rows(arma::find(data_mat.col(0) == i + 1));
           arma::mat mat_cur_t = data_mat.rows(arma::find(data_mat.col(0) == i + 1 + lag));
           arma::vec Tnvec_s = mat_cur_s.col(1);
           arma::vec Tnvec_t = mat_cur_t.col(1);
           int Mn = Tnvec_s.size();
           int Mn_lag = Tnvec_t.size();

           // Compute the weight vectors for each s, t and for each bw
           // and replace replace non-finite values with 0
           ////  For the argument s
           arma::vec Tn_s_diff_over_bw = (Tnvec_s - s(k)) / bw;
           arma::vec wvec_s = kernel_func(Tn_s_diff_over_bw);
           wvec_s /= arma::accu(wvec_s);
           wvec_s.replace(arma::datum::nan, 0);
           wvec_s.replace(arma::datum::inf, 0);
           wvec_s.replace(-arma::datum::inf, 0);

           //// For the argument t
           arma::vec Tn_t_diff_over_bw = (Tnvec_t - t(k)) / bw;
           arma::vec wvec_t = kernel_func(Tn_t_diff_over_bw);
           wvec_t /= arma::accu(wvec_t);
           wvec_t.replace(arma::datum::nan, 0);
           wvec_t.replace(arma::datum::inf, 0);
           wvec_t.replace(-arma::datum::inf, 0);

           // Compute the maximums and c_n(t;h) and c_n(s,h)
           double wmax_s = wvec_s.max();
           double wmax_t = wvec_t.max();

           // Compute and store the vector \pi_n(s;h)
           arma::uvec idx_is_one_pi_s = arma::find(abs(Tn_s_diff_over_bw) <= 1);
           double pn_s = 0;
           if (!idx_is_one_pi_s.is_empty()) {
             pn_s = 1;
           }
           pn_s_vec(i) = pn_s;

           // Compute and store the vector \pi_n(t;h)
           arma::uvec idx_is_one_pi_t = arma::find(abs(Tn_t_diff_over_bw) <= 1);
           double pn_lag_t = 0;
           if (!idx_is_one_pi_t.is_empty()) {
             pn_lag_t = 1;
           }
           pn_lag_t_vec(i) = pn_lag_t;

           // Compute bias term numerator
           //// Extract moment
           arma::uvec cur_idx_mom_s = arma::find(mat_mom_s.col(0) == s(k));
           arma::uvec cur_idx_mom_t = arma::find(mat_mom_t.col(0) == t(k));
           double mom_s_square = mom_vec_s(cur_idx_mom_s(0));
           double mom_t_square = mom_vec_t(cur_idx_mom_t(0));
           //// other terms
           arma::vec Tn_s_2Hs = arma::pow(arma::abs(Tn_s_diff_over_bw), 2 * Hs);
           arma::vec Tn_t_2Ht = arma::pow(arma::abs(Tn_t_diff_over_bw), 2 * Ht);
           double bn_s = arma::sum(Tn_s_2Hs % arma::abs(wvec_s));
           double bn_t = arma::sum(Tn_t_2Ht % arma::abs(wvec_t));
           bias_term_num += 3 * mom_t_square * Ls * std::pow(bw, 2 * Hs) * pn_s * pn_lag_t * bn_s +
             3 * mom_s_square * Lt * std::pow(bw, 2 * Ht) * pn_s * pn_lag_t * bn_t;

           // Compute variance term numerator
           arma::uvec cur_idx_sig_s = arma::find(mat_sig_s.col(0) == s(k));
           arma::uvec cur_idx_sig_t = arma::find(mat_sig_t.col(0) == t(k));
           double sig_s_square = sig2_vec_s(cur_idx_sig_s(0));
           double sig_t_square = sig2_vec_t(cur_idx_sig_t(0));
           variance_term_num += 3 * sig_s_square * mom_t_square * pn_s * pn_lag_t * wmax_s +
             3 * sig_t_square * mom_s_square * pn_s * pn_lag_t * wmax_t +
             3 * sig_s_square * sig_t_square * pn_s * pn_lag_t * wmax_s * wmax_t;
         }

         // Compute P_N(t;h)
         double PN_lag = arma::accu(pn_s_vec % pn_lag_t_vec);

         // Compute bias term
         double bias_term = bias_term_num / PN_lag;

         // Compute variance term
         double variance_term = variance_term_num / (PN_lag * PN_lag);

         // Compute dependence term
         // Add the lag-0 autocovariance to the vector
         arma::uvec idx_lag0 = arma::find(mat_emp_autocov.col(0) == s(k) && mat_emp_autocov.col(1) == t(k) && mat_emp_autocov.col(3) == 0);
         arma::uvec idx_lag = arma::find(mat_emp_autocov.col(0) == s(k) && mat_emp_autocov.col(1) == t(k) && (mat_emp_autocov.col(3) != 0));
         double XsXt_var = mat_emp_autocov(idx_lag0(0), 5);
         arma::mat XsXt_mat_lr_var = mat_emp_autocov.rows(idx_lag);
         double dependence_term_num = XsXt_var + arma::accu(arma::abs(2 * XsXt_mat_lr_var.col(5)));
         double dependence_term = dependence_term_num  / PN_lag;

         // Autocovariance risk
         double autocov_risk = 2 * (bias_term + variance_term + dependence_term);
         mat_res_risk.row(idx_bw * n + k) = {s(k), t(k), bw, bw, PN_lag, h(0), Hs, Ls, Ht, Lt, bias_term, variance_term, dependence_term, autocov_risk};
       }
     }

   } else {
     // If two bandwidth are used
     // arma::mat mat_res_risk(bw_size * bw_size * n, 8);
     for (int idx_bw_s = 0; idx_bw_s < bw_size; ++idx_bw_s) {
       arma::mat mat_res_risk_cur_bw(bw_size * n, 14);
       for (int idx_bw_t = 0; idx_bw_t < bw_size; ++idx_bw_t) {
         // extract the bandwidth
         double bw_s = bw_grid_to_use[idx_bw_s];
         double bw_t = bw_grid_to_use[idx_bw_s];

         // Do the computation for each s, t
         for (int k = 0; k < n; ++k) {
           // Extract local regularity parameters
           arma::uvec idx_locreg_cur_s = arma::find(mat_locreg_s.col(0) == s(k));
           arma::uvec idx_locreg_cur_t = arma::find(mat_locreg_t.col(0) == t(k));
           double Hs = mat_locreg_s(idx_locreg_cur_s(0), 4);
           double Ls = mat_locreg_s(idx_locreg_cur_s(0), 5);
           double Ht = mat_locreg_t(idx_locreg_cur_t(0), 4);
           double Lt = mat_locreg_t(idx_locreg_cur_t(0), 5);

           // Init. output
           arma::vec pn_s_vec(n_curve - lag);
           arma::vec pn_lag_t_vec(n_curve - lag);
           double bias_term_num = 0;
           double variance_term_num = 0;

           for(int i = 0 ; i < n_curve - lag; ++i){
             // Extrat the current curve index data
             arma::mat mat_cur_s = data_mat.rows(arma::find(data_mat.col(0) == i + 1));
             arma::mat mat_cur_t = data_mat.rows(arma::find(data_mat.col(0) == i + 1 + lag));
             arma::vec Tnvec_s = mat_cur_s.col(1);
             arma::vec Tnvec_t = mat_cur_t.col(1);
             int Mn = Tnvec_s.size();
             int Mn_lag = Tnvec_t.size();

             // Compute the weight vectors for each s, t and for each bw
             // At the we replace replace non-finite values with 0
             ////  For the argument s
             arma::vec Tn_s_diff_over_bw = (Tnvec_s - s(k)) / bw_s;
             arma::vec wvec_s = kernel_func(Tn_s_diff_over_bw);
             wvec_s /= arma::accu(wvec_s);
             wvec_s.replace(arma::datum::nan, 0);
             wvec_s.replace(arma::datum::inf, 0);
             wvec_s.replace(-arma::datum::inf, 0);

             //// For the argument t
             arma::vec Tn_t_diff_over_bw = (Tnvec_t - t(k)) / bw_t;
             arma::vec wvec_t = kernel_func(Tn_t_diff_over_bw);
             wvec_t /= arma::accu(wvec_t);
             wvec_t.replace(arma::datum::nan, 0);
             wvec_t.replace(arma::datum::inf, 0);
             wvec_t.replace(-arma::datum::inf, 0);

             // Compute the maximums and c_n(t;h) and c_n(s,h)
             double wmax_s = wvec_s.max();
             double wmax_t = wvec_t.max();

             // Compute and store the vector \pi_n(s;h)
             arma::uvec idx_is_one_pi_s = arma::find(abs(Tn_s_diff_over_bw) <= 1);
             double pn_s = 0;
             if (!idx_is_one_pi_s.is_empty()) {
               pn_s = 1;
             }
             pn_s_vec(i) = pn_s;

             // Compute and store the vector \pi_n(t;h)
             arma::uvec idx_is_one_pi_t = arma::find(abs(Tn_t_diff_over_bw) <= 1);
             double pn_lag_t = 0;
             if (!idx_is_one_pi_t.is_empty()) {
               pn_lag_t = 1;
             }
             pn_lag_t_vec(i) = pn_lag_t;

             // Compute bias term numerator
             //// Extract moment
             arma::uvec cur_idx_mom_s = arma::find(mat_mom_s.col(0) == s(k));
             arma::uvec cur_idx_mom_t = arma::find(mat_mom_t.col(0) == t(k));
             double mom_s_square = mom_vec_s(cur_idx_mom_s(0));
             double mom_t_square = mom_vec_t(cur_idx_mom_t(0));
             //// other terms
             arma::vec Tn_s_2Hs = arma::pow(arma::abs(Tn_s_diff_over_bw), 2 * Hs);
             arma::vec Tn_t_2Ht = arma::pow(arma::abs(Tn_t_diff_over_bw), 2 * Ht);
             double bn_s = arma::sum(Tn_s_2Hs % arma::abs(wvec_s));
             double bn_t = arma::sum(Tn_t_2Ht % arma::abs(wvec_t));
             bias_term_num += 4 * mom_t_square * Ls * std::pow(bw_s, 2 * Hs) * pn_s * pn_lag_t * bn_s +
               4 * mom_s_square * Lt * std::pow(bw_t, 2 * Ht) * pn_s * pn_lag_t * bn_t;

             // Compute variance term numerator
             arma::uvec cur_idx_sig_s = arma::find(mat_sig_s.col(0) == s(k));
             arma::uvec cur_idx_sig_t = arma::find(mat_sig_t.col(0) == t(k));
             double sig_s_square = sig2_vec_s(cur_idx_sig_s(0));
             double sig_t_square = sig2_vec_t(cur_idx_sig_t(0));
             variance_term_num += 4 * sig_s_square * mom_t_square * pn_s * pn_lag_t * wmax_s +
               4 * sig_t_square * mom_s_square * pn_s * pn_lag_t * wmax_t +
               4 * sig_s_square * sig_t_square * pn_s * pn_lag_t * wmax_s * wmax_t;
           }

           // Compute P_N(t;h)
           double PN_lag = arma::accu(pn_s_vec % pn_lag_t_vec);

           // Compute bias term
           double bias_term = bias_term_num / PN_lag;

           // Compute variance term
           double variance_term = variance_term_num / (PN_lag * PN_lag);

           // Compute dependence term
           // Add the lag-0 autocovariance to the vector
           arma::uvec idx_lag0 = arma::find(mat_emp_autocov.col(0) == s(k) && mat_emp_autocov.col(1) == t(k) && mat_emp_autocov.col(3) == 0);
           arma::uvec idx_lag = arma::find(mat_emp_autocov.col(0) == s(k) && mat_emp_autocov.col(1) == t(k) && mat_emp_autocov.col(3) != 0);
           double XsXt_var = mat_emp_autocov(idx_lag0(0), 5);
           arma::mat XsXt_mat_lr_var = mat_emp_autocov.rows(idx_lag);
           double dependence_term_num = XsXt_var + arma::accu(2 * arma::abs(XsXt_mat_lr_var.col(5)));
           double dependence_term = dependence_term_num  / PN_lag;

           // Autocovariance risk
           double autocov_risk = 2 * (bias_term + variance_term + dependence_term);
           mat_res_risk_cur_bw.row(idx_bw_t * n + k) = {s(k), t(k), bw_s, bw_t, PN_lag, h(0), Hs, Ls, Ht, Lt, bias_term, variance_term, dependence_term, autocov_risk};
         }
       }
       mat_res_risk(span(idx_bw_s * bw_size * n, (idx_bw_s + 1) * bw_size * n - 1), span(0, 13)) = mat_res_risk_cur_bw;
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
 //'  by the first and then the second column.
 //'
 //' @param s A numeric vector.
 //' @param t A numeric vector of the same length as `s`.
 //'
 //' @return A matrix where each row contains a unique pair (s(k), t(k)) in
 //' upper triangular form, sorted by the first and then
 //' the second column.
 //'
 //' @examples
 //' \dontrun{
 //'   library(RcppArmadillo)
 //'   s <- c(0.2, 0.4, 0.6, 0.4)
 //'   t <- c(0.4, 0.2, 0.6, 0.6)
 //'   get_upper_tri_couple(s, t)
 //' }
 //'
 arma::mat get_upper_tri_couple(const arma::vec& s, const arma::vec& t) {
   int n = t.size();
   arma::mat mat_st(n, 2);
   for (int k = 0; k < n; ++k) {
     mat_st.row(k) = {std::min(s(k), t(k)), std::max(s(k), t(k))};
   }

   // Remove duplicated rows
   std::unordered_set<std::string> seen;
   std::vector<arma::uword> unique_indices;

   for (arma::uword i = 0; i < mat_st.n_rows; ++i) {
     std::string row_str = "";
     for (arma::uword j = 0; j < mat_st.n_cols; ++j) {
       row_str += std::to_string(mat_st(i, j)) + ",";
     }
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
   // Step 1: Create a combined key
   arma::uvec combined_key = arma::conv_to<arma::uvec>::from(unique_mat.col(0) * 1e6 + unique_mat.col(1));

   // Step 2: Get the sort indices based on the combined key
   arma::uvec sort_indices = arma::sort_index(combined_key);

   // Step 3: Reorder the rows of the matrix using the sort indices
   arma::mat sorted_mat = unique_mat.rows(sort_indices);

   return sorted_mat;
 }

 //' Diagonal Correction of a Matrix
 //'
 //' This function corrects the diagonal elements of a matrix based on specified columns.
 //'
 //' @param mat Input matrix.
 //' @param idx_col_s Column index for 's'.
 //' @param idx_col_t Column index for 't'.
 //' @param idx_col_diff_st Column index for the difference between 's' and 't'.
 //' @param idx_col_bw Column index for the bandwidth.
 //' @param idx_col_cov Column index for the covariance.
 //' @return A matrix with corrected diagonal elements.
 //' @export
 //'
 arma::mat diagonal_correct(const arma::mat& mat,
                            const arma::uword idx_col_s,
                            const arma::uword idx_col_t,
                            const arma::uword idx_col_diff_st,
                            const arma::uword idx_col_bw,
                            const arma::uword idx_col_cov) {
   // Initialize the result matrix
   arma::mat result = mat;

   // Function to process each unique value in a given column
   auto process_unique_column = [&](const arma::vec& col, const arma::vec& unique_values) {
     for (arma::uword v = 0; v < unique_values.n_elem; ++v) {
       double current_value = unique_values(v);

       // Subset the matrix for the current group
       arma::uvec indices = arma::find(col == current_value);
       arma::mat group_mat = mat.rows(indices);

       // Sort group_mat in decreasing order of the idx_col_diff_st column
       arma::uvec sort_indices = arma::sort_index(group_mat.col(idx_col_diff_st), "descend");
       group_mat = group_mat.rows(sort_indices);

       // Find the last occurrence where the idx_col_diff_st column is greater than or equal to the idx_col_bw column
       arma::uword last_occurrence = group_mat.n_rows;
       for (arma::uword i = 0; i < group_mat.n_rows; ++i) {
         float d = group_mat(i, idx_col_diff_st);
         float bw_d = group_mat(i, idx_col_bw);
         if (d >= bw_d) {
           last_occurrence = i;
         }
       }

       // If such an occurrence is found, replace the idx_col_cov column values for the group
       if (last_occurrence < group_mat.n_rows) {
         double replacement_value = group_mat(last_occurrence, idx_col_cov);
         arma::uvec indices_to_be_replaced = indices.elem(arma::find(mat(indices, arma::uvec({idx_col_diff_st})) < mat(indices, arma::uvec({idx_col_bw}))));
         result(indices_to_be_replaced, arma::uvec({idx_col_cov})).fill(replacement_value);
       }
     }
   };

   // Upper triangular correction for idx_col_t
   arma::vec tvec_unique = arma::unique(mat.col(idx_col_t));
   process_unique_column(mat.col(idx_col_t), tvec_unique);

   // Lower triangular correction for idx_col_s
   arma::vec svec_unique = arma::unique(mat.col(idx_col_s));
   process_unique_column(mat.col(idx_col_s), svec_unique);

   return result;
 }

 //' Build a full grid of pairs
 //'
 //' This function takes two vectors of parameters and returns a matrix containing
 //' all possible pairs of the parameters. Each row in the resulting matrix represents
 //' a unique pair where the first element is from the first vector and the second element
 //' is from the second vector.
 //'
 //' @param u A numeric vector of parameters for the first dimension.
 //' @param v A numeric vector of parameters for the second dimension.
 //' @return A matrix where each row is a pair from the two parameter grids.
 //'
 //' @examples
 //' \dontrun{
 //' u <- c(0.2, 0.4, 0.5)
 //' v <- c(0.7, 0.8)
 //' build_grid(u, v)
 //' }
 //'
 //' @export
 //'
 // [[Rcpp::export]]
 arma::mat build_grid(const arma::vec& u, const arma::vec& v) {
   int nu = u.n_elem;
   int nv = v.n_elem;
   arma::mat grid(nu * nv, 2);

   for (int i = 0; i < nu; ++i) {
     for (int j = 0; j < nv; ++j) {
       grid(i * nv + j, 0) = u(i);
       grid(i * nv + j, 1) = v(j);
     }
   }

   return grid;
 }


 //' Estimate the risk of the covariance or autocovariance function
 //'
 //' Estimate the risk function of the lag-\eqn{\ell}, \eqn{\ell} = 0, 1,..., autocovariance function estimator of Maissoro et al. (2024).
 //'
 //' @param data A \code{DataFrame} containing the data with columns \code{"id_curve"}, \code{"tobs"}, and \code{"X"}.
 //' @param s A numeric vector specifying time points \code{s} for which to estimate autocovariance.
 //' @param t A numeric vector specifying time points \code{t} for which to estimate autocovariance.
 //' @param lag An integer specifying the lag value for autocovariance.
 //' @param param_grid A numeric vector of grid points for parameter estimation.
 //' @param optbw_s Optional numeric vector specifying optimal bandwidths for \code{s}. Default is \code{NULL}.
 //' @param optbw_t Optional numeric vector specifying optimal bandwidths for \code{t}. Default is \code{NULL}.
 //' @param bw_grid Optional numeric vector of bandwidth grid values. Default is \code{NULL}.
 //' @param use_same_bw A logical value indicating if the same bandwidth should be used for \code{s} and \code{t}. Default is \code{FALSE}.
 //' @param center A logical value indicating if the data should be centered before estimation. Default is \code{TRUE}.
 //' @param kernel_name A string specifying the kernel to use for estimation. Supported values are \code{"epanechnikov"}, \code{"biweight"}, \code{"triweight"}, \code{"tricube"}, \code{"triangular"}, \code{"uniform"}. Default is \code{"epanechnikov"}.
 //'
 //' @return A numeric matrix with six columns: \code{s}, \code{t}, \code{optbw_s}, \code{optbw_t}, the estimated autocovariance, and additional local regression results.
 //' @return A \code{matrix} containing the following fourteen columns in order:
 //'          \itemize{
 //'            \item{s : The first argument of the autocovariance function.}
 //'            \item{t : The second argument of the autocovariance function.}
 //'            \item{optbw_s : The optimal bandwidth for the first argument of the autocovariance function. If \code{use_same_bw = TRUE}, the same bandwidth candidate is used for \code{s} and for \code{t}, so the 3rd and 4th columns contain the same values.}
 //'            \item{optbw_t : The optimal bandwidth for the second argument of the autocovariance function.}
 //'            \item{Hs : The estimates of the local exponent for each t. It corresponds to \eqn{H_s}.}
 //'            \item{Ls : The estimates of the Hölder constant for each t. It corresponds to \eqn{L_s^2}.}
 //'            \item{Ht : The estimates of the local exponent for each t. It corresponds to \eqn{H_t}.}
 //'            \item{Lt : The estimates of the Hölder constant for each t. It corresponds to \eqn{L_t^2}.}
 //'            \item{PNs : The number of curves used to estimate the mean at s. It corresponds to \eqn{P_N(s;h)}.}
 //'            \item{muhat_s : The estimates of the mean at s. It corresponds to \eqn{\widehat{\mu}_N(s;h)}.}
 //'            \item{PNt : The number of curves used to estimate the mean at t. It corresponds to \eqn{P_N(t;h)}.}
 //'            \item{muhat_t : The estimates of the mean at t. It corresponds to \eqn{\widehat{\mu}_N(t;h)}.}
 //'            \item{PNl : The number of curves used to estimate the autocovariance at (s,t). It corresponds to \eqn{P_{N,\ell}(s,t;h_s, h_t)}.}
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
 //' param_grid <- seq(0.1, 1, length.out = 5)
 //' result <- estimate_autocov_cpp(data, s, t, lag = 0, param_grid = param_grid)
 //' }
 //'
 //' @export
 //'
 // [[Rcpp::export]]
 arma::mat estimate_autocov_cpp(const Rcpp::DataFrame data,
                                const arma::vec s,
                                const arma::vec t,
                                const int lag,
                                const arma::vec param_grid,
                                const Rcpp::Nullable<arma::vec> optbw_s = R_NilValue,
                                const Rcpp::Nullable<arma::vec> optbw_t = R_NilValue,
                                const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                const bool use_same_bw = false,
                                const bool center = true,
                                const std::string kernel_name = "epanechnikov"){
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

   // Convert data to matrix
   arma::mat data_mat(data.nrows(), 3);
   data_mat.col(0) = as<arma::vec>(data["id_curve"]);
   data_mat.col(1) = as<arma::vec>(data["tobs"]);
   data_mat.col(2) = as<arma::vec>(data["X"]);
   arma::vec unique_id_curve = arma::unique(data_mat.col(0));
   double n_curve = unique_id_curve.n_elem;

   // the vector s and t to be used in the estimation procedure
   arma::vec svec;
   arma::vec tvec;

   // If lag = 0, the covariance is symmetric
   // We can use this property to speed up the computation time
   // When we have to estimate the bandwidth parameters
   if (lag == 0 && (optbw_s.isNull() || optbw_t.isNull())) {
     arma::mat mat_st_upper = get_upper_tri_couple(s, t);
     // Set observation points size
     svec.set_size(mat_st_upper.n_rows);
     tvec.set_size(mat_st_upper.n_rows);
     // Affect the vector
     svec = mat_st_upper.col(0);
     tvec = mat_st_upper.col(1);
   } else {
     svec = s;
     tvec = t;
   }
   // Take number of couples (s,t)
   int n = tvec.n_elem;

   // Estimate optimal risk function / Vérify the
   arma::vec optbw_s_to_use(n);
   arma::vec optbw_t_to_use(n);
   arma::mat mat_locreg(n, 6);

   if (optbw_s.isNull() || optbw_t.isNull()) {
     // Estimate the parameters on a grid if tvec.size() > param_grid.size() * param_grid.size()
     int n_grid = param_grid.size() * param_grid.size();

     if (n >= n_grid) {
       arma::mat grid = build_grid(param_grid, param_grid);
       arma::vec sgrid = grid.col(0);
       arma::vec tgrid = grid.col(1);

       // Estimate the risk function on a grid
       // And return the optimum
       arma::mat mat_risk_grid = estimate_autocov_risk_cpp(data, sgrid, tgrid, lag, bw_grid, use_same_bw, center, kernel_name);
       arma::mat mat_risk_min_grid(n_grid, mat_risk_grid.n_cols);

       for (int g = 0; g < n_grid; ++g) {
         // Find rows in mat_risk where the first column equals t(k)
         arma::uvec idx_risk_cur = arma::find(mat_risk_grid.col(0) == sgrid(g) && mat_risk_grid.col(1) == tgrid(g));

         // Extract the risk column for those rows
         arma::vec risk = mat_risk_grid(idx_risk_cur, arma::uvec({13}));

         // Find the minimum index in the risk column, ignoring non-finite values
         arma::uword idx_min = arma::index_min(risk.elem(arma::find_finite(risk)));

         // return the minimun matrix
         mat_risk_min_grid.row(g) = mat_risk_grid.row(idx_risk_cur(idx_min));
       }

       // Matching using the nearest neighbour strategy
       for (int k = 0; k < n; ++k) {
         // Find rows in mat_risk where the first column equals t(k)
         arma::vec dist = arma::square(mat_risk_min_grid.col(0) - svec(k)) + arma::square(mat_risk_min_grid.col(1) - tvec(k));

         // Find the minimum index in the risk column
         arma::uword idx_min_dist = arma::index_min(dist);

         // Extract values corresponding to the minimum risk index
         optbw_s_to_use(k) = mat_risk_min_grid(idx_min_dist, 2);
         optbw_t_to_use(k) = mat_risk_min_grid(idx_min_dist, 3);
         mat_locreg.row(k) = {svec(k), tvec(k),
                              mat_risk_min_grid(idx_min_dist, 6), mat_risk_min_grid(idx_min_dist, 7),
                              mat_risk_min_grid(idx_min_dist, 8), mat_risk_min_grid(idx_min_dist, 9)};
       }

     } else {

       // Estimate risk function
       arma::mat mat_risk = estimate_autocov_risk_cpp(data, svec, tvec, lag, bw_grid, use_same_bw, center, kernel_name);

       for (int k = 0; k < n; ++k) {
         // Find rows in mat_risk where the first column equals t(k)
         arma::uvec idx_risk_cur = arma::find(mat_risk.col(0) == svec(k) && mat_risk.col(1) == tvec(k));

         // Extract the risk column for those rows
         arma::vec risk = mat_risk(idx_risk_cur, arma::uvec({13}));

         // Find the minimum index in the risk column, ignoring non-finite values
         arma::uword idx_min = arma::index_min(risk.elem(arma::find_finite(risk)));

         // Extract values corresponding to the minimum risk index
         optbw_s_to_use(k) = mat_risk(idx_risk_cur(idx_min), 2);
         optbw_t_to_use(k) = mat_risk(idx_risk_cur(idx_min), 3);
         mat_locreg.row(k) = {mat_risk(idx_risk_cur(idx_min), 0), mat_risk(idx_risk_cur(idx_min), 1),
                              mat_risk(idx_risk_cur(idx_min), 6), mat_risk(idx_risk_cur(idx_min), 7),
                              mat_risk(idx_risk_cur(idx_min), 8), mat_risk(idx_risk_cur(idx_min), 9)};
       }
     }
   } else {
     arma::vec optbw_s_cur = as<arma::vec>(optbw_s);
     arma::vec optbw_t_cur = as<arma::vec>(optbw_t);
     if (optbw_s_cur.size() != n || optbw_t_cur.size() != n) {
       stop("If 'optbw_s' and 'optbw_t' are not NULL, they must be the same length as 's' and as 't'.");
     } else {
       optbw_s_to_use = optbw_s_cur;
       optbw_t_to_use = optbw_t_cur;
       mat_locreg.col(0) = svec;
       mat_locreg.col(1) = tvec;
       mat_locreg.col(2) = arma::zeros(tvec.size());
       mat_locreg.col(3) = arma::zeros(tvec.size());
       mat_locreg.col(5) = arma::zeros(tvec.size());
       mat_locreg.col(5) = arma::zeros(tvec.size());
     }
   }

   // Estimate mean function
   arma::mat mat_mean_s = estimate_mean_cpp(data, svec, Rcpp::wrap(optbw_s_to_use), R_NilValue, kernel_name);
   arma::mat mat_mean_t = estimate_mean_cpp(data, tvec, Rcpp::wrap(optbw_t_to_use), R_NilValue, kernel_name);
   arma::vec muhat_s = mat_mean_s.col(5);
   arma::vec muhat_t = mat_mean_t.col(5);

   // Compute the autocovariance function
   arma::vec PN_lag(n, fill::zeros);
   arma::vec autocov_numerator(n, fill::zeros);

   for(int i = 0 ; i < n_curve - lag; ++i){
     // Extrat the current curve index data
     arma::uvec indices_cur_s = arma::find(data_mat.col(0) == i + 1);
     arma::uvec indices_cur_t = arma::find(data_mat.col(0) == i + 1 + lag);
     arma::vec Tnvec_s = data_mat(indices_cur_s, arma::uvec({1}));
     arma::vec Ynvec_s = data_mat(indices_cur_s, arma::uvec({2}));
     arma::vec Tnvec_t = data_mat(indices_cur_t, arma::uvec({1}));
     arma::vec Ynvec_t = data_mat(indices_cur_t, arma::uvec({2}));

     // Smooth using Nadaraya-Watson estimator
     arma::vec Xhat_s = estimate_nw_cpp(Ynvec_s, Tnvec_s, svec, optbw_s_to_use, kernel_name);
     arma::vec Xhat_t = estimate_nw_cpp(Ynvec_t, Tnvec_t, tvec, optbw_t_to_use, kernel_name);
     Xhat_s.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
     Xhat_t.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);

     // Compute the vector p_n(t;h) and update P_N(t;h)
     arma::vec pn_s = arma::regspace(0, n - 1);
     arma::vec pn_t = arma::regspace(0, n - 1);
     pn_s.transform([svec, optbw_s_to_use, Tnvec_s](int j) { return  arma::find(abs((Tnvec_s - svec(j))) <= optbw_s_to_use(j)).is_empty() ? 0 : 1 ;});
     pn_t.transform([tvec, optbw_t_to_use, Tnvec_t](int j) { return  arma::find(abs((Tnvec_t - tvec(j))) <= optbw_t_to_use(j)).is_empty() ? 0 : 1 ;});
     PN_lag += pn_s % pn_t;

     // Estimate the mean function numerator
     if (center) {
       autocov_numerator += pn_s % pn_t % (Xhat_s - muhat_s) % (Xhat_t - muhat_t);
     } else {
       autocov_numerator += pn_s % pn_t % Xhat_s % Xhat_t;
     }
   }

   // Compute \widehat \Gamma_{N, \ell}(s,t;h_s, h_t)
   arma::vec autocovhat;
   if (center) {
     autocovhat = autocov_numerator / PN_lag;
   } else {
     autocovhat = autocov_numerator / PN_lag - muhat_s % muhat_t;
   }

   // Init the output
   int n_couple = t.size();
   arma::mat mat_res_autocov(n_couple, 14);

   // Diagonal correction
   if (lag == 0) {

     // Put covariance and max(h_s, h_t) in a matrix
     arma::mat mat_cov(n, 5);
     mat_cov.col(0) = svec;
     mat_cov.col(1) = tvec;
     mat_cov.col(2) = arma::abs(svec - tvec);
     mat_cov.col(3) = 0.5 * (optbw_s_to_use + optbw_t_to_use + arma::abs(optbw_s_to_use - optbw_t_to_use));
     mat_cov.col(4) = autocovhat;

     // Correct the diagnal
     arma::mat mat_cov_corrected = diagonal_correct(mat_cov, 0, 1, 2, 3, 4);
     int n_upper_tri = mat_cov_corrected.n_rows;

     // Initialize the result matrix with double the number of rows
     arma::mat res_cov(n_upper_tri * 2, 5);

     // Fill the upper triangular elements
     res_cov(arma::span(0, n_upper_tri - 1), 0) = mat_cov_corrected.col(0);  // s
     res_cov(arma::span(0, n_upper_tri - 1), 1) = mat_cov_corrected.col(1);  // t
     res_cov(arma::span(0, n_upper_tri - 1), 2) = optbw_s_to_use;            // bw_s
     res_cov(arma::span(0, n_upper_tri - 1), 3) = optbw_t_to_use;            // bw_t
     res_cov(arma::span(0, n_upper_tri - 1), 4) = mat_cov_corrected.col(4);  // cov

     // Fill the lower triangular elements
     res_cov(arma::span(n_upper_tri, 2 * n_upper_tri - 1), 0) = mat_cov_corrected.col(1);  // t
     res_cov(arma::span(n_upper_tri, 2 * n_upper_tri - 1), 1) = mat_cov_corrected.col(0);  // s
     res_cov(arma::span(n_upper_tri, 2 * n_upper_tri - 1), 2) = optbw_t_to_use;            // bw_t
     res_cov(arma::span(n_upper_tri, 2 * n_upper_tri - 1), 3) = optbw_s_to_use;            // bw_s
     res_cov(arma::span(n_upper_tri, 2 * n_upper_tri - 1), 4) = mat_cov_corrected.col(4);  // cov

     // return the covariance result
     //// P_{N, \ell}(s,t; h_s, h_t)
     arma::mat mat_PN_lag(n, 3);
     mat_PN_lag.col(0) = svec;
     mat_PN_lag.col(1) = tvec;
     mat_PN_lag.col(2) = PN_lag;

     //// Mean matrix
     arma::mat mat_mean(n, 6);
     mat_mean.col(0) = svec;
     mat_mean.col(1) = tvec;
     mat_mean.col(2) = mat_mean_s.col(4);
     mat_mean.col(3) = mat_mean_s.col(5);
     mat_mean.col(4) = mat_mean_t.col(4);
     mat_mean.col(5) = mat_mean_t.col(5);

     for (int k = 0; k < n_couple ; ++k) {
       arma::uvec idx_cov_st = arma::find(res_cov.col(0) == s(k) && res_cov.col(1) == t(k));
       arma::uvec idx_locreg_st = arma::find((mat_locreg.col(0) == s(k) && mat_locreg.col(1) == t(k)) ||
         (mat_locreg.col(0) == t(k) && mat_locreg.col(1) == s(k))); // For the transpose
       arma::uvec idx_PN_st = arma::find((mat_PN_lag.col(0) == s(k) && mat_PN_lag.col(1) == t(k)) ||
         (mat_PN_lag.col(0) == t(k) && mat_PN_lag.col(1) == s(k)));  // For the transpose
       arma::uvec idx_muhat_st = arma::find((mat_mean.col(0) == s(k) && mat_mean.col(1) == t(k)) ||
         (mat_mean.col(0) == t(k) && mat_mean.col(1) == s(k)));

       mat_res_autocov.row(k) = {res_cov(idx_cov_st(0), 0),
                                 res_cov(idx_cov_st(0), 1),
                                 res_cov(idx_cov_st(0), 2),
                                 res_cov(idx_cov_st(0), 3),
                                 mat_locreg(idx_locreg_st(0), 2),
                                 mat_locreg(idx_locreg_st(0), 3),
                                 mat_locreg(idx_locreg_st(0), 4),
                                 mat_locreg(idx_locreg_st(0), 5),
                                 mat_mean(idx_muhat_st(0), 2),
                                 mat_mean(idx_muhat_st(0), 3),
                                 mat_mean(idx_muhat_st(0), 4),
                                 mat_mean(idx_muhat_st(0), 5),
                                 mat_PN_lag(idx_PN_st(0), 2),
                                 res_cov(idx_cov_st(0), 4)};
     }
   } else {
     // return the result
     mat_res_autocov.col(0) = svec;
     mat_res_autocov.col(1) = tvec;
     mat_res_autocov.col(2) = optbw_s_to_use;
     mat_res_autocov.col(3) = optbw_t_to_use;
     mat_res_autocov.col(4) = mat_locreg.col(2);
     mat_res_autocov.col(5) = mat_locreg.col(3);
     mat_res_autocov.col(6) = mat_locreg.col(4);
     mat_res_autocov.col(7) = mat_locreg.col(5);
     mat_res_autocov.col(8) = mat_mean_s.col(4);
     mat_res_autocov.col(9) = mat_mean_s.col(5);
     mat_res_autocov.col(10) = mat_mean_t.col(4);
     mat_res_autocov.col(11) = mat_mean_t.col(5);
     mat_res_autocov.col(12) = PN_lag;
     mat_res_autocov.col(13) = autocovhat;
   }

   return mat_res_autocov;
 }


