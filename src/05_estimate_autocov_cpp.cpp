#include <RcppArmadillo.h>
#include <omp.h>

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
   int n_rows_res = use_same_bw ? bw_size * n : bw_size * bw_size * n;
   arma::mat mat_res_risk(n_rows_res, 14);

   // Compute the risk for each (s,t) and each bandwidth in bw_grid
   if (use_same_bw) {
#pragma omp parallel for
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

         // Compute the weight vectors for each s, t and for each bw
         // and replace replace non-finite values with 0
         //// For the argument s
         arma::vec Tn_s_diff_over_bw = (data_mat.col(1) - s(k)) / bw;
         arma::vec wvec_s = kernel_func(Tn_s_diff_over_bw);
         wvec_s /= arma::accu(wvec_s);
         wvec_s.replace(arma::datum::nan, 0);
         wvec_s.replace(arma::datum::inf, 0);
         wvec_s.replace(-arma::datum::inf, 0);

         //// For the argument t
         arma::vec Tn_t_diff_over_bw = (data_mat.col(1) - t(k)) / bw;
         arma::vec wvec_t = kernel_func(Tn_t_diff_over_bw);
         wvec_t /= arma::accu(wvec_t);
         wvec_t.replace(arma::datum::nan, 0);
         wvec_t.replace(arma::datum::inf, 0);
         wvec_t.replace(-arma::datum::inf, 0);

         // Extract moment
         arma::uvec cur_idx_mom_s = arma::find(mat_mom_s.col(0) == s(k));
         arma::uvec cur_idx_mom_t = arma::find(mat_mom_t.col(0) == t(k));
         double mom_s_square = mom_vec_s(cur_idx_mom_s(0));
         double mom_t_square = mom_vec_t(cur_idx_mom_t(0));

         // Extract the estimates of the sd
         arma::uvec cur_idx_sig_s = arma::find(mat_sig_s.col(0) == s(k));
         arma::uvec cur_idx_sig_t = arma::find(mat_sig_t.col(0) == t(k));
         double sig_s_square = sig2_vec_s(cur_idx_sig_s(0));
         double sig_t_square = sig2_vec_t(cur_idx_sig_t(0));

         // Compute some element of the bias term
         arma::vec bn_s_vec = arma::pow(arma::abs(Tn_s_diff_over_bw), 2 * Hs) % arma::abs(wvec_s);
         arma::vec bn_t_vec = arma::pow(arma::abs(Tn_t_diff_over_bw), 2 * Ht) % arma::abs(wvec_t);

         // Init. output
         arma::vec pn_s_vec = arma::zeros(n_curve - lag);
         arma::vec pn_lag_t_vec = arma::zeros(n_curve - lag);
         double bias_term_num = 0;
         double variance_term_num = 0;

         for(int i = 0 ; i < n_curve - lag; ++i){
           // Exact the indexes : current and forward
           arma::uvec idx_i = arma::find(data_mat.col(0) == i + 1);
           arma::uvec idx_i_lag = arma::find(data_mat.col(0) == i + 1 + lag);

           // Extract weight
           arma::vec wvec_s_i = wvec_s(idx_i) / arma::accu(wvec_s(idx_i));
           arma::vec wvec_t_i_lag = wvec_t(idx_i_lag) / arma::accu(wvec_t(idx_i_lag));
           wvec_s_i.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
           wvec_t_i_lag.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);

           // Compute and store the vector \pi_n(s;h)
           pn_s_vec(i) = arma::find(abs(Tn_s_diff_over_bw(idx_i)) <= 1).is_empty() ? 0 : 1;

           // Compute and store the vector \pi_n(t;h)
           pn_lag_t_vec(i) = arma::find(abs(Tn_t_diff_over_bw(idx_i_lag)) <= 1).is_empty() ? 0 : 1;

           // Compute bias term numerator
           bias_term_num += 3 * mom_t_square * Ls * std::pow(bw, 2 * Hs) * pn_s_vec(i) * pn_lag_t_vec(i) * arma::sum(bn_s_vec(idx_i)) +
             3 * mom_s_square * Lt * std::pow(bw, 2 * Ht) * pn_s_vec(i) * pn_lag_t_vec(i) * arma::sum(bn_t_vec(idx_i_lag));

           // Compute variance term numerator
           variance_term_num += 3 * sig_s_square * mom_t_square * pn_s_vec(i) * pn_lag_t_vec(i) * wvec_s_i.max() +
             3 * sig_t_square * mom_s_square * pn_s_vec(i) * pn_lag_t_vec(i) * wvec_t_i_lag.max() +
             3 * sig_s_square * sig_t_square * pn_s_vec(i) * pn_lag_t_vec(i) * wvec_s_i.max() * wvec_t_i_lag.max();
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
#pragma omp parallel for
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

           // Compute the weight vectors for each s, t and for each bw
           // and replace replace non-finite values with 0
           //// For the argument s
           arma::vec Tn_s_diff_over_bw = (data_mat.col(1) - s(k)) / bw_s;
           arma::vec wvec_s = kernel_func(Tn_s_diff_over_bw);
           wvec_s /= arma::accu(wvec_s);
           wvec_s.replace(arma::datum::nan, 0);
           wvec_s.replace(arma::datum::inf, 0);
           wvec_s.replace(-arma::datum::inf, 0);

           //// For the argument t
           arma::vec Tn_t_diff_over_bw = (data_mat.col(1) - t(k)) / bw_t;
           arma::vec wvec_t = kernel_func(Tn_t_diff_over_bw);
           wvec_t /= arma::accu(wvec_t);
           wvec_t.replace(arma::datum::nan, 0);
           wvec_t.replace(arma::datum::inf, 0);
           wvec_t.replace(-arma::datum::inf, 0);

           // Extract moment
           arma::uvec cur_idx_mom_s = arma::find(mat_mom_s.col(0) == s(k));
           arma::uvec cur_idx_mom_t = arma::find(mat_mom_t.col(0) == t(k));
           double mom_s_square = mom_vec_s(cur_idx_mom_s(0));
           double mom_t_square = mom_vec_t(cur_idx_mom_t(0));

           // Extract the estimates of the sd
           arma::uvec cur_idx_sig_s = arma::find(mat_sig_s.col(0) == s(k));
           arma::uvec cur_idx_sig_t = arma::find(mat_sig_t.col(0) == t(k));
           double sig_s_square = sig2_vec_s(cur_idx_sig_s(0));
           double sig_t_square = sig2_vec_t(cur_idx_sig_t(0));

           // Compute some element of the bias term
           arma::vec bn_s_vec = arma::pow(arma::abs(Tn_s_diff_over_bw), 2 * Hs) % arma::abs(wvec_s);
           arma::vec bn_t_vec = arma::pow(arma::abs(Tn_t_diff_over_bw), 2 * Ht) % arma::abs(wvec_t);

           // Init. output
           arma::vec pn_s_vec = arma::zeros(n_curve - lag);
           arma::vec pn_lag_t_vec = arma::zeros(n_curve - lag);
           double bias_term_num = 0;
           double variance_term_num = 0;

           for(int i = 0 ; i < n_curve - lag; ++i){
             // Exact the indexes : current and forward
             arma::uvec idx_i = arma::find(data_mat.col(0) == i + 1);
             arma::uvec idx_i_lag = arma::find(data_mat.col(0) == i + 1 + lag);

             // Extract weight
             arma::vec wvec_s_i = wvec_s(idx_i) / arma::accu(wvec_s(idx_i));
             arma::vec wvec_t_i_lag = wvec_t(idx_i_lag) / arma::accu(wvec_t(idx_i_lag));
             wvec_s_i.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
             wvec_t_i_lag.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);

             // Compute and store the vector \pi_n(s;h)
             pn_s_vec(i) = arma::find(abs(Tn_s_diff_over_bw(idx_i)) <= 1).is_empty() ? 0 : 1;

             // Compute and store the vector \pi_n(t;h)
             pn_lag_t_vec(i) = arma::find(abs(Tn_t_diff_over_bw(idx_i_lag)) <= 1).is_empty() ? 0 : 1;

             // Compute bias term numerator
             bias_term_num += 4 * mom_t_square * Ls * std::pow(bw_s, 2 * Hs) * pn_s_vec(i) * pn_lag_t_vec(i) * arma::sum(bn_s_vec(idx_i)) +
               4 * mom_s_square * Lt * std::pow(bw_t, 2 * Ht) * pn_s_vec(i) * pn_lag_t_vec(i) * arma::sum(bn_t_vec(idx_i_lag));

             // Compute variance term numerator
             variance_term_num += 4 * sig_s_square * mom_t_square * pn_lag_t_vec(i) * pn_lag_t_vec(i) * wvec_s_i.max() +
               4 * sig_t_square * mom_s_square * pn_lag_t_vec(i) * pn_lag_t_vec(i) * wvec_t_i_lag.max() +
               4 * sig_s_square * sig_t_square * pn_lag_t_vec(i) * pn_lag_t_vec(i) * wvec_s_i.max() * wvec_t_i_lag.max();
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

 //' Get Upper Triangular Couples
 //'
 //' This function constructs a matrix of upper triangular couples from vectors s and t,
 //' ensuring that if max(s) < max(t), the missing values from t are included in s.
 //'
 //' @param s A vector of values.
 //' @param t A vector of values.
 //' @return A matrix with unique (s, t) pairs where s >= t, sorted by s in descending order and t in ascending order.
 //' @export
 // [[Rcpp::export]]
 arma::mat get_upper_tri_couple(const arma::vec& s, const arma::vec& t) {
   // Sort unique values of s in descending order and t in ascending order
   arma::vec s_unique = arma::sort(arma::unique(s), "descend");
   arma::vec t_unique = arma::sort(arma::unique(t), "ascend");

   // If max(s) < max(t), include the missing values from t_unique in s_unique
   double max_s = s_unique.max();
   double max_t = t_unique.max();
   if (max_s < max_t) {
     arma::uvec t_indices = arma::find(t_unique > max_s);
     s_unique = arma::join_vert(s_unique, t_unique.elem(t_indices));
     s_unique = arma::sort(s_unique, "descend");
   }

   // Generate pairs
   arma::mat mat_st(s_unique.n_elem * t_unique.n_elem, 2);
   int idx = 0;
   for (arma::uword i = 0; i < s_unique.n_elem; ++i) {
     for (arma::uword j = 0; j < t_unique.n_elem; ++j) {
       if (s_unique(i) >= t_unique(j)) {
         mat_st(idx, 0) = s_unique(i);
         mat_st(idx, 1) = t_unique(j);
         ++idx;
       }
     }
   }

   // Resize matrix to remove unused rows
   mat_st.resize(idx, 2);

   return mat_st;
 }

 //' Diagonal Correction of a Matrix
 //'
 //' This function corrects the diagonal elements of a matrix based on specified columns.
 //'
 //' @param mat_cov Input matrix.
 //' @param idx_col_s Column index for 's'.
 //' @param idx_col_t Column index for 't'.
 //' @param idx_col_diff_st Column index for the difference between 's' and 't'.
 //' @param idx_col_bw Column index for the bandwidth.
 //' @param idx_col_cov Column index for the covariance.
 //' @return A matrix with corrected diagonal elements.
 //' @export
 // [[Rcpp::export]]
 arma::mat diagonal_correct(const arma::mat& mat_cov,
                            const arma::uword idx_col_s,
                            const arma::uword idx_col_t,
                            const arma::uword idx_col_diff_st,
                            const arma::uword idx_col_bw,
                            const arma::uword idx_col_cov) {
   // Get unique sorted values of s and t
   arma::vec s_unique = arma::sort(arma::unique(mat_cov.col(idx_col_s)), "descend");
   arma::vec t_unique = arma::sort(arma::unique(mat_cov.col(idx_col_t)), "ascend");

   int ns = s_unique.size();
   int nt = t_unique.size();
   arma::mat mat_diff = arma::zeros(ns, nt);
   arma::mat mat_bw = arma::zeros(ns, nt);
   arma::mat mat_cov_reshape = arma::zeros(ns, nt);

   // Fill matrices with values from mat_cov
   for (int i = 0; i < ns; ++i) {
     for (int j = 0; j < nt; ++j) {
       arma::uvec idx = arma::find(mat_cov.col(idx_col_s) == s_unique(i) && mat_cov.col(idx_col_t) == t_unique(j));
       if (!idx.is_empty()) {
         mat_diff(i, j) = mat_cov(idx(0), idx_col_diff_st);
         mat_bw(i, j) = mat_cov(idx(0), idx_col_bw);
         mat_cov_reshape(i, j) = mat_cov(idx(0), idx_col_cov);
       }
     }
   }

   // Correct the matrix values based on the given conditions
   for (int ids = 0; ids < ns; ++ids) {
     int idt_replace = 0;
     for (int idt = 0; idt < nt - ids; ++idt) {
       if (mat_diff(ids, idt) > mat_bw(ids, idt)) {
         idt_replace = idt;
       }
     }
     double replace_cov = mat_cov_reshape(ids, idt_replace);

     // Replace diagonal and above diagonal values
     int n_diag_step = nt - ids - idt_replace;
     for (int idx_step = 0; idx_step < n_diag_step; ++idx_step) {
       if ((ids + idx_step) < ns && (idt_replace + idx_step) < nt) {
         mat_cov_reshape(ids + idx_step, idt_replace + idx_step) = replace_cov;
       }
     }

     // Handle case when ids is 0
     if (ids == 0) {
       for (int idx_rep = idt_replace; idx_rep < nt - ids; ++idx_rep) {
         if (idx_rep < nt) {
           mat_cov_reshape(ids, idx_rep) = replace_cov;
         }
       }
     }

     // Handle case when ids is ns - 1
     if (ids == ns - 1) {
       int ids_replace_backward = 0;
       for (int ids_backward = 0; ids_backward < ns; ++ids_backward) {
         if (mat_diff(ids_backward, 0) > mat_bw(ids_backward, 0)) {
           ids_replace_backward++;
         }
       }
       for (int idx_rep = ids_replace_backward; idx_rep < ns; ++idx_rep) {
         if (idx_rep < ns) {
           mat_cov_reshape(idx_rep, 0) = mat_cov_reshape(ids_replace_backward, 0);
         }
       }
     }
   }

   // Remove the data after the anti-diagonal
   arma::mat mat_cov_res = arma::zeros(ns, nt);
   for (int i = 0; i < ns; ++i) {
     for (int j = 0; j < nt - i; ++j) {
       if (j < nt) {
         mat_cov_res(i, j) = mat_cov_reshape(i, j);
       }
     }
   }

   // Initialize the result matrix
   arma::mat result = mat_cov;

   // Fill the result matrix
   for (arma::uword i = 0; i < ns; ++i) {
     for (arma::uword j = 0; j < nt - i; ++j) {
       arma::uvec idx_to_set = arma::find(mat_cov.col(idx_col_s) == s_unique(i) && mat_cov.col(idx_col_t) == t_unique(j));
       if (!idx_to_set.is_empty()) {
         result(idx_to_set, arma::uvec({idx_col_cov})).fill(mat_cov_res(i, j));
       }
     }
   }

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


 //' Estimate the autocovariance function
 //'
 //' This function estimates the lag-\eqn{\ell}, \eqn{\ell} = 0, 1,..., autocovariance function based on the methodology of Maissoro et al. (2024).
 //'
 //' @param data A \code{DataFrame} containing the data with columns \code{"id_curve"}, \code{"tobs"}, and \code{"X"}.
 //' @param s A numeric vector specifying time points \code{s} for which to estimate autocovariance.
 //' @param t A numeric vector specifying time points \code{t} for which to estimate autocovariance.
 //' @param lag An integer specifying the lag value for autocovariance.
 //' @param optbw_s Optional numeric vector specifying optimal bandwidths for \code{s}. Default is \code{NULL}.
 //' @param optbw_t Optional numeric vector specifying optimal bandwidths for \code{t}. Default is \code{NULL}.
 //' @param bw_grid Optional numeric vector of bandwidth grid values. Default is \code{NULL}.
 //' @param use_same_bw A logical value indicating if the same bandwidth should be used for \code{s} and \code{t}. Default is \code{FALSE}.
 //' @param center A logical value indicating if the data should be centered before estimation. Default is \code{TRUE}.
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
   arma::vec unique_id_curve = arma::unique(data_mat.col(0));
   double n_curve = unique_id_curve.n_elem;

   // Vector to be used in the estimation procedure
   arma::vec svec = s, tvec = t;
   if (lag == 0 && (optbw_s.isNull() || optbw_t.isNull())) {
     arma::mat mat_st_upper = get_upper_tri_couple(s, t);
     svec = mat_st_upper.col(0);
     tvec = mat_st_upper.col(1);
   }
   int n = tvec.n_elem;

   arma::vec optbw_s_to_use(n), optbw_t_to_use(n);
   arma::mat mat_locreg(n, 6);

   if (optbw_s.isNull() || optbw_t.isNull()) {
     arma::mat mat_risk = estimate_autocov_risk_cpp(data, svec, tvec, lag, bw_grid, use_same_bw, center, kernel_name);
     for (int k = 0; k < n; ++k) {
       arma::uvec idx_risk_cur = arma::find(mat_risk.col(0) == svec(k) && mat_risk.col(1) == tvec(k));
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

     if (optbw_s_to_use_temp.size() == n && optbw_t_to_use_temp.size() == n) {
       optbw_s_to_use = optbw_s_to_use_temp;
       optbw_t_to_use = optbw_t_to_use_temp;
     } else if (optbw_s_to_use_temp.size() == 1 && optbw_t_to_use_temp.size() == 1){
       optbw_s_to_use = arma::ones(n) * (optbw_s_to_use_temp(0));
       optbw_t_to_use = arma::ones(n) * (optbw_t_to_use_temp(0));
     } else {
       stop("If 'optbw_s' and 'optbw_t' are not NULL, they must be the same length as 's' and as 't' or of length 1.");
     }

     mat_locreg.col(0) = svec;
     mat_locreg.col(1) = tvec;
     mat_locreg.cols(2, 5).zeros();
   }

   // Estimate mean function
   arma::mat mat_mean_s = estimate_mean_cpp(data, svec, Rcpp::wrap(optbw_s_to_use), R_NilValue, kernel_name);
   arma::mat mat_mean_t = estimate_mean_cpp(data, tvec, Rcpp::wrap(optbw_t_to_use), R_NilValue, kernel_name);
   arma::vec muhat_s = mat_mean_s.col(5);
   arma::vec muhat_t = mat_mean_t.col(5);

   // Compute the autocovariance function
   arma::vec PN_lag(n, arma::fill::zeros);
   arma::vec autocov_numerator(n, arma::fill::zeros);

   for (int i = 0; i < n_curve - lag; ++i) {
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
     arma::vec pn_s = arma::regspace<arma::vec>(0, n - 1).transform([&](int j) { return arma::any(arma::abs(Tnvec_s - svec(j)) <= optbw_s_to_use(j)) ? 1.0 : 0.0; });
     arma::vec pn_t = arma::regspace<arma::vec>(0, n - 1).transform([&](int j) { return arma::any(arma::abs(Tnvec_t - tvec(j)) <= optbw_t_to_use(j)) ? 1.0 : 0.0; });
     PN_lag += pn_s % pn_t;

     // Estimate the numerator function numerator
     if (center) {
       autocov_numerator += pn_s % pn_t % (Xhat_s - muhat_s) % (Xhat_t - muhat_t);
     } else {
       autocov_numerator += pn_s % pn_t % Xhat_s % Xhat_t;
     }
   }

   // Compute autocovariance estimate
   arma::vec autocovhat;
   if (center) {
     autocovhat = autocov_numerator / PN_lag;
   } else {
     autocovhat = autocov_numerator / PN_lag - muhat_s % muhat_t;
   }

   // Initialize the result matrix
   int n_couple = t.size();
   arma::mat mat_res_autocov(n_couple, 14);

   // Diagonal correction if lag == 0
   if (lag == 0) {
     arma::mat mat_cov_res(n, 5);
     if (correct_diagonal) {
       arma::mat mat_cov(n, 5);
       mat_cov.col(0) = svec;
       mat_cov.col(1) = tvec;
       mat_cov.col(2) = arma::abs(svec - tvec);
       mat_cov.col(3) = optbw_s_to_use + optbw_t_to_use; // h_s h_t
       mat_cov.col(4) = autocovhat;
       Rcout << "Before matrix cov ---> \n";
       Rcout << mat_cov;
       // Diagonal correction
       arma::mat mat_cov_corrected = diagonal_correct(mat_cov, 0, 1, 2, 3, 4);
       Rcout << "Ok diagonal correction ---> \n";
       // replace the column containing arma::abs(svec - tvec) and max(h_s, h_t) by the bw_s and bw_t
       mat_cov_corrected.col(2) = optbw_s_to_use;
       mat_cov_corrected.col(3) = optbw_t_to_use;
       Rcout << "Ok optbw_t_to_use -> \n";
       mat_cov_res = mat_cov_corrected;
       Rcout << "Ok mat_cov_res -> \n";

     } else {
       mat_cov_res.col(0) = svec;
       mat_cov_res.col(1) = tvec;
       mat_cov_res.col(2) = optbw_s_to_use;
       mat_cov_res.col(3) = optbw_t_to_use;
       mat_cov_res.col(4) = autocovhat;
     }


     // Initialize the output matrix
     int n_upper_tri = mat_cov_res.n_rows;
     arma::mat res_cov(n_upper_tri * 2, 5);

     // Fill the upper triangular elements
     res_cov(arma::span(0, n_upper_tri - 1), arma::span::all) = mat_cov_res;

     // Fill the lower triangular elements
     res_cov(arma::span(n_upper_tri, 2 * n_upper_tri - 1), 0) = mat_cov_res.col(1);
     res_cov(arma::span(n_upper_tri, 2 * n_upper_tri - 1), 1) = mat_cov_res.col(0);
     res_cov(arma::span(n_upper_tri, 2 * n_upper_tri - 1), 2) = mat_cov_res.col(3);
     res_cov(arma::span(n_upper_tri, 2 * n_upper_tri - 1), 3) = mat_cov_res.col(2);
     res_cov(arma::span(n_upper_tri, 2 * n_upper_tri - 1), 4) = mat_cov_res.col(4);
     Rcout << "Ok res_cov -> \n";
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
       // Rcout << "Ok index ------> " << k << " - idx_cov_st = " << idx_cov_st << " \n";
       // Rcout << "Ok index ------> " << k << " - idx_locreg_st = " << idx_locreg_st << " \n";
       // Rcout << "Ok index ------> " << k << " - idx_PN_st = " << idx_PN_st << " \n";
       // Rcout << "Ok index ------> " << k << " - idx_muhat_st = " << idx_muhat_st << " \n";
       // Rcout << "Ok index ------> " << k << " - n_couple = " << n_couple << " \n";
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
     Rcout << "Ok for the loop ---> \n";
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




