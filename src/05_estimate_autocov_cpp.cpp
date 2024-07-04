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
         wvec_s.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);

         //// For the argument t
         arma::vec Tn_t_diff_over_bw = (data_mat.col(1) - t(k)) / bw;
         arma::vec wvec_t = kernel_func(Tn_t_diff_over_bw);
         wvec_t.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);

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


         // Init. output
         arma::vec pn_s_vec = arma::zeros(n_curve - lag);
         arma::vec pn_lag_t_vec = arma::zeros(n_curve - lag);
         double bias_term_num = 0;
         double variance_term_num = 0;

         for(int i = 0 ; i < n_curve - lag; ++i){
           // Exact the indexes : current and forward
           arma::uvec idx_i = arma::find(data_mat.col(0) == unique_id_curve(i));
           arma::uvec idx_i_lag = arma::find(data_mat.col(0) == unique_id_curve(i + lag));

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
           // // compute b_n(s;h_s) and b_n(t;h_t)
           arma::vec bn_s_vec = arma::pow(arma::abs(Tn_s_diff_over_bw(idx_i)), 2 * Hs) % arma::abs(wvec_s_i);
           arma::vec bn_t_vec = arma::pow(arma::abs(Tn_t_diff_over_bw(idx_i_lag)), 2 * Ht) % arma::abs(wvec_t_i_lag);

           // // Compute bias term numerator
           bias_term_num += mom_t_square * Ls * std::pow(bw, 2 * Hs) * pn_s_vec(i) * pn_lag_t_vec(i) * arma::sum(bn_s_vec) +
             mom_s_square * Lt * std::pow(bw, 2 * Ht) * pn_s_vec(i) * pn_lag_t_vec(i) * arma::sum(bn_t_vec);

           // Compute variance term numerator
           variance_term_num += sig_s_square * mom_t_square * pn_s_vec(i) * pn_lag_t_vec(i) * wvec_s_i.max() +
             sig_t_square * mom_s_square * pn_s_vec(i) * pn_lag_t_vec(i) * wvec_t_i_lag.max() +
             sig_s_square * sig_t_square * pn_s_vec(i) * pn_lag_t_vec(i) * wvec_s_i.max() * wvec_t_i_lag.max();
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
         double dependence_term = (dependence_term_num  / PN_lag) / 3;

         // Autocovariance risk
         double autocov_risk = bias_term + variance_term + dependence_term;
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
         double bw_t = bw_grid_to_use[idx_bw_t];

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
           wvec_s.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);

           //// For the argument t
           arma::vec Tn_t_diff_over_bw = (data_mat.col(1) - t(k)) / bw_t;
           arma::vec wvec_t = kernel_func(Tn_t_diff_over_bw);
           wvec_t.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);

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

           // Init. output
           arma::vec pn_s_vec = arma::zeros(n_curve - lag);
           arma::vec pn_lag_t_vec = arma::zeros(n_curve - lag);
           double bias_term_num = 0;
           double variance_term_num = 0;

           for(int i = 0 ; i < n_curve - lag; ++i){
             // Exact the indexes : current and forward
             arma::uvec idx_i = arma::find(data_mat.col(0) == unique_id_curve(i));
             arma::uvec idx_i_lag = arma::find(data_mat.col(0) == unique_id_curve(i + lag));

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
             // // compute b_n(s;h_s) and b_n(t;h_t)
             arma::vec bn_s_vec = arma::pow(arma::abs(Tn_s_diff_over_bw(idx_i)), 2 * Hs) % arma::abs(wvec_s_i);
             arma::vec bn_t_vec = arma::pow(arma::abs(Tn_t_diff_over_bw(idx_i_lag)), 2 * Ht) % arma::abs(wvec_t_i_lag);

             // // Compute bias term numerator
             bias_term_num += mom_t_square * Ls * std::pow(bw_s, 2 * Hs) * pn_s_vec(i) * pn_lag_t_vec(i) * arma::sum(bn_s_vec) +
               mom_s_square * Lt * std::pow(bw_t, 2 * Ht) * pn_s_vec(i) * pn_lag_t_vec(i) * arma::sum(bn_t_vec);

             // Compute variance term numerator
             variance_term_num += sig_s_square * mom_t_square * pn_lag_t_vec(i) * pn_lag_t_vec(i) * wvec_s_i.max() +
               sig_t_square * mom_s_square * pn_lag_t_vec(i) * pn_lag_t_vec(i) * wvec_t_i_lag.max() +
               sig_s_square * sig_t_square * pn_lag_t_vec(i) * pn_lag_t_vec(i) * wvec_s_i.max() * wvec_t_i_lag.max();
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
           double dependence_term = (dependence_term_num  / PN_lag) / 4;

           // Autocovariance risk
           double autocov_risk = bias_term + variance_term + dependence_term;
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
 arma::mat get_upper_tri_couple_old(const arma::vec& s, const arma::vec& t) {
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

 //' Remove Duplicate Rows from an Armadillo Matrix
 //'
 //' This function removes duplicate rows from an Armadillo matrix.
 //'
 //' @param mat An Armadillo matrix.
 //' @return A matrix with duplicate rows removed.
 //'
 //' @export
 // [[Rcpp::export]]
 arma::mat remove_duplicates(const arma::mat& mat) {
   // Step 1: Create a set of unique rows
   std::set<std::vector<double>> unique_rows_set;

   // Step 2: Populate the set with rows of mat
   for (arma::uword i = 0; i < mat.n_rows; ++i) {
     std::vector<double> row_vec(mat.row(i).begin(), mat.row(i).end());
     unique_rows_set.insert(row_vec);
   }

   // Step 3: Convert set back to Armadillo matrix
   arma::mat unique_mat(unique_rows_set.size(), mat.n_cols);
   arma::uword row_idx = 0;
   for (auto it = unique_rows_set.begin(); it != unique_rows_set.end(); ++it) {
     unique_mat.row(row_idx) = arma::rowvec(*it);
     ++row_idx;
   }

   return unique_mat;
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
     int optbw_s_to_use_temp_size = optbw_s_to_use_temp.size();
     int optbw_t_to_use_temp_size = optbw_t_to_use_temp.size();

     if (optbw_s_to_use_temp_size == n && optbw_t_to_use_temp_size == n) {
       optbw_s_to_use = optbw_s_to_use_temp;
       optbw_t_to_use = optbw_t_to_use_temp;
     } else if (optbw_s_to_use_temp_size == 1 && optbw_t_to_use_temp_size == 1){
       optbw_s_to_use = arma::ones(n) * (optbw_s_to_use_temp(0));
       optbw_t_to_use = arma::ones(n) * (optbw_t_to_use_temp(0));
     } else if (lag == 0) {
       arma::vec optbw_svec(n);
       arma::vec optbw_tvec(n);
       for (int j = 0; j < n; ++j) {
         arma::uvec idx_optbw_svecj = arma::find(s  == svec(j));
         arma::uvec idx_optbw_tvecj = arma::find(t  == tvec(j));
         optbw_svec(j) = optbw_s_to_use_temp(idx_optbw_svecj(0));
         optbw_tvec(j) = optbw_t_to_use_temp(idx_optbw_tvecj(0));
       }
       optbw_s_to_use = optbw_svec;
       optbw_t_to_use = optbw_tvec;
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

   // // Estimate the error sd for diagonal correction
   arma::vec vec_sig_s(n, arma::fill::zeros);
   arma::vec vec_sig_t(n, arma::fill::zeros);
   if (lag == 0 && correct_diagonal) {
     arma::mat mat_sig_s = estimate_sigma_cpp(data, svec);
     arma::mat mat_sig_t = estimate_sigma_cpp(data, tvec);
     vec_sig_s = mat_sig_s.col(1);
     vec_sig_t = mat_sig_t.col(1);
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

     // Diagonal correction
     if (lag == 0 && correct_diagonal) {
       // Compute the weight product
       arma::vec weight_product = arma::regspace<arma::vec>(0, n - 1).transform([&](int k) {
         // W_{n,k}(t;h_t)
         arma::vec weight_sk = kernel_func((Tnvec_s - svec(k)) / optbw_s_to_use(k));
         weight_sk /= arma::accu(weight_sk);
         weight_sk.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);

         // W_{n,k}(s;h_s)
         arma::vec weight_tk = kernel_func((Tnvec_t - tvec(k)) / optbw_t_to_use(k));
         weight_tk /= arma::accu(weight_tk);
         weight_tk.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
         return arma::accu(weight_sk % weight_tk);
         });

       // Estimate the numerator of the diagonal part
       diag_correct_numerator += vec_sig_s % vec_sig_t % pn_s % pn_t % weight_product;
     }
   }

   // Compute autocovariance estimate
   // // Note that if lag = 0 or correct_diagonal = FALSE, Then diag_correct_numerator = 0
   arma::vec autocovhat;
   if (center) {
     autocovhat = autocov_numerator / PN_lag - diag_correct_numerator / PN_lag;
   } else {
     autocovhat = autocov_numerator / PN_lag - muhat_s % muhat_t - diag_correct_numerator / PN_lag;
   }

   // Return the result
   // // Init output
   arma::mat mat_res_autocov(n, 14);

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
       arma::uvec idx_cov_st_upper = arma::find((res_autocov.col(0) == s(k)) && (res_autocov.col(1) == t(k)));
       if (! idx_cov_st_upper.is_empty()) {
         mat_res_autocov.row(k) = res_autocov.row(idx_cov_st_upper(0));
       } else {
         arma::uvec idx_cov_st_lower = arma::find((res_autocov.col(1) == s(k)) && (res_autocov.col(0) == t(k)));
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




