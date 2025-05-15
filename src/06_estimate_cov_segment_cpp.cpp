#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "02_smoothing_rcpp.h"
#include "03_estimate_locreg_rcpp.h"
#include "07_estimate_constants_cpp.h"
#include "04_estimate_mean_cpp.h"
using namespace Rcpp;
using namespace arma;

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
 //'
 //' @return A matrix containing the following ten columns:
 //'          \itemize{
 //'            \item{t : The points at which the risk function is estimated.}
 //'            \item{h : The candidate bandwidth.}
 //'            \item{PN : The number of curves used to estimate the mean at t. It corresponds to \eqn{P_N(t;h)}.}
 //'            \item{locreg_bw : The bandwidth used to estimate the local regularity parameters.}
 //'            \item{Ht : The estimates of the local exponent for each t. It corresponds to \eqn{H_t}}
 //'            \item{Lt : The estimates of the Hölder constant for each t. It corresponds to \eqn{L_t^2}.}
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
 //' @export
 //' @import Rdpack
 //' @references
 //' \insertAllCited{}

 //'
 // [[Rcpp::export]]
 arma::mat estimate_cov_segment_risk_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                         const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                         const bool center = true,
                                         const std::string kernel_name = "epanechnikov") {
   if (t.size() == 0) {
     Rcpp::stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.");
   }

   // Check if kernel_name is supported
   std::vector<std::string> supported_kernels = {"epanechnikov", "biweight", "triweight",
                                                 "tricube", "triangular", "uniform"};
   if (std::find(supported_kernels.begin(), supported_kernels.end(), kernel_name) == supported_kernels.end()) {
     Rcpp::stop("Unsupported kernel name. Choose from: epanechnikov, biweight, triweight, tricube, triangular, uniform.");
   }

   // Select the kernel function
   std::function<arma::vec(const arma::vec)> kernel_func = select_kernel(kernel_name);

   // Take the length of the vector t
   int n = t.size();

   // Convert data to matrix
   arma::mat data_mat(data.nrows(), 3);
   data_mat.col(0) = Rcpp::as<arma::vec>(data["id_curve"]);
   data_mat.col(1) = Rcpp::as<arma::vec>(data["tobs"]);
   data_mat.col(2) = Rcpp::as<arma::vec>(data["X"]);
   arma::vec unique_id_curve = arma::unique(data_mat.col(0));
   double n_curve = unique_id_curve.n_elem;

   // Set the bandwidth grid
   arma::vec bw_grid_to_use;
   if (bw_grid.isNull()) {
     // Estimate lambda
     double lambdahat = arma::mean(arma::hist(data_mat.col(0), unique_id_curve));
     double b0 = std::pow(n_curve * lambdahat, -1 / (2 * 0.05 + 1)); // rate with minimum local exponent = 0.05
     // double bK = 4 * std::pow(n_curve * lambdahat, -1 / (2 * 1 + 1)); // rate with maximum local exponent = 1
     double bK = 0.2; // rate with maximum local exponent = 1
     bw_grid_to_use = arma::logspace(log10(b0), log10(bK), 20);
   } else {
     bw_grid_to_use = Rcpp::as<arma::vec>(bw_grid);
   }
   int bw_size = bw_grid_to_use.size();

   // Estimate local regularity
   arma::mat mat_locreg = estimate_locreg_cpp(data, t, true, kernel_name, R_NilValue, R_NilValue);
   arma::vec h(n_curve, arma::fill::value(mat_locreg(0, 1))); // extract the presmoothing bandwidth

   // Estimate the error sd
   arma::mat mat_sig = estimate_sigma_cpp(data, t);
   arma::vec sig_vec = mat_sig.col(1);

   // Estimate moment
   arma::mat mat_mom_t = estimate_empirical_mom_cpp(data, arma::unique(t), h, 2, center, kernel_name);
   arma::vec mom_vec_t = mat_mom_t.col(2);

   // Estimate the autocovariance
   arma::mat mat_emp_autocov = estimate_empirical_XsXt_autocov_cpp(data, arma::unique(t), arma::unique(t), 0, arma::regspace(0, n_curve - 1), h, kernel_name, center);

   // Compute the numerator of \mathbb{D}(t;h_t)
   arma::mat mat_num_DD_t = estimate_numerator_dependence_term_DD_cpp(data, arma::unique(t), bw_grid_to_use, h, 3, kernel_name, center);

   // Compute the risk for each t and each bandwidth in bw_grid
   arma::mat mat_res_risk(bw_size * n, 10);

   // Parallelize outer loop
#pragma omp parallel for
   for (int idx_bw = 0; idx_bw < bw_size; ++idx_bw) {
     for (int k = 0; k < n; ++k) {
       // Extract local regularity parameters for each t
       arma::uvec idx_locreg_cur = arma::find(mat_locreg.col(0) == t(k));
       double Ht = mat_locreg(idx_locreg_cur(0), 4);
       double Lt = mat_locreg(idx_locreg_cur(0), 5);

       // Compute the weight vectors for each t and for each bw
       // and replace replace non-finite values with 0
       arma::vec Tn_t_diff_over_bw = (data_mat.col(1) - t(k)) / bw_grid_to_use(idx_bw);
       arma::vec wvec = kernel_func(Tn_t_diff_over_bw);
       wvec.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);

       // Extract the sd
       arma::uvec cur_idx_sig = arma::find(mat_sig.col(0) == t(k));
       double sig_square = sig_vec(cur_idx_sig(0)) * sig_vec(cur_idx_sig(0));

       // Extract moment
       arma::uvec cur_idx_mom_t = arma::find(mat_mom_t.col(0) == t(k));
       double mom_t_square = mom_vec_t(cur_idx_mom_t(0));

       // Init output
       arma::vec pn_vec(n_curve);
       double bias_term_num = 0;
       double variance_term_num = 0;

       for (int i = 0; i < n_curve; ++i) {
         // Extract the current curve index data
         arma::uvec idx_i = arma::find(data_mat.col(0) == i + 1);

         // Extract weight
         arma::vec wvec_i = wvec(idx_i) / arma::accu(wvec(idx_i));
         wvec_i.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);

         // Compute and store the vector \pi_n(t;h)
         pn_vec(i) = arma::find(abs(Tn_t_diff_over_bw(idx_i)) <= 1).is_empty() ? 0 : 1;

         // Compute bias term numerator
         // // compute and b_n(t;h)
         arma::vec bn_vec = arma::pow(arma::abs(Tn_t_diff_over_bw(idx_i)), 2 * Ht) % arma::abs(wvec_i);

         // // Compute bias term numerator
         bias_term_num += Lt * std::pow(bw_grid_to_use(idx_bw), 2 * Ht) * pn_vec(i) * arma::sum(bn_vec);

         // Compute variance term numerator
         variance_term_num += sig_square * pn_vec(i) * wvec_i.max();
       }

       // Compute P_N(t;h)
       double PN = arma::accu(pn_vec);

       // Compute bias term
       double bias_term = 20 * mom_t_square * bias_term_num / PN;

       // Compute variance term
       double variance_term = 20 * mom_t_square * variance_term_num / (PN * PN);


       // ::::::::: Compute dependence term :::::::
       // Extract the numerator of \mathbb{D}(t;h_t)
       arma::uvec idx_num_DD_t = arma::find( (mat_num_DD_t.col(0) == t(k)) % (mat_num_DD_t.col(1) == bw_grid_to_use(idx_bw)) );
       double num_Dt = mat_num_DD_t(idx_num_DD_t(0), 2);
       double PNt = mat_num_DD_t(idx_num_DD_t(0), 3);
       double Dt = num_Dt / std::pow(PNt, 3) ;
       double first_dependence_term = 25 * Dt / PN;

       // Compute \mathbb{D}(s;h_s)
       // Add the lag-0 autocovariance to the vector
       arma::uvec idx_lag0 = arma::find( (mat_emp_autocov.col(0) == t(k)) % (mat_emp_autocov.col(1) == t(k)) % (mat_emp_autocov.col(3) == 0) );
       arma::uvec idx_lag = arma::find( (mat_emp_autocov.col(0) == t(k)) % (mat_emp_autocov.col(1) == t(k)) % (mat_emp_autocov.col(3) != 0) );

       double XsXt_var = mat_emp_autocov(idx_lag0(0), 5);
       arma::mat XsXt_mat_lr_var = mat_emp_autocov.rows(idx_lag);
       double second_dependence_term_num = XsXt_var + arma::accu(arma::abs(2 * XsXt_mat_lr_var.col(5))) / PN;
       double second_dependence_term = second_dependence_term_num / PN;

       double dependence_term = first_dependence_term + second_dependence_term;
       // ::::::::::: End dependence term ::::::::

       // Covariance segment risk function
       double cov_segment_risk = bias_term + variance_term + dependence_term;

#pragma omp critical
{
  mat_res_risk.row(idx_bw * n + k) = {t(k), bw_grid_to_use(idx_bw), PN, h(0), Ht, Lt, bias_term, variance_term, dependence_term, cov_segment_risk};
}
     }
   }
   return mat_res_risk;
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
 //' @param optbw An optional numeric vector of bandwidths to be used at each evaluation point in \code{t}.
 //' If \code{NULL}, bandwidths are selected based on risk minimization.
 //' @param bw_grid An optional numeric vector of candidate bandwidths used for optimization when \code{optbw} is \code{NULL}.
 //' @param center A logical value indicating if the data should be centered before estimation. Default is \code{true}.
 //' @param kernel_name A string indicating the kernel to use. Supported kernels are:
 //' \code{"epanechnikov"}, \code{"biweight"}, \code{"triweight"}, \code{"tricube"},
 //' \code{"triangular"}, and \code{"uniform"}.
 //'
 //' @return A matrix with 8 columns and \code{length(t)} rows:
 //' \describe{
 //'   \item{\code{t}}{Evaluation point}
 //'   \item{\code{optbw}}{Bandwidth used at \code{t}}
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
 //' @export
 //' @import Rdpack
 //' @references
 //' \insertAllCited{}
 //'
 // [[Rcpp::export]]
 arma::mat estimate_cov_segment_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                    const Rcpp::Nullable<arma::vec> optbw = R_NilValue,
                                    const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                    const bool center = true,
                                    const std::string kernel_name = "epanechnikov"){
   if (t.size() == 0) {
     stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.");
   }
   // Check if kernel_name is one of the supported kernels
   std::vector<std::string> supported_kernels = {"epanechnikov", "biweight", "triweight",
                                                 "tricube", "triangular", "uniform"};
   if (std::find(supported_kernels.begin(), supported_kernels.end(), kernel_name) == supported_kernels.end()) {
     stop("Unsupported kernel name. Choose from: epanechnikov, biweight, triweight, tricube, triangular, uniform.");
   }

   // Select the kernel function
   std::function<arma::vec(const arma::vec)> kernel_func = select_kernel(kernel_name);

   // Take the length of the vector t
   int n = t.size();

   // Convert data to matrix
   arma::mat data_mat(data.nrows(), 3);
   data_mat.col(0) = as<arma::vec>(data["id_curve"]);
   data_mat.col(1) = as<arma::vec>(data["tobs"]);
   data_mat.col(2) = as<arma::vec>(data["X"]);
   arma::vec unique_id_curve = arma::unique(data_mat.col(0));
   double n_curve = unique_id_curve.n_elem;

   // Estimate the error sd
   arma::mat mat_sig = estimate_sigma_cpp(data, t);
   arma::vec sig_square_vec = arma::square(mat_sig.col(1));

   // Estimate / Verify optimal risk function
   arma::vec optbw_to_use(n);
   arma::vec Ht_used(n);
   arma::vec Lt_used(n);
   if (optbw.isNull()) {
     arma::mat mat_risk = estimate_cov_segment_risk_cpp(data, t, bw_grid, center, kernel_name);
     for (int k = 0; k < n; ++k) {
       // Find rows in mat_risk where the first column equals t(k)
       arma::uvec idx_risk_cur = arma::find(mat_risk.col(0) == t(k));

       // Extract the risk column for those rows
       arma::vec risk = mat_risk(idx_risk_cur, arma::uvec({9}));

       // Find the minimum index in the risk column, ignoring non-finite values
       arma::uword idx_min = arma::index_min(risk.elem(arma::find_finite(risk)));

       // Extract values corresponding to the minimum risk index
       optbw_to_use(k) = mat_risk(idx_risk_cur(idx_min), 1);
       Ht_used(k) = mat_risk(idx_risk_cur(idx_min), 4);
       Lt_used(k) = mat_risk(idx_risk_cur(idx_min), 5);
     }
   } else {
     arma::vec optbw_cur = as<arma::vec>(optbw);
     int optbw_cur_size = optbw_cur.size();
     if (optbw_cur_size != n) {
       stop("If 'optbw' is not NULL, it must be the same length as 't'.");
     } else {
       optbw_to_use = optbw_cur;
     }
   }

   // Compute the mean function
   arma::mat mat_res_mean = estimate_mean_cpp(data, t, Rcpp::wrap(optbw_to_use), R_NilValue, kernel_name);
   arma::vec muhat = mat_res_mean.col(5);

   // Compute the covariance segment function
   // Init output
   arma::mat mat_res_covseg(n, 8);
   arma::vec PN(n, fill::zeros);
   arma::vec covseg_numerator(n, fill::zeros);
   arma::vec covseg_correction_numerator(n, fill::zeros);

   for(int i = 0 ; i < n_curve; ++i){
     // Extrat the current curve index data
     double idx_cur_curve = unique_id_curve(i);
     arma::uvec indices_cur = arma::find(data_mat.col(0) == idx_cur_curve);
     arma::vec Tnvec = data_mat(indices_cur, arma::uvec({1}));
     arma::vec Ynvec = data_mat(indices_cur, arma::uvec({2}));

     // Smooth using Nadaraya-Watson estimator
     arma::vec Xhat = estimate_nw_cpp(Ynvec, Tnvec, t, optbw_to_use, kernel_name);
     Xhat.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);

     // Compute the vector p_n(t;h) and update P_N(t;h)
     arma::vec pn = arma::regspace(0, n - 1);
     pn.transform([t, optbw_to_use, Tnvec](int j) { return  arma::find(abs((Tnvec - t(j))) <= optbw_to_use(j)).is_empty() ? 0 : 1 ;});
     PN += pn;

     // Estimate the mean function numerator
     covseg_numerator += pn % arma::square(Xhat - muhat);

     // The correction term
     arma::vec sum_weight_square(n);
     for(int k = 0 ; k < n; ++k) {
       // Compute the weight vectors for each t and for each optbw
       // and replace replace non-finite values with 0
       arma::vec Tn_t_diff_over_bw = (Tnvec - t(k)) / optbw_to_use(k);
       arma::vec num_wvec = kernel_func(Tn_t_diff_over_bw);
       num_wvec.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
       arma::vec wvec = num_wvec / arma::accu(num_wvec);
       wvec.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
       sum_weight_square(k) = arma::accu(arma::square(wvec));
     }

     covseg_correction_numerator += pn % sum_weight_square;


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


