#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "02_smoothing_rcpp.h"
#include "03_estimate_locreg_rcpp.h"
#include "07_estimate_constants_cpp.h"
using namespace Rcpp;
using namespace arma;

//' Estimate the risk of the mean function
 //'
 //' This function estimates the risk function \eqn{R_\mu(t;h)} of the mean function estimation proposed by Maissoro et al. (2024).
 //'
 //' @param data A DataFrame containing the columns "id_curve", "tobs", and "X", typically the output of the function \link{.format_data}.
 //' @param t A numeric vector of observation points at which we want to estimate the mean function of the underlying process.
 //' @param bw_grid A numeric vector representing the bandwidth grid in which the best smoothing parameter is selected for each t.
 //' It can be NULL, in which case it will be defined as an exponential grid of \eqn{N \times \lambda}.
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
 //' @export
 //'
 // [[Rcpp::export]]
 arma::mat estimate_mean_risk_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                  const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
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

   // Set the bandwidth grid
   arma::vec bw_grid_to_use;
   if (bw_grid.isNull()) {
     // Estimate lambda
     double lambdahat = arma::mean(hist(data_mat.col(0), unique_id_curve));
     double b0 =  std::pow(n_curve * lambdahat, - 1 / (2 * 0.05 + 1)); // rate with minimum local exponent = 0.05
     // double bK = 4 * std::pow(n_curve * lambdahat, -1 / (2 * 1 + 1)); // rate with maximum local exponent = 1
     double bK = 0.5; // rate with maximum local exponent = 1
     bw_grid_to_use = arma::logspace(log10(b0), log10(bK), 30);
   } else {
     bw_grid_to_use = as<arma::vec>(bw_grid);
   }
   int bw_size = bw_grid_to_use.size();

   // Estimate local regularity
   arma::mat mat_locreg = estimate_locreg_cpp(data, t, true, kernel_name, R_NilValue, R_NilValue);
   arma::vec h(n_curve, fill::value(mat_locreg(0, 1))); // extract the presmoothing bandwidth

   // Estimate the error sd
   arma::mat mat_sig = estimate_sigma_cpp(data, t);
   arma::vec sig_vec = mat_sig.col(1);

   // Estimate the autocovariance
   arma::mat mat_emp_autocov = estimate_empirical_autocov_cpp(data, t, h, arma::regspace(0, n_curve-1), kernel_name);

   // Compute the risk for each t and each bandwidth in bw_grid
   arma::mat mat_res_risk(bw_size * n, 10);
   for (int idx_bw = 0; idx_bw < bw_size; ++idx_bw) {
     for (int k = 0; k < n; ++k) {
       // Extract local regularity parameters for each t
       arma::uvec idx_locreg_cur = arma::find(mat_locreg.col(0) == t(k));
       double Ht = mat_locreg(idx_locreg_cur(0), 4);
       double Lt = mat_locreg(idx_locreg_cur(0), 5);
       // Init output
       arma::vec pn_vec(n_curve);
       double bias_term_num = 0;
       double variance_term_num = 0;

       for(int i = 0 ; i < n_curve; ++i){
         // Extrat the current curve index data
         double idx_cur_curve = unique_id_curve(i);
         arma::uvec indices_cur = arma::find(data_mat.col(0) == idx_cur_curve);
         arma::vec Tnvec = data_mat(indices_cur, arma::uvec({1}));
         int Mn = Tnvec.size();

         // Compute the weight vector for each t and for each bw
         // At the we replace replace non-finite values with 0
         arma::vec wvec = arma::zeros(Mn);
         arma::vec Tn_t_diff_over_bw = (Tnvec - t(k)) / bw_grid_to_use(idx_bw);
         wvec = kernel_func(Tn_t_diff_over_bw);
         wvec /= arma::accu(wvec);
         wvec.replace(arma::datum::nan, 0);
         wvec.replace(arma::datum::inf, 0);
         wvec.replace(-arma::datum::inf, 0);

         // Compute the maximum and c_n(t;h)
         double wmax = wvec.max();
         double cn = arma::accu(wvec);

         // Compute the vector p_n(t;h)
         arma::uvec idx_is_one_pi = arma::find(abs(Tn_t_diff_over_bw) <= 1);
         double pn = 0;
         if (!idx_is_one_pi.is_empty()) {
           pn = 1;
         }

         // Store \pi_n(t;h)
         pn_vec(i) = pn;

         // Compute bias term numerator
         arma::vec Tn_t_2H = arma::pow(arma::abs(Tn_t_diff_over_bw), 2 * Ht);
         double bn = arma::sum(Tn_t_2H % arma::abs(wvec));
         bias_term_num += Lt * std::pow(bw_grid_to_use(idx_bw), 2 * Ht) * pn * cn * bn;

         // Compute variance term numerator
         arma::uvec cur_idx_sig = arma::find(mat_sig.col(0) == t(k));
         double sig_square = sig_vec(cur_idx_sig(0)) * sig_vec(cur_idx_sig(0));
         variance_term_num += sig_square * pn * cn * wmax;
       }

       // Compute P_N(t;h)
       double PN = arma::accu(pn_vec);

       // Compute bias term
       double bias_term = 2 * bias_term_num / PN;

       // Compute variance term
       double variance_term = 2 * variance_term_num / (PN * PN);

       // Compute dependence term
       arma::vec dep_term_all_lag(n_curve);
       // Add the lag-0 autocovariance to the vector
       arma::uvec idx_lag0 = arma::find(mat_emp_autocov.col(0) == t(k) && mat_emp_autocov.col(1) == 0);
       dep_term_all_lag(0) = mat_emp_autocov(idx_lag0(0), 2);
       for (int ell = 1; ell < n_curve; ++ell) {
         arma::vec pi_i = pn_vec.subvec(0, n_curve - 1 - ell);
         arma::vec pi_i_plus_ell = pn_vec.subvec(ell, n_curve - 1);
         double rho = sum(pi_i % pi_i_plus_ell) / PN;
         arma::uvec idx_lag = arma::find(mat_emp_autocov.col(0) == t(k) && mat_emp_autocov.col(1) == ell);
         double autocov = mat_emp_autocov(idx_lag(0), 2);
         dep_term_all_lag(ell) =  2 * autocov * rho;
       }
       double dependence_term = 2 * std::abs(arma::accu(dep_term_all_lag))  / PN;
       double mean_risk = bias_term + variance_term + dependence_term;
       mat_res_risk.row(idx_bw * n + k) = {t(k), bw_grid_to_use(idx_bw), PN, h(0), Ht, Lt, bias_term, variance_term, dependence_term, mean_risk};
     }
   }
   return mat_res_risk;
 }

 //' Estimate Mean Function for Functional Data
 //'
 //' This function estimates the mean function for functional data using the
 //' Nadaraya-Watson estimator with a specified kernel. The function takes a
 //' data frame of observations, a vector of evaluation points, and optional
 //' bandwidth parameters.
 //'
 //' @param data A data frame containing the functional data. It should have three columns:
 //' \code{id_curve} (curve identifiers), \code{tobs} (observation times), and \code{X} (observations).
 //' @param t A numeric vector of evaluation points between 0 and 1.
 //' @param optbw An optional numeric vector of optimal bandwidths. If \code{NULL}, the function
 //' will estimate the optimal bandwidths.
 //' @param bw_grid An optional numeric vector of bandwidths to be used if \code{optbw} is \code{NULL}.
 //' @param kernel_name A string specifying the kernel to be used. Supported kernels are:
 //' \code{"epanechnikov"}, \code{"biweight"}, \code{"triweight"}, \code{"tricube"},
 //' \code{"triangular"}, and \code{"uniform"}.
 //' @return A matrix where each row contains the following columns:
 //' \itemize{
 //'   \item \code{t} - The evaluation point.
 //'   \item \code{optbw} - The optimal bandwidth used.
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
 // [[Rcpp::export]]
 arma::mat estimate_mean_cpp(const Rcpp::DataFrame data, const arma::vec t,
                             const Rcpp::Nullable<arma::vec> optbw = R_NilValue,
                             const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
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

   // Estimate / Verify optimal risk function
   arma::vec optbw_to_use(n);
   arma::vec Ht_used(n);
   arma::vec Lt_used(n);
   if (optbw.isNull()) {
     arma::mat mat_risk = estimate_mean_risk_cpp(data, t, bw_grid, kernel_name);
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
     if (optbw_cur.size() != n) {
       stop("If 'optbw' is not NULL, it must be the same length as 't'.");
     } else {
       optbw_to_use = optbw_cur;
     }
   }

   // Compute the mean function
   // Init output
   arma::mat mat_res_mean(n, 6);
   arma::vec PN(n, fill::zeros);
   arma::vec mean_numerator(n, fill::zeros);

   for(int i = 0 ; i < n_curve; ++i){
     // Extrat the current curve index data
     double idx_cur_curve = unique_id_curve(i);
     arma::uvec indices_cur = arma::find(data_mat.col(0) == idx_cur_curve);
     arma::vec Tnvec = data_mat(indices_cur, arma::uvec({1}));
     arma::vec Ynvec = data_mat(indices_cur, arma::uvec({2}));
     int Mn = Tnvec.size();

     // Smooth using Nadaraya-Watson estimator
     arma::vec Xhat = estimate_nw_cpp(Ynvec, Tnvec, t, optbw_to_use, kernel_name);
     Xhat.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);

     // Compute the vector p_n(t;h) and update P_N(t;h)
     arma::vec pn = arma::regspace(0, n - 1);
     pn.transform([t, optbw_to_use, Tnvec](int j) { return  arma::find(abs((Tnvec - t(j))) <= optbw_to_use(j)).is_empty() ? 0 : 1 ;});
     PN += pn;

     // Estimate the mean function numerator
     mean_numerator += pn % Xhat;

   }

   // Compute \widehat \mu_N(t;h)
   arma::vec muhat = mean_numerator / PN;

   // return the result
   mat_res_mean.col(0) = t;
   mat_res_mean.col(1) = optbw_to_use;
   mat_res_mean.col(2) = Ht_used;
   mat_res_mean.col(3) = Lt_used;
   mat_res_mean.col(4) = PN;
   mat_res_mean.col(5) = muhat;
   return mat_res_mean;
 }


