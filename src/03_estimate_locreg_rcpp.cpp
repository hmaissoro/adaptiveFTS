#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "02_smoothing_rcpp.h"
using namespace Rcpp;
using namespace arma;

//' Local Regularity Parameters Estimation
 //'
 //' @param data A DataFrame containing the columns "id_curve", "tobs", and "X". Typically, the output of the function \link{format_data}.
 //' @param t Numeric vector. Observation points at which we want to estimate the local regularity parameters of the underlying process.
 //' @param center Logical. If \code{TRUE}, the curves are centered.
 //' @param kernel_name String. The kernel function of the Nadaraya-Watson estimator. Default is "epanechnikov".
 //' @param h Numeric (positive vector or scalar). The bandwidth of the Nadaraya-Watson estimator for the local regularity estimation.
 //' Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
 //' If \code{h} is a scalar, then all curves will be smoothed with the same bandwidth.
 //' Otherwise, if \code{h} is a vector, its length must be equal to the number of curves in \code{data}
 //' and each element of the vector must correspond to a curve given in the same order as in \code{data}.
 //' @param Delta Numeric (positive). The length of the neighborhood of \code{t} around which the local regularity is to be estimated.
 //' Default \code{Delta = NULL} and thus it will be estimated from the data.
 //'
 //' @return A \code{matrix} containing the following six columns in order:
 //' \enumerate{
 //'   \item t: The points around which the local regularity parameters are estimated.
 //'   \item locreg_bw: The presmoothing bandwidth.
 //'   \item Delta: The length of the neighborhood of \code{t} around which the local regularity is to be estimated.
 //'   \item Nused: The number of curves that give non-degenerate estimates around \code{t}.
 //'   \item Ht: The local exponent estimates for each \code{t}.
 //'   \item Lt: The Hölder constant estimates for each \code{t}.
 //' }
 //'
 //' @examples
 //' \dontrun{
 //' data <- data.frame(id_curve = rep(1:5, each = 10), tobs = rep(seq(0, 1, length.out = 10), 5), X = rnorm(50))
 //' estimate_locreg_cpp(data, seq(0, 1, length.out = 20), center = TRUE)
 //' }
 //'
 //' @export
 // [[Rcpp::export]]
 arma::mat estimate_locreg_cpp(const Rcpp::DataFrame data, const arma::vec t,
                               const bool center,
                               const std::string kernel_name = "epanechnikov",
                               const Rcpp::Nullable<arma::vec> h = R_NilValue,
                               const Rcpp::Nullable<double> Delta = R_NilValue){
   if (t.size() == 0) {
     stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.");
   }
   // Check if kernel_name is one of the supported kernels
   std::vector<std::string> supported_kernels = {"epanechnikov", "biweight", "triweight",
                                                 "tricube", "triangular", "uniform"};
   if (std::find(supported_kernels.begin(), supported_kernels.end(), kernel_name) == supported_kernels.end()) {
     stop("Unsupported kernel name. Choose from: epanechnikov, biweight, triweight, tricube, triangular, uniform.");
   }

   // Smooth curve
   arma::mat data_mat(data.nrows(), 3);
   data_mat.col(0) = as<arma::vec>(data["id_curve"]);
   data_mat.col(1) = as<arma::vec>(data["tobs"]);
   data_mat.col(2) = as<arma::vec>(data["X"]);
   arma::vec unique_id_curve = arma::unique(data_mat.col(0));
   int n_curve = unique_id_curve.n_elem;

   // Control the presmoothing bandwidth
   arma::vec hvec_to_use(n_curve);
   if (h.isNull()) {
     // Get the bandwidth grid and estimate best bandwidth
     double hbest = get_nw_optimal_bw_cpp(data, R_NilValue, R_NilValue, kernel_name);
     hvec_to_use.fill(hbest);
   } else {
     arma::vec hvec = Rcpp::as<arma::vec>(h);
     if (hvec.size() != n_curve) {
       if (hvec.size() != 1) {
         Rcpp::stop("'h' must be a numeric vector with length equal to the number of curves in 'data' or a scalar value.");
       } else {
         hvec_to_use.fill(hvec(0));
       }
     } else {
       hvec_to_use = hvec;
     }
   }

   // Control on Delta
   double Delta_to_use = 0;
   if (Delta.isNull() || Delta_to_use < 0 || Delta_to_use > 1) {
     // Estimate lambda
     double lambdahat = arma::mean(hist(data_mat.col(0), unique_id_curve));
     Delta_to_use = std::min(exp(- std::pow(log(lambdahat), 1.0 / 3)), 0.2);
   } else {
     Delta_to_use = as<double>(Delta);
   }

   // Take into account the cases where t-Delta/2 < 0 and t + Delta/2
   int n = t.size();
   arma::vec t1(n), t2(n), t3(n);

   for (int i = 0; i < n; ++i) {
     double ti = t(i);

     // Add small uniform value to avoid ties
     if ((ti - Delta_to_use / 2) <= 0) {
       t1(i) = ti;
       t2(i) = ti + Delta_to_use / 2 + arma::randu(distr_param(1e-06, 2e-06));
       t3(i) = ti + Delta_to_use + arma::randu(distr_param(1e-06, 2e-06));
     } else if ((ti + Delta_to_use / 2) >= 1) {
       t3(i) = ti;
       t2(i) = ti - Delta_to_use / 2 + arma::randu(distr_param(1e-06, 2e-06));;
       t1(i) = ti - Delta_to_use + arma::randu(distr_param(1e-06, 2e-06));;
     } else {
       t1(i) = ti - Delta_to_use / 2 + arma::randu(distr_param(1e-06, 2e-06));;
       t2(i) = ti;
       t3(i) = ti + Delta_to_use / 2 + arma::randu(distr_param(1e-06, 2e-06));;
     }
   }

   // Initialize a matrix to put
   arma::mat mat_res_nw(n_curve * n, 7);

   for(int i = 0 ; i < n_curve; ++i){
     // Extrat the current curve index data
     double idx_cur_curve = unique_id_curve(i);
     arma::uvec indices_cur = arma::find(data_mat.col(0) == idx_cur_curve);
     arma::mat mat_cur = data_mat.rows(indices_cur);

     // Smooth using Nadaraya-Watson estimator
     arma::vec h_to_use(1, fill::value(hvec_to_use(idx_cur_curve - 1)));
     arma::vec Xhat_t1 = estimate_nw_cpp(mat_cur.col(2), mat_cur.col(1), t1, h_to_use, kernel_name);
     arma::vec Xhat_t2 = estimate_nw_cpp(mat_cur.col(2), mat_cur.col(1), t2, h_to_use, kernel_name);
     arma::vec Xhat_t3 = estimate_nw_cpp(mat_cur.col(2), mat_cur.col(1), t3, h_to_use, kernel_name);

     // Store the data
     arma::mat cur_mat(n, 7);
     cur_mat.col(0) = idx_cur_curve * arma::ones(t.n_elem);
     cur_mat.col(1) = t1;
     cur_mat.col(2) = Xhat_t1;
     cur_mat.col(3) = t2;
     cur_mat.col(4) = Xhat_t2;
     cur_mat.col(5) = t3;
     cur_mat.col(6) = Xhat_t3;
     mat_res_nw(span(i * n, (i + 1) * n - 1), span(0, 6)) = cur_mat;
   }

   // Estimate locreg
   arma::mat mat_locreg(n, 6);
   for(int j = 0 ; j < n; ++j){
     // Extract all smooth data for each t
     arma::uvec cur_idx = arma::find((mat_res_nw.col(1) == t1(j)) % (mat_res_nw.col(3) == t2(j)) % (mat_res_nw.col(5) == t3(j)));
     arma::mat mat_res_nw_by_t = mat_res_nw.rows(cur_idx);
     arma::vec xt1 = mat_res_nw_by_t.col(2);
     arma::vec xt2 = mat_res_nw_by_t.col(4);
     arma::vec xt3 = mat_res_nw_by_t.col(6);

     // Remove NaN values
     arma::uvec finiteIndices = arma::find_finite(xt1 + xt2 + xt3);
     arma::vec filtered_xt1 = xt1.elem(finiteIndices);
     arma::vec filtered_xt2 = xt2.elem(finiteIndices);
     arma::vec filtered_xt3 = xt3.elem(finiteIndices);

     // Remove extreme values
     arma::vec p0025(finiteIndices.size(), fill::value(0.025));
     arma::vec p0975(finiteIndices.size(), fill::value(0.975));
     arma::uvec rxt1 = (filtered_xt1 >= arma::quantile(filtered_xt1, p0025)) % (filtered_xt1 <= arma::quantile(filtered_xt1, p0975));
     arma::uvec rxt2 = (filtered_xt2 >= arma::quantile(filtered_xt2, p0025)) % (filtered_xt2 <= arma::quantile(filtered_xt2, p0975));
     arma::uvec rxt3 = (filtered_xt3 >= arma::quantile(filtered_xt3, p0025)) % (filtered_xt3 <= arma::quantile(filtered_xt3, p0975));

     arma::vec xt1_final = filtered_xt1.elem(arma::find(rxt1 % rxt2 % rxt3));
     arma::vec xt2_final = filtered_xt2.elem(arma::find(rxt1 % rxt2 % rxt3));
     arma::vec xt3_final = filtered_xt3.elem(arma::find(rxt1 % rxt2 % rxt3));

     // Center data if necessary
     if (center) {
       xt1_final -= arma::mean(xt1_final);
       xt2_final -= arma::mean(xt2_final);
       xt3_final -= arma::mean(xt3_final);
     }

     // Compute local regularity parameters
     double theta_t1_t3 = arma::mean(arma::square(xt1_final - xt3_final));
     double theta_t1_t2 = arma::mean(arma::square(xt1_final - xt2_final));
     double theta_t2_t3 = arma::mean(arma::square(xt2_final - xt3_final));
     double Ht_row = (std::log(theta_t1_t3) - std::log(theta_t2_t3)) / (2 * std::log(2));
     double Ht = (Ht_row <= 0.1) ? 0.1 : (Ht_row >=1) ? 1 : Ht_row;
     double Lt = theta_t2_t3 / std::pow(std::abs(t2(j) - t3(j)), 2 * Ht);

     // Store the result in matrix mat_locreg initialized above
     double Nused = xt1_final.size();
     double locreg_bw = arma::median(hvec_to_use);
     mat_locreg.row(j) = {t(j), locreg_bw, Delta_to_use, Nused, Ht, Lt};
   }

   return mat_locreg;
 }
