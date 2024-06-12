#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "02_smoothing_rcpp.h"
#include "03_estimate_locreg_rcpp.h"
#include "07_estimate_constants_cpp.h"
#include "04_estimate_mean_cpp.h"
#include "05_estimate_autocov_cpp.h"
using namespace Rcpp;
using namespace arma;


//' Reshape Matrix from Long Format Data
 //'
 //'
 //' This function reshapes a matrix from long format data where each row
 //' represents a value corresponding to specific (s, t) coordinates.
 //'
 //' @param A A matrix with at least three columns. The columns represent
 //'          the s coordinates, t coordinates, and the values, respectively.
 //' @param idx_col_s The index of the column in `A` that contains the s coordinates.
 //' @param idx_col_t The index of the column in `A` that contains the t coordinates.
 //' @param idx_col_value The index of the column in `A` that contains the values.
 //'
 //' @details
 //' The function first sorts the input matrix `A` based on the s and t coordinates.
 //' It then extracts the unique s and t coordinates to determine the dimensions of
 //' the reshaped matrix. Finally, it populates the reshaped matrix by
 //' matching the s and t coordinates to the corresponding values.
 //'
 //' @return
 //' A reshaped matrix where each element at position (i, j) represents the value
 //' at the corresponding (s, t) coordinates.
 //'
 //' @examples
 //' \dontrun{
 //' library(Rcpp)
 //' library(RcppArmadillo)
 //'
 //' # Source the C++ code
 //' sourceCpp("path/to/your/file.cpp")
 //'
 //' # Example data matrix
 //' A <- matrix(c(1, 1, 10,
 //'               1, 2, 20,
 //'               2, 1, 30,
 //'               2, 2, 40,
 //'               3, 1, 50,
 //'               3, 2, 60), ncol = 3, byrow = TRUE)
 //'
 //' # Reshape to wide format
 //' wide_matrix <- reshape_matrix(A, 0, 1, 2)
 //' print(wide_matrix)
 //' }
 // [[Rcpp::export]]
 arma::mat reshape_matrix(const arma::mat& A,
                          const arma::uword idx_col_s,
                          const arma::uword idx_col_t,
                          const arma::uword idx_col_value) {
   if (A.n_cols < 3) {
     Rcpp::stop("Input matrix A must have at least 3 columns.");
   }

   // Unique s and t values
   arma::vec svec = arma::unique(A.col(idx_col_s));
   arma::vec tvec = arma::unique(A.col(idx_col_t));
   int ns = svec.size();
   int nt = tvec.size();

   // Initialize the reshaped matrix
   arma::mat reshaped(ns, nt, arma::fill::zeros);

   // Populate the reshaped matrix
   for (arma::uword i = 0; i < A.n_rows; ++i) {
     arma::uword row_idx = arma::index_min(arma::abs(svec - A(i, idx_col_s)));
     arma::uword col_idx = arma::index_min(arma::abs(tvec - A(i, idx_col_t)));
     reshaped(row_idx, col_idx) = A(i, idx_col_value);
   }

   return reshaped;
 }



 //' Combine Four Matrices into One Block Matrix
 //'
 //' This function combines four matrices into one block matrix.
 //'
 //' @param A First matrix.
 //' @param B Second matrix.
 //' @param C Third matrix.
 //' @param D Fourth matrix.
 //' @return A block matrix combining A, B, C, and D.
 //'
 //' @examples
 //' \dontrun{
 //' A <- matrix(1:4, 2, 2)
 //' B <- matrix(5:8, 2, 2)
 //' C <- matrix(9:12, 2, 2)
 //' D <- matrix(13:16, 2, 2)
 //' combine_matrices(A, B, C, D)
 //' }
 //' @export
 // [[Rcpp::export]]
 arma::mat combine_matrices(const arma::mat& A, const arma::mat& B, const arma::mat& C, const arma::mat& D) {
   // Check if the dimensions match for horizontal and vertical concatenation
   if (A.n_rows != B.n_rows || C.n_rows != D.n_rows) {
     Rcpp::stop("Matrices to be combined in the same row must have the same number of rows.");
   }
   if (A.n_cols != C.n_cols || B.n_cols != D.n_cols) {
     Rcpp::stop("Matrices to be combined in the same column must have the same number of columns.");
   }

   // Combine horizontally
   arma::mat top = arma::join_rows(A, B);
   arma::mat bottom = arma::join_rows(C, D);

   // Combine vertically
   arma::mat result = arma::join_cols(top, bottom);

   return result;
 }

 //' Get the optimal bandwidth for given time points
 //'
 //' This function selects the optimal bandwidth for each pair of time points (s, t)
 //' by finding the minimum risk value in the autocovariance risk matrix.
 //'
 //' @param mat_autocov_risk A numeric matrix containing the autocovariance risk estimates. The columns should include:
 //' \itemize{
 //'   \item \code{0}: s - The first argument of the autocovariance function.
 //'   \item \code{1}: t - The second argument of the autocovariance function.
 //'   \item \code{2}: hs - The candidate bandwidth for the first argument of the autocovariance function.
 //'   \item \code{3}: ht - The candidate bandwidth for the second argument of the autocovariance function.
 //'   \item \code{4}: PNl - The number of curves used in the estimation of the autocovariance at (s,t).
 //'   \item \code{5}: locreg_bw - The bandwidth used to estimate the local regularity parameters.
 //'   \item \code{6}: Hs - The estimates of the local exponent for each t, corresponding to H_s.
 //'   \item \code{7}: Ls - The estimates of the Hölder constant for each t, corresponding to L_s^2.
 //'   \item \code{8}: Ht - The estimates of the local exponent for each t, corresponding to H_t.
 //'   \item \code{9}: Lt - The estimates of the Hölder constant for each t, corresponding to L_t^2.
 //'   \item \code{10}: bias_term - The bias term of the risk function.
 //'   \item \code{11}: variance_term - The variance term of the risk function.
 //'   \item \code{12}: dependence_term - The dependence term of the risk function.
 //'   \item \code{13}: autocov_risk - The estimates of the risk function of the covariance/autocovariance function.
 //' }
 //' @param s A numeric vector of time points for which the bandwidth needs to be optimized.
 //' @param t A numeric vector of time points for which the bandwidth needs to be optimized.
 //' @return A numeric matrix with 8 columns:
 //' \itemize{
 //'   \item \code{s}: The input time point \code{s}.
 //'   \item \code{t}: The input time point \code{t}.
 //'   \item \code{Hs}: The local regularity estimate for \code{s}.
 //'   \item \code{Ls^2}: The local Hölder constant estimate for \code{s}.
 //'   \item \code{Ht}: The local regularity estimate for \code{t}.
 //'   \item \code{Lt^2}: The local Hölder constant estimate for \code{t}.
 //'   \item \code{optbw_s}: The optimal bandwidth for \code{s}.
 //'   \item \code{optbw_t}: The optimal bandwidth for \code{t}.
 //' }
 //'
 //' @examples
 //' \dontrun{
 //'   mat_autocov_risk <- matrix(runif(14 * 100), ncol = 14)
 //'   s <- seq(0, 1, length.out = 10)
 //'   t <- seq(0, 1, length.out = 10)
 //'   best_bw <- get_best_bw(mat_autocov_risk, s, t)
 //' }
 //' @export
 // [[Rcpp::export]]
 arma::mat get_best_bw(const arma::mat& mat_autocov_risk, const arma::vec s, const arma::vec t) {
   int n = s.size();
   arma::mat mat_res(n, 8);
   mat_res.col(0) = s;
   mat_res.col(1) = t;

   for (int k = 0; k < n; ++k) {
     arma::uvec idx_risk_cur = arma::find(mat_autocov_risk.col(0) == s(k) && mat_autocov_risk.col(1) == t(k));
     arma::vec risk = mat_autocov_risk(idx_risk_cur, arma::uvec({13}));
     arma::uword idx_min = arma::index_min(risk.elem(arma::find_finite(risk)));

     mat_res(k, 2) = mat_autocov_risk(idx_risk_cur(idx_min), 6); // H_s
     mat_res(k, 3) = mat_autocov_risk(idx_risk_cur(idx_min), 7); // L_s^2
     mat_res(k, 4) = mat_autocov_risk(idx_risk_cur(idx_min), 8); // H_t
     mat_res(k, 5) = mat_autocov_risk(idx_risk_cur(idx_min), 9); // L_t^2
     mat_res(k, 6) = mat_autocov_risk(idx_risk_cur(idx_min), 2); // optbw_s
     mat_res(k, 7) = mat_autocov_risk(idx_risk_cur(idx_min), 3); // optbw_t
   }

   return mat_res;
 }

 //' Get Nearest Best Bandwidth
 //'
 //' This function finds the nearest optimal bandwidth parameters for new time points using a nearest neighbor strategy.
 //'
 //' @param mat_opt_param A \code{matrix} containing the optimal parameters, typically the output of the \code{get_best_bw} function.
 //' The matrix should have columns representing \code{s}, \code{t}, \code{Hs}, \code{Ls^2}, \code{Ht}, \code{Lt^2}, \code{optbw_s}, and \code{optbw_t}.
 //' @param snew A numeric vector specifying the new \code{s} time points.
 //' @param tnew A numeric vector specifying the new \code{t} time points.
 //' @return A \code{matrix} with four columns: \code{snew}, \code{tnew}, \code{optbw_s}, and \code{optbw_t} for the nearest optimal bandwidths.
 //' @export
 // [[Rcpp::export]]
 arma::mat get_nearest_best_bw(const arma::mat& mat_opt_param, const arma::vec snew, const arma::vec tnew) {
   int n = snew.size();
   arma::mat mat_res(n, 4);

   mat_res.col(0) = snew;
   mat_res.col(1) = tnew;

   for (int k = 0; k < n; ++k) {
     // Calculate squared distances from (snew[k], tnew[k]) to all (s, t) in mat_opt_param
     arma::vec dist = arma::square(mat_opt_param.col(0) - snew(k)) + arma::square(mat_opt_param.col(1) - tnew(k));

     // Find the index of the minimum distance
     arma::uword idx_min_dist = arma::index_min(dist);

     // Extract optimal bandwidths corresponding to the minimum distance index
     mat_res(k, 2) = mat_opt_param(idx_min_dist, 6);
     mat_res(k, 3) = mat_opt_param(idx_min_dist, 7);
   }

   return mat_res;
 }




 // [[Rcpp::export]]
 Rcpp::List estimate_curve(const Rcpp::DataFrame data,
                           const arma::vec t,
                           const Rcpp::Nullable<int> id_curve = R_NilValue,
                           const Rcpp::Nullable<arma::vec> optbw_s = R_NilValue,
                           const Rcpp::Nullable<arma::vec> optbw_t = R_NilValue,
                           const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                           const bool use_same_bw = false,
                           const bool center = true,
                           const bool correct_diagonal = true,
                           const std::string kernel_name = "epanechnikov"){
   // Take unique observation points t
   int nt_pred = t.size();
   if ( nt_pred == 0) {
     stop("'t' must be a numeric vectors or scalar value(s) between 0 and 1.");
   }
   arma::vec tvec = arma::unique(t);
   int n = tvec.size();

   // Check if kernel_name is one of the supported kernels
   // and select the kernel function
   std::vector<std::string> supported_kernels = {"epanechnikov", "biweight", "triweight",
                                                 "tricube", "triangular", "uniform"};
   if (std::find(supported_kernels.begin(), supported_kernels.end(), kernel_name) == supported_kernels.end()) {
     stop("Unsupported kernel name. Choose from: epanechnikov, biweight, triweight, tricube, triangular, uniform.");
   }
   std::function<arma::vec(const arma::vec)> kernel_func = select_kernel(kernel_name);

   // Convert data to matrix
   arma::mat data_mat(data.nrows(), 3);
   data_mat.col(0) = as<arma::vec>(data["id_curve"]);
   data_mat.col(1) = as<arma::vec>(data["tobs"]);
   data_mat.col(2) = as<arma::vec>(data["X"]);
   arma::vec unique_id_curve = arma::unique(data_mat.col(0));
   double n_curve = unique_id_curve.n_elem;

   // Init the index of the curve to be predicted
   double idx_curve_pred = 0;
   double idx_curve_lag = 0;
   if (id_curve.isNull()){
     idx_curve_pred = n_curve;
     idx_curve_lag = idx_curve_pred - 1;
   } else {
     double idx_temp = Rcpp::as<int>(id_curve);
     if (idx_curve_pred > n_curve || idx_curve_pred > 1) {
       stop("If 'id_curve' is not NULL, then it must be an integer between 2 and the total number of curves.");
     } else {
       idx_curve_pred = idx_temp;
       idx_curve_lag = idx_curve_pred - 1;
     }
   }

   // Next, we extract idx_curve_pred and idx_curve_pred - 1
   // Extract the observation point, take unique values and sort them
   arma::uvec idx_cur_pred = arma::find(data_mat.col(0) == idx_curve_pred);
   arma::uvec idx_cur_lag = arma::find(data_mat.col(0) == idx_curve_lag);
   arma::vec Tvec_pred = data_mat(idx_cur_pred, arma::uvec({1}));
   arma::vec Tvec_lag = data_mat(idx_cur_lag, arma::uvec({1}));
   arma::vec Yvec_pred = data_mat(idx_cur_pred, arma::uvec({2}));
   arma::vec Yvec_lag = data_mat(idx_cur_lag, arma::uvec({2}));

   Rcout << "--> Set up : ok \n ";
   // Estimate the observation error standard deviation
   arma::mat mat_sig_pred = estimate_sigma_cpp(data, Tvec_pred);
   arma::mat mat_sig_lag = estimate_sigma_cpp(data, Tvec_lag);
   arma::mat Sigma_pred = arma::diagmat(arma::pow(mat_sig_pred.col(1), 2));
   arma::mat Sigma_lag = arma::diagmat(arma::pow(mat_sig_lag.col(1), 2));
   Rcout << "--> Sigma estimation : ok \n ";

   // Estimate bandwidth parameters on a grid
   // // Estimate cov and autocovariance on the grid
   arma::vec vec_grid = arma::vec({0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95});
   arma::mat grid_fixe = build_grid(vec_grid, vec_grid);
   arma::mat mat_cov_risk = estimate_autocov_risk_cpp(data, grid_fixe.col(0), grid_fixe.col(1), 0, R_NilValue, use_same_bw, center, kernel_name);
   arma::mat mat_autocov_risk = estimate_autocov_risk_cpp(data, grid_fixe.col(0), grid_fixe.col(1), 1, R_NilValue, use_same_bw, center, kernel_name);

   // // get optimal cov and autocov variance parameters
   mat mat_opt_cov_param = get_best_bw(mat_cov_risk, grid_fixe.col(0), grid_fixe.col(1));
   mat mat_opt_autocov_param = get_best_bw(mat_autocov_risk, grid_fixe.col(0), grid_fixe.col(1));


   // Estimate the covariances
   // // Get optimal bandwidth parameter
   arma::mat grid_pred_pred = build_grid(Tvec_pred, Tvec_pred);
   arma::mat grid_lag_lag = build_grid(Tvec_lag, Tvec_lag);
   arma::mat grid_pred_lag = build_grid(Tvec_pred, Tvec_lag);
   arma::mat grid_lag_pred = build_grid(Tvec_lag, Tvec_pred);
   arma::mat grid_lag_tvec = build_grid(Tvec_lag, tvec);
   arma::mat grid_pred_tvec = build_grid(Tvec_pred, tvec);

   arma::mat grid_pred_pred_optbw = get_nearest_best_bw(mat_opt_cov_param, grid_pred_pred.col(0), grid_pred_pred.col(1));
   arma::mat grid_lag_lag_optbw = get_nearest_best_bw(mat_opt_cov_param, grid_lag_lag.col(0), grid_lag_lag.col(1));
   arma::mat grid_pred_lag_optbw = get_nearest_best_bw(mat_opt_autocov_param, grid_pred_lag.col(0), grid_pred_lag.col(1));
   arma::mat grid_lag_pred_optbw = get_nearest_best_bw(mat_opt_autocov_param, grid_lag_pred.col(0), grid_lag_pred.col(1));
   arma::mat grid_lag_tvec_optbw = get_nearest_best_bw(mat_opt_autocov_param, grid_lag_tvec.col(0), grid_lag_tvec.col(1));
   arma::mat grid_pred_tvec_optbw = get_nearest_best_bw(mat_opt_cov_param, grid_pred_tvec.col(0), grid_pred_tvec.col(1));

   // // estimate covariances and autocovariances
   arma::mat mat_cov_pred_pred_all = estimate_autocov_cpp(data, grid_pred_pred_optbw.col(0), grid_pred_pred_optbw.col(1), 0,
                                                          Rcpp::wrap(grid_pred_pred_optbw.col(2)), Rcpp::wrap(grid_pred_pred_optbw.col(3)),
                                                          bw_grid, use_same_bw, center, correct_diagonal, kernel_name);
   arma::mat mat_cov_lag_lag_all = estimate_autocov_cpp(data, grid_lag_lag_optbw.col(0), grid_lag_lag_optbw.col(1), 0,
                                                        Rcpp::wrap(grid_lag_lag_optbw.col(2)), Rcpp::wrap(grid_lag_lag_optbw.col(3)),
                                                        bw_grid, use_same_bw, center, correct_diagonal, kernel_name);
   arma::mat mat_autocov_pred_lag_all = estimate_autocov_cpp(data, grid_pred_lag_optbw.col(0), grid_pred_lag_optbw.col(1), 0,
                                                             Rcpp::wrap(grid_pred_lag_optbw.col(2)), Rcpp::wrap(grid_pred_lag_optbw.col(3)),
                                                             bw_grid, use_same_bw, center, correct_diagonal, kernel_name);
   arma::mat mat_autocov_lag_pred_all = estimate_autocov_cpp(data, grid_lag_pred_optbw.col(0), grid_lag_pred_optbw.col(1), 0,
                                                             Rcpp::wrap(grid_lag_pred_optbw.col(2)), Rcpp::wrap(grid_lag_pred_optbw.col(3)),
                                                             bw_grid, use_same_bw, center, correct_diagonal, kernel_name);
   arma::mat mat_autocov_lag_tvec_all = estimate_autocov_cpp(data, grid_lag_tvec_optbw.col(0), grid_lag_tvec_optbw.col(1), 0,
                                                             Rcpp::wrap(grid_lag_tvec_optbw.col(2)), Rcpp::wrap(grid_lag_tvec_optbw.col(3)),
                                                             bw_grid, use_same_bw, center, correct_diagonal, kernel_name);
   arma::mat mat_cov_pred_tvec_all = estimate_autocov_cpp(data, grid_pred_tvec_optbw.col(0), grid_pred_tvec_optbw.col(1), 0,
                                                          Rcpp::wrap(grid_pred_tvec_optbw.col(2)), Rcpp::wrap(grid_pred_tvec_optbw.col(3)),
                                                          bw_grid, use_same_bw, center, correct_diagonal, kernel_name);
   Rcout << "--> autocov estimation : ok \n ";


   // // Build the matrix VarY_mat
   arma::mat mat_cov_lag_lag = reshape_matrix(mat_cov_lag_lag_all, 0, 1, 13);
   arma::mat G0_lag_lag = mat_cov_lag_lag + Sigma_lag;
   arma::mat mat_autocov_lag_pred = reshape_matrix(mat_autocov_lag_pred_all, 0, 1, 13);
   arma::mat mat_autocov_pred_lag = reshape_matrix(mat_autocov_pred_lag_all, 0, 1, 13);
   arma::mat mat_cov_pred_pred = reshape_matrix(mat_cov_pred_pred_all, 0, 1, 13);
   arma::mat G0_pred_pred = mat_cov_pred_pred + Sigma_pred;
   arma::mat mat_VarY = combine_matrices(G0_lag_lag, mat_autocov_lag_pred, mat_autocov_pred_lag, G0_pred_pred) ;
   Rcout << "--> mat_VarY build : ok \n ";
   // Build the matrix covY_Xn0
   arma::mat mat_autocov_lag_tvec = reshape_matrix(mat_autocov_lag_tvec_all, 0, 1, 13);
   arma::mat mat_cov_pred_tvec = reshape_matrix(mat_cov_pred_tvec_all, 0, 1, 13);
   arma::mat covY_Xn0 = arma::join_cols(mat_autocov_lag_tvec, mat_cov_pred_tvec);
   Rcout << "--> covY_Xn0 build : ok \n ";

   // Build the vector Y_{n_0, 1} - M_{n_0, 1}
   //// Estimate mean function
   arma::mat mat_mean_n0 = estimate_mean_cpp(data, Tvec_pred, R_NilValue, R_NilValue, kernel_name);
   arma::mat mat_mean_n0_lag = estimate_mean_cpp(data, Tvec_lag, R_NilValue, R_NilValue, kernel_name);
   arma::vec Yn_minus_mean_n0 = Yvec_pred - mat_mean_n0.col(5);
   arma::vec Yn_minus_mean_n0_lag = Yvec_lag - mat_mean_n0_lag.col(5);
   arma::vec vec_Yn0_lag = arma::join_vert(Yn_minus_mean_n0_lag, Yn_minus_mean_n0);

   // Build the BLUP
   //// Estimate the mean function
   arma::mat mat_mean_tvec = estimate_mean_cpp(data, tvec, R_NilValue, R_NilValue, kernel_name);
   Rcout << "--> BULP : mean : ok \n ";
   //// Estimate Bn0
   arma::mat Bn0 = arma::solve(mat_VarY, covY_Xn0);
   Rcout << "--> BULP : Bn0 : ok \n ";
   //// Set the BLUP
   arma::vec vec_blup = mat_mean_tvec.col(5) + arma::trans(Bn0) * vec_Yn0_lag;
   Rcout << "--> BULP : vec_blup : ok \n ";
   // Return
   int n_res = tvec.n_elem;
   arma::mat mat_res(n_res, 3);
   mat_res.col(0) = tvec;
   mat_res.col(1) = mat_mean_tvec.col(5);
   mat_res.col(2) = vec_blup;

   Rcpp::List result;

   result["opt_cov_param"] = mat_opt_cov_param;
   result["opt_autocov_param"] = mat_opt_autocov_param;
   result["cov_pred_pred"] = mat_cov_pred_pred;
   result["Sigma_pred"] = Sigma_pred;
   result["cov_pred_pred_all"] = mat_cov_pred_pred_all;
   result["cov_lag_lag"] = mat_cov_lag_lag;
   result["Sigma_lag"] = Sigma_lag;
   result["cov_lag_lag_all"] = mat_cov_lag_lag_all;
   result["autocov_lag_pred"] = mat_autocov_lag_pred;
   result["autocov_lag_pred_all"] = mat_autocov_lag_pred_all;
   result["autocov_pred_lag"] = mat_autocov_pred_lag;
   result["autocov_pred_lag_all"] = mat_autocov_pred_lag_all;
   result["autocov_lag_tvec"] = mat_autocov_lag_tvec;
   result["autocov_lag_tvec_all"] = mat_autocov_lag_tvec_all;
   result["cov_pred_tvec"] = mat_cov_pred_tvec;
   result["cov_pred_tvec_all"] = mat_cov_pred_tvec_all;
   result["mat_VarY"] = mat_VarY;
   result["covY_Xn0"] = covY_Xn0;
   result["muhat"] = mat_mean_tvec.col(5);
   result["vec_Yn0_lag"] = vec_Yn0_lag;
   result["res_blup"] = mat_res;

   return result;
 }
