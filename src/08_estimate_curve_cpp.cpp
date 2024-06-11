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
   arma::mat Sigma_n0 = arma::diagmat(arma::pow(mat_sig_pred.col(1), 2));
   arma::mat Sigma_n0_lag = arma::diagmat(arma::pow(mat_sig_lag.col(1), 2));
   Rcout << "--> Sigma estimation : ok \n ";
   // Estimate the covariances
   arma::mat grid_n0 = build_grid(Tvec_pred, Tvec_pred);
   arma::mat grid_n0_lag = build_grid(Tvec_lag, Tvec_lag);
   arma::mat mat_cov_n0_all = estimate_autocov_cpp(data, grid_n0.col(0), grid_n0.col(1), 0, optbw_s, optbw_t, bw_grid, use_same_bw, center, correct_diagonal, kernel_name);
   arma::mat mat_cov_n0_lag_all = estimate_autocov_cpp(data, grid_n0_lag.col(0), grid_n0_lag.col(1), 0, optbw_s, optbw_t, bw_grid, use_same_bw, center, correct_diagonal, kernel_name);
   Rcout << "--> cov estimation : ok \n ";
   // Estimate the autocovariances
   arma::mat grid_n0__n0_lag = build_grid(Tvec_pred, Tvec_lag);
   arma::mat grid_n0_lag__n0 = build_grid(Tvec_lag, Tvec_pred);
   arma::mat mat_autocov_n0__n0_lag_all = estimate_autocov_cpp(data, grid_n0__n0_lag.col(0), grid_n0__n0_lag.col(1), 1, optbw_s, optbw_t, bw_grid, use_same_bw, center, correct_diagonal, kernel_name);
   arma::mat mat_autocov_n0_lag__n0_all = estimate_autocov_cpp(data, grid_n0_lag__n0.col(0), grid_n0_lag__n0.col(1), 1, optbw_s, optbw_t, bw_grid, use_same_bw, center, correct_diagonal, kernel_name);
   Rcout << "--> autocov estimation : ok \n ";
   // Estimate covariance and autocovariance related to prediction point
   arma::mat grid_n0_lag__tvec = build_grid(Tvec_lag, tvec);
   arma::mat grid_n0__tvec = build_grid(Tvec_pred, tvec);
   Rcout << "--> tvec cov and autocov estimation -- grid : ok \n ";
   arma::mat mat_autocov_n0_lag_tvec_all = estimate_autocov_cpp(data, grid_n0_lag__tvec.col(0), grid_n0_lag__tvec.col(1), 1, optbw_s, optbw_t, bw_grid, use_same_bw, center, correct_diagonal, kernel_name);
   Rcout << "--> tvec cov and autocov estimation -- autocov : ok \n ";
   arma::mat mat_cov_n0_tvec_all = estimate_autocov_cpp(data, grid_n0__tvec.col(0), grid_n0__tvec.col(1), 0, optbw_s, optbw_t, bw_grid, use_same_bw, center, correct_diagonal, kernel_name);
   Rcout << "--> tvec cov and autocov estimation : ok \n ";
   // Build the matrix VarY_mat
   arma::mat mat_cov_n0_lag = reshape_matrix(mat_cov_n0_lag_all, 0, 1, 13);
   arma::mat mat_cov_n0 = reshape_matrix(mat_cov_n0_all, 0, 1, 13);
   arma::mat mat_autocov_n0__n0_lag = reshape_matrix(mat_autocov_n0__n0_lag_all, 0, 1, 13);
   arma::mat mat_autocov_n0_lag__n0 = reshape_matrix(mat_autocov_n0_lag__n0_all, 0, 1, 13);
   arma::mat G_n0_lag = mat_cov_n0_lag + Sigma_n0_lag;
   arma::mat G_n0 = mat_cov_n0 + Sigma_n0;
   arma::mat mat_VarY = combine_matrices(G_n0_lag, mat_autocov_n0_lag__n0, mat_autocov_n0__n0_lag, G_n0) ;
   Rcout << "--> mat_VarY build : ok \n ";
   // Build the matrix covY_Xn0
   arma::mat mat_autocov_n0_lag_tvec = reshape_matrix(mat_autocov_n0_lag_tvec_all, 0, 1, 13);
   arma::mat mat_cov_n0_tvec = reshape_matrix(mat_cov_n0_tvec_all, 0, 1, 13);
   arma::mat covY_Xn0 = arma::join_cols(mat_autocov_n0_lag_tvec, mat_cov_n0_tvec);
   Rcout << "--> covY_Xn0 build : ok \n ";
   // Build the vector Y_{n_0, 1} - M_{n_0, 1}
   //// Estimate mean function
   arma::mat mat_mean_n0 = estimate_mean_cpp(data, Tvec_pred, R_NilValue, R_NilValue, kernel_name);
   arma::mat mat_mean_n0_lag = estimate_mean_cpp(data, Tvec_lag, R_NilValue, R_NilValue, kernel_name);
   arma::vec Yn_minus_mean_n0 = Yvec_pred - mat_mean_n0.col(5);
   arma::vec Yn_minus_mean_n0_lag = Yvec_lag - mat_mean_n0_lag.col(5);
   Rcout << "Print Yvec_pred : " << arma::trans(Yvec_pred) << "\n";
   Rcout << "Print mean Yvec_pred : " << arma::trans(mat_mean_n0.col(5)) << "\n";
   Rcout << "Print Yvec_pred - mean : " << arma::trans(Yn_minus_mean_n0) << "\n";
   Rcout << "Print Yvec_lag : " << arma::trans(Yvec_lag) << "\n";
   Rcout << "Print mean Yvec_lag : " << arma::trans(mat_mean_n0_lag.col(5)) << "\n";
   Rcout << "Print Yvec_lag - mean : " << arma::trans(Yn_minus_mean_n0_lag) << "\n";

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

   result["Cov_n0"] = mat_cov_n0;
   result["Sigma_n0"] = Sigma_n0;
   result["cov_n0_all"] = mat_cov_n0_all;
   result["cov_n0_lag"] = mat_cov_n0_lag;
   result["Sigma_n0_lag"] = Sigma_n0_lag;
   result["cov_n0_lag_all"] = mat_cov_n0_lag_all;
   result["autocov_n0_lag__n0"] = mat_autocov_n0_lag__n0;
   result["autocov_n0_lag__n0_all"] = mat_autocov_n0_lag__n0_all;
   result["autocov_n0__n0_lag"] = mat_autocov_n0__n0_lag;
   result["autocov_n0__n0_lag_all"] = mat_autocov_n0__n0_lag_all;
   result["autocov_n0_lag_tvec"] = mat_autocov_n0_lag_tvec;
   result["autocov_n0_lag_tvec_all"] = mat_autocov_n0_lag_tvec_all;
   result["cov_n0_tvec"] = mat_cov_n0_tvec;
   result["cov_n0_tvec_all"] = mat_cov_n0_tvec_all;
   result["mat_VarY"] = mat_VarY;
   result["covY_Xn0"] = covY_Xn0;
   result["muhat"] = mat_mean_tvec.col(5);
   result["vec_Yn0_lag"] = vec_Yn0_lag;
   result["res_blup"] = mat_res;

   return result;
 }
