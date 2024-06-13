#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

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
 // [[Rcpp::export]]
 arma::mat diagonal_correct(const arma::mat& mat,
                            const arma::uword idx_col_s,
                            const arma::uword idx_col_t,
                            const arma::uword idx_col_diff_st,
                            const arma::uword idx_col_bw,
                            const arma::uword idx_col_cov) {
   // Initialize the result matrix
   arma::mat result = mat;

   // Get unique values of s and t
   arma::vec s_unique = arma::sort(arma::unique(mat.col(idx_col_s)), "descend");
   arma::vec t_unique = arma::sort(arma::unique(mat.col(idx_col_t)), "ascend");

   // Loop through unique values of s
   for (arma::uword idx_s = 0; idx_s < s_unique.size(); ++idx_s) {
     // Find the rows where s equals s_unique(idx_s)
     arma::uvec idx_s_rows = arma::find(mat.col(idx_col_s) == s_unique(idx_s));

     // Find the rows where s and t match
     for (arma::uword idx_t = 0; idx_t < t_unique.size(); ++idx_t) {
       arma::uword idx_t_replace = idx_t;
       for (arma::uword idx_row : idx_s_rows) {
         if (mat(idx_row, idx_col_diff_st) > mat(idx_row, idx_col_bw) && mat(idx_row, idx_col_t) == t_unique(idx_t)) {
           idx_t_replace = idx_row;
           break;
         }
       }
       // Replace covariance values for the selected rows
       arma::uvec idx_cov_replace = arma::find(mat.col(idx_col_s) == s_unique(idx_s) && mat.col(idx_col_t) == t_unique(idx_t_replace));
       if (!idx_cov_replace.empty()) {
         double cov_replace = mat(idx_cov_replace(0), idx_col_cov);
         arma::uvec idx_to_replace = arma::find(mat.col(idx_col_s) == s_unique(idx_s) && mat.col(idx_col_t) >= t_unique(idx_t_replace));
         result(idx_to_replace(0), idx_col_cov) = cov_replace;
       }
     }
   }

   return result;
 }


 arma::mat diagonal_correct(const arma::mat& mat,
                            const arma::uword idx_col_s,
                            const arma::uword idx_col_t,
                            const arma::uword idx_col_diff_st,
                            const arma::uword idx_col_bw,
                            const arma::uword idx_col_cov) {
   // Initialize the result matrix
   arma::mat result = mat;

   arma::vec s_unique = arma::sort(arma::unique(mat.col(idx_col_s)), "descend");
   arma::vec t_unique = arma::sort(arma::unique(mat.col(idx_col_s)), "ascend");
   arma::uword ns = s_unique.size();
   arma::uword nt = t_unique.size();

   for (arma::uword idx_s = 0; idx_s < ns; ++idx_s) {
     arma::uword idx_t_replace = 0;
     for (arma::uword idx_t = 0; idx_t < nt; ++idx_t) {
       arma::uword idx_srow_tcol = arma::find(result.col(idx_col_s) == s_unique(idx_s) && result.col(idx_col_t) == t_unique(idx_t));
       if (mat(idx_srow_tcol, idx_col_diff_st) > result(idx_srow_tcol, idx_col_bw)) {
         idx_t_replace = idx_srow_tcol;
       }
     }
     arma::uword idx_cov_replace = arma::find(result.col(idx_col_s) == s_unique(idx_s) && result.col(idx_col_t) == t_unique(idx_t_replace));
     double cov_replace = result(idx_cov_replace, idx_col_cov);
     arma::uword n_replace = nt - idx_t_replace;
     for (arma::uword step = 0; step < n_replace; ++step){
       arma::uword idx_to_replace = arma::find(result.col(idx_col_s) == s_unique(idx_s + step) && result.col(idx_col_t) == t_unique(idx_t_replace + step));
       if (! (arma::Datum<arma::uword>::nan == idx_to_replace)) {
         result(idx_to_replace, idx_col_cov) = cov_replace;
       }
     }
   }

   return result;
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
















