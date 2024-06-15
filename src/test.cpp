#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Corrects the covariance matrix by ensuring the diagonal and above-diagonal values meet specific conditions.
 //'
 //' This function reshapes the covariance matrix based on unique values of s and t, fills matrices for differences
 //' and bandwidths, and performs corrections on diagonal and above-diagonal values. The final corrected matrix
 //' replaces the original values based on these conditions.
 //'
 //' @param mat_cov The input covariance matrix.
 //' @param idx_col_s The index of the column representing 's' values.
 //' @param idx_col_t The index of the column representing 't' values.
 //' @param idx_col_diff_st The index of the column representing differences between s and t.
 //' @param idx_col_bw The index of the column representing bandwidth values.
 //' @param idx_col_cov The index of the column representing covariance values.
 //' @return A corrected covariance matrix.
 //' @export
 // [[Rcpp::export]]
 arma::mat diagonal_correct_test(const arma::mat& mat_cov,
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
     for (int j = 0; j < ns - i; ++j) {
       arma::uvec idx = arma::find(mat_cov.col(idx_col_s) == s_unique(i) && mat_cov.col(idx_col_t) == t_unique(j));
       if (!idx.is_empty() && idx(0) < mat_cov.n_rows) {
         mat_diff(i, j) = mat_cov(idx(0), idx_col_diff_st);
         mat_bw(i, j) = mat_cov(idx(0), idx_col_bw);
         mat_cov_reshape(i, j) = mat_cov(idx(0), idx_col_cov);
       }
     }
   }

   for (int ids = 0; ids < ns; ++ids) {
     int idt_replace = 0;
     for (int idt = 0; idt < ns - ids; ++idt) {
       if (mat_diff(ids, idt) > mat_bw(ids, idt)) {
         idt_replace = idt;
       }
     }
     double replace_cov = mat_cov_reshape(ids, idt_replace);

     // Replace diagonal and above diagonal values
     int n_diag_step = ns - ids - idt_replace;
     for (int idx_step = 0; idx_step < n_diag_step; ++idx_step) {
       if ((ids + idx_step) < ns && (idt_replace + idx_step) < nt) {
         mat_cov_reshape(ids + idx_step, idt_replace + idx_step) = replace_cov;
       }
     }

     // Handle case when ids is 0
     if (ids == 0) {
       for (int idx_rep = idt_replace; idx_rep < ns - ids; ++idx_rep) {
         if (idx_rep < nt) {
           mat_cov_reshape(ids, idx_rep) = replace_cov;
         }
       }
     } else {
       for (int idx_row = idt_replace + 1; idx_row < ns - ids; ++idx_row) {
         if ((ids - 1) < ns && (idx_row - 1) < nt) {
           mat_cov_reshape(ids, idx_row) = mat_cov_reshape(ids - 1, idx_row - 1);
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
         if (ids_replace_backward < ns) {
           mat_cov_reshape(idx_rep, 0) = mat_cov_reshape(ids_replace_backward, 0);
         }
       }
     }
   }

   // Remove the data after the anti-diagonal
   arma::mat mat_cov_res = arma::zeros(ns, nt);
   for (int i = 0; i < ns; ++i) {
     for (int j = 0; j < ns - i; ++j) {
       if (j < nt) {
         mat_cov_res(i, j) = mat_cov_reshape(i, j);
       }
     }
   }

   // Initialize the result matrix
   arma::mat result = mat_cov;

   // Fill the result matrix
   for (arma::uword i = 0; i < ns; ++i) {
     for (arma::uword j = 0; j < ns - i; ++j) {
       arma::uvec idx_to_set = arma::find((mat_cov.col(0) == s_unique(i) && mat_cov.col(1) == t_unique(j)) ||
         (mat_cov.col(1) == s_unique(i) && mat_cov.col(0) == t_unique(j)));
       if (!idx_to_set.is_empty() && idx_to_set(0) < mat_cov.n_rows) {
         result(idx_to_set, arma::uvec({idx_col_cov})).fill(mat_cov_res(i, j));
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
















