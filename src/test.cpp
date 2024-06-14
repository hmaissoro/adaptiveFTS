#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat diagonal_correct(const arma::mat& mat_cov,
                           const arma::uword idx_col_s,
                           const arma::uword idx_col_t,
                           const arma::uword idx_col_diff_st,
                           const arma::uword idx_col_bw,
                           const arma::uword idx_col_cov) {
  // Initialize the result matrix
  arma::mat result = mat_cov;

  arma::vec s_unique = arma::sort(arma::unique(mat_cov.col(idx_col_s)), "descend");
  arma::vec t_unique = arma::sort(arma::unique(mat_cov.col(idx_col_s)), "ascend");
  int ns = s_unique.size();
  int nt = t_unique.size();

  for (int idx_s = 0; idx_s < ns; ++idx_s) {
    int idx_t_replace = 0;
    for (int idx_t = 0; idx_t < nt; ++idx_t) {
      arma::uvec idx_srow_tcol = arma::find(result.col(idx_col_s) == s_unique(idx_s) && result.col(idx_col_t) == t_unique(idx_t));
      if (! idx_srow_tcol.is_empty()) {
        double d = mat_cov(idx_srow_tcol(0), idx_col_diff_st);
        double bw = mat_cov(idx_srow_tcol(0), idx_col_bw);
        if (d > bw) {
          idx_t_replace += 1;
        }
      }
    }
    arma::uvec idx_cov_replace = arma::find(mat_cov.col(idx_col_s) == s_unique(idx_s) && mat_cov.col(idx_col_t) == t_unique(idx_t_replace));
    double cov_replace = mat_cov(idx_cov_replace(0), idx_col_cov);

    int n_replace = std::min(nt - idx_t_replace, ns - idx_s);
    for (int step = 0; step < n_replace; ++step){
      // Replace in the diagonal way
      arma::uvec idx_to_replace = arma::find(mat_cov.col(idx_col_s) == s_unique(idx_s + step) && mat_cov.col(idx_col_t) == t_unique(idx_t_replace + step));
      if (!idx_to_replace.is_empty()) {
        result(idx_to_replace(0), idx_col_cov) = cov_replace;
      }

      // For the case where idx_s = 0
      if (idx_s == 0) {
        arma::uvec idx_to_replace = arma::find(mat_cov.col(idx_col_s) == s_unique(idx_s) && mat_cov.col(idx_col_t) == t_unique(idx_t_replace + step));
        if (!idx_to_replace.is_empty()) {
          result(idx_to_replace(0), idx_col_cov) = cov_replace;
        }
      }
    }
  }
  // For the case where idx_s = ns - 1
  arma::uvec idx_tcol = arma::sort(arma::find(result.col(idx_col_t) == t_unique(0)));

  int idx_repalce = 0;
  for (int i = 0; i < idx_tcol.size(); ++i) {
    double d = mat_cov(idx_tcol(i), idx_col_diff_st);
    double bw = mat_cov(idx_tcol(i), idx_col_bw);
    if (d > bw) {
      idx_repalce = idx_tcol(i);
    }
  }
  double cov_replace = mat_cov(idx_repalce, idx_col_cov);
  for (int i = 0; i < idx_tcol.size(); ++i) {
    if (idx_tcol(i) > idx_repalce) {
      result(idx_tcol(i), idx_col_cov) = cov_replace;
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
















