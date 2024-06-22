#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


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

















