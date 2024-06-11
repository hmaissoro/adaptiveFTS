#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

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








