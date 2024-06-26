#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat ensure_positive_definite(arma::mat A, double c = 1e-6) {

  // Step 1: Eigenvalue decomposition
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, A, "std");

  // Step 2: Replace negative eigenvalues with the minimum positive eigenvalue multiplied by a constant
  // Find the indexes of positive eigenvalues
  arma::uvec idx_positive_eigenval = arma::find(eigval > 0);
  if (idx_positive_eigenval.size() == eigval.size()) {
    return A;
  } else {
    // Find of minimum positive eigenvalue
    double min_positive = arma::min(eigval(idx_positive_eigenval));
    if (min_positive == arma::datum::inf) {
      throw std::runtime_error("All eigenvalues are non-positive.");
    }

    eigval.transform([&](double val){ return val <= 0 ? min_positive * c : val; });

    // Step 3: Reconstruct the matrix
    arma::mat A_recons = eigvec * arma::diagmat(eigval) * eigvec.t();

    // Ensure the matrix is symmetric
    A_recons = 0.5 * (A_recons + A_recons.t());

    return A_recons;
  }
}
















