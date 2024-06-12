#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat get_nearest_best_bw(const arma::mat& mat_mat_opt_param,
                      const arma::vec snew,
                      const arma::vec tnew) {

  int n = snew.size();
  arma::mat mat_res(n, 4);
  mat_res.col(0) = snew;
  mat_res.col(1) = tnew;

  // Matching using the nearest neighbour strategy
  for (int k = 0; k < n; ++k) {
    // Find rows in mat_risk where the first column equals t(k)
    arma::vec dist = arma::square(mat_mat_opt_param.col(0) - snew(k)) + arma::square(mat_mat_opt_param.col(1) - tnew(k));

    // Find the minimum index in the risk column
    arma::uword idx_min_dist = arma::index_min(dist);

    // Extract values corresponding to the minimum risk index
    mat_res(k, 2) = mat_mat_opt_param(idx_min_dist, 6);
    mat_res(k, 3) = mat_mat_opt_param(idx_min_dist, 7);
  }

  return mat_res;
}








