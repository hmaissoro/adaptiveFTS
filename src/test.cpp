#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat get_best_mean_bw(const arma::mat& mat_mean_risk, const arma::vec t) {
  int n = t.size();
  arma::mat mat_res(n, 4);
  mat_res.col(0) = t;

  for (int k = 0; k < n; ++k) {
    arma::uvec idx_risk_cur = arma::find(mat_mean_risk.col(0) == t(k));
    arma::vec risk = mat_mean_risk(idx_risk_cur, arma::uvec({9}));
    arma::uword idx_min = arma::index_min(risk.elem(arma::find_finite(risk)));

    mat_res(k, 1) = mat_autocov_risk(idx_risk_cur(idx_min), 4); // H_t
    mat_res(k, 2) = mat_autocov_risk(idx_risk_cur(idx_min), 5); // L_t^2
    mat_res(k, 3) = mat_autocov_risk(idx_risk_cur(idx_min), 1); // optbw_t
  }
  return mat_res;
}








