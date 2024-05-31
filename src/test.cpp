#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat get_pn_vec(const arma::vec svec,
                     const arma::vec Tnvec,
                     const arma::vec optbw_s_to_use) {
  // Compute the vector p_n(t;h) and update P_N(t;h)
  arma::vec pn_s = arma::regspace(0, svec.size() - 1);
  pn_s.transform([svec, optbw_s_to_use, Tnvec](int j) { return  arma::find(abs((Tnvec - svec(j))) <= optbw_s_to_use(j)).is_empty() ? 0 : 1 ;});
  arma::vec pn_s_bis(svec.size(), fill::zeros);
  int pn_s_temp = 0;
  for (int g = 0; g < svec.size(); ++g) {
    arma::vec diff_vec = abs((Tnvec - svec(g)) / optbw_s_to_use(g));
    arma::uvec idx_vec = arma::find(diff_vec <= 1);
    if (idx_vec.is_empty()) {
      pn_s_temp = 0;
    } else {
      pn_s_temp = 1;
    }
    pn_s_bis(g) = pn_s_temp;
  }
  arma::mat mat_res(svec.size(), 2);
  mat_res.col(0) = pn_s;
  mat_res.col(1) = pn_s_bis;
  return mat_res;
}









