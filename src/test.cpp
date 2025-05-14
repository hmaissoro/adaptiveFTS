#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "02_smoothing_rcpp.h"
using namespace Rcpp;
using namespace arma;


// Estimate \mathbb{D}(t) for a single t
// [[Rcpp::export]]
double estimate_numerator_DD_single_t(const arma::vec& Znvec_raw, const int max_lag) {
  arma::vec Znvec = Znvec_raw;
  Znvec.replace(datum::nan, 0.0);  // Replace NaN with zero

  int N = Znvec.n_elem;
  double D_val = 0.0;

  // Quartic term
  for (int l1 = 1; l1 <= max_lag - 2; ++l1) {
    for (int l2 = l1 + 1; l2 <= max_lag - 1; ++l2) {
      for (int l3 = l2 + 1; l3 <= max_lag; ++l3) {
        int max_k = N - l3;
        if (max_k <= 0) continue;

        double sum_k = 0.0;
        for (int k = 0; k < max_k; ++k) {
          sum_k += Znvec(k) * Znvec(k + l1) * Znvec(k + l2) * Znvec(k + l3);
        }
        D_val += 24.0 * std::abs(sum_k / max_k);
      }
    }
  }

  // Cubic-square term
  for (int l1 = 1; l1 <= max_lag - 1; ++l1) {
    for (int l2 = l1 + 1; l2 <= max_lag; ++l2) {
      int max_k = N - l2;
      if (max_k <= 0) continue;

      double sum_k = 0.0;
      for (int k = 0; k < max_k; ++k) {
        double z1 = Znvec(k);
        double z2 = Znvec(k + l1);
        double z3 = Znvec(k + l2);
        sum_k += z1 * z1 * z2 * z3 + z1 * z2 * z2 * z3 + z1 * z2 * z3 * z3;
      }
      D_val += 6.0 * std::abs(sum_k / max_k);
    }
  }

  return D_val;
}


// Function to estimate empirical autocovariance for curves X(s)X(t)
// [[Rcpp::export]]
arma::mat estimate_numerator_dependence_term_DD_cpp(const Rcpp::DataFrame data,
                                                    const arma::vec t,
                                                    const int max_lag,
                                                    const arma::vec h,
                                                    const std::string kernel_name,
                                                    const bool center) {
  int n = t.size();

  // Smooth curves
  // Implement curve smoothing here
  arma::mat data_mat(data.nrows(), 3);
  data_mat.col(0) = as<arma::vec>(data["id_curve"]);
  data_mat.col(1) = as<arma::vec>(data["tobs"]);
  data_mat.col(2) = as<arma::vec>(data["X"]);
  arma::vec unique_id_curve = arma::unique(data_mat.col(0));

  // Initialize a matrix for output
  double n_curve = unique_id_curve.n_elem;
  arma::mat mat_res_nw(n_curve * n, 3);
  for(int i = 0 ; i < n_curve; ++i){
    // Extrat the current curve index data
    double idx_cur_curve = unique_id_curve(i);
    arma::uvec indices_cur = arma::find(data_mat.col(0) == idx_cur_curve);
    arma::mat mat_cur = data_mat.rows(indices_cur);

    // Smooth using Nadaraya-Watson estimator
    arma::vec h_to_use(1, fill::value(h(idx_cur_curve - 1)));
    arma::vec Xhat_t = estimate_nw_cpp(mat_cur.col(2), mat_cur.col(1), t, h_to_use, kernel_name);

    // Store the data
    arma::mat cur_mat(n, 3);
    cur_mat.col(0) = idx_cur_curve * arma::ones(n);
    cur_mat.col(1) = t;
    cur_mat.col(2) = Xhat_t;
    mat_res_nw(span(i * n, (i + 1) * n - 1), span(0, 2)) = cur_mat;
  }

  // Implement the \gamma_{cross_lag}(s,t) estimation
  // Estimate autocovariance
  // Declare a matrix to store the computed autocovariance for all lags
  arma::mat mat_DD_res(n, 2);
  for(int j = 0; j < n; ++j){
    // Extract all smooth data for each t
    arma::uvec cur_idx = arma::find(mat_res_nw.col(1) == t(j));
    arma::mat mat_res_nw_by_st = mat_res_nw.rows(cur_idx);
    arma::vec Xhat_t = mat_res_nw_by_st.col(2);
    Xhat_t.replace(datum::nan, 0.0);
    Xhat_t.replace(datum::inf, 0.0);

    // Estimate the mean
    double muhat_t = arma::mean(Xhat_t.elem(arma::find_finite(Xhat_t)));

    // Center data if necessary
    if (center) {
      Xhat_t -= muhat_t;
    }

    // Estimate dependence coefficient
    double DD_by_t = estimate_DD_single_t(Xhat_t, max_lag);

    // Store the result
    mat_DD_res(j, 0) = t(j);
    mat_DD_res(j, 1) = DD_by_t;

  }

  // Return results

  return mat_DD_res ;
}

