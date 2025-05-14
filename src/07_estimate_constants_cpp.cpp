#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "02_smoothing_rcpp.h"
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat estimate_sigma_cpp(const Rcpp::DataFrame data, const arma::vec t){
  int n = t.size();
  // Smooth curve
  arma::mat data_mat(data.nrows(), 3);
  data_mat.col(0) = as<arma::vec>(data["id_curve"]);
  data_mat.col(1) = as<arma::vec>(data["tobs"]);
  data_mat.col(2) = as<arma::vec>(data["X"]);
  arma::vec unique_id_curve = arma::unique(data_mat.col(0));
  double n_curve = unique_id_curve.n_elem;

  // Initialize a matrix for output
  arma::mat mat_res_sig(n, 2);
  for(int j = 0; j < n; ++j){
    // Initialize a vector to store Z_n = 1/2 * (Y_{n,i_t} - Y_{n,j_t})^2
    arma::vec zn_vec(n_curve);
    for(int i = 0 ; i < n_curve; ++i){
      // Extrat the current curve index data
      double idx_cur_curve = unique_id_curve(i);
      arma::uvec indices_cur = arma::find(data_mat.col(0) == idx_cur_curve);
      arma::mat mat_cur = data_mat.rows(indices_cur);

      // Compute absolute differences between elements of tobs and t(j)
      arma::vec abs_diff = arma::abs(mat_cur.col(1) - t(j));

      // Sort indices based on the absolute differences
      arma::uvec idx = arma::sort_index(abs_diff);
      arma::vec Yn = mat_cur.col(2);
      double Zn = 0.5 * (Yn(idx(0)) - Yn(idx(1))) * (Yn(idx(0)) - Yn(idx(1)));

      // Store the value
      zn_vec(i) = Zn;
    }
    double sig = sqrt(arma::mean(zn_vec));
    mat_res_sig.row(j) = {t(j), sig};
  }
  return mat_res_sig;
}

// [[Rcpp::export]]
arma::mat estimate_empirical_mom_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                     const arma::vec h, const double mom_order,
                                     const double center,
                                     const std::string kernel_name = "epanechnikov"){
  int n = t.size();

  // Smooth curve
  arma::mat data_mat(data.nrows(), 3);
  data_mat.col(0) = as<arma::vec>(data["id_curve"]);
  data_mat.col(1) = as<arma::vec>(data["tobs"]);
  data_mat.col(2) = as<arma::vec>(data["X"]);
  arma::vec unique_id_curve = arma::unique(data_mat.col(0));

  // Initialize a matrix for output
  double n_curve = unique_id_curve.n_elem;
  arma::mat mat_res_nw(n_curve * n, 3);
  for (int i = 0 ; i < n_curve; ++i) {
    // Extrat the current curve index data
    double idx_cur_curve = unique_id_curve(i);
    arma::uvec indices_cur = arma::find(data_mat.col(0) == idx_cur_curve);
    arma::mat mat_cur = data_mat.rows(indices_cur);

    // Smooth using Nadaraya-Watson estimator
    arma::vec h_to_use(1, fill::value(h(idx_cur_curve - 1)));
    arma::vec Xhat = estimate_nw_cpp(mat_cur.col(2), mat_cur.col(1), t, h_to_use, kernel_name);

    // Store the data
    arma::mat cur_mat(n, 3);
    cur_mat.col(0) = idx_cur_curve * arma::ones(n);
    cur_mat.col(1) = t;
    cur_mat.col(2) = Xhat;
    mat_res_nw(span(i * n, (i + 1) * n - 1), span(0, 2)) = cur_mat;
  }

  // Declare a matrix to store the computed moment for all lags
  arma::mat mat_mom_res(n, 3);
  for (int idx_t = 0; idx_t < n; ++idx_t) {
    // Extract all smooth data for each t
    arma::uvec cur_idx = arma::find(mat_res_nw.col(1) == t(idx_t));
    arma::mat mat_res_nw_by_t = mat_res_nw.rows(cur_idx);
    arma::vec x = mat_res_nw_by_t.col(2);

    // Estimate the mean and center if center = true
    if (center) {
      double muhat = arma::mean(x.elem(arma::find_finite(x)));
      x -= muhat;
    }
    arma::vec x_pow = arma::pow(x, mom_order);

    // Remove NaN values and compute mom_order th moment
    arma::uvec finiteIndices = arma::find_finite(x_pow);
    double Ex_mom_order = 0;
    if (finiteIndices.is_empty()) {
      // If autocov_vec contains only NaN values
      x_pow.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
      Ex_mom_order = arma::mean(x_pow);
    } else {
      Ex_mom_order = arma::mean(x_pow.elem(finiteIndices));
    }

    // Store the result for one t
    mat_mom_res.row(idx_t) = {t(idx_t), mom_order, Ex_mom_order};
  }
  // Replace non-finite value by zero
  mat_mom_res.replace(datum::nan, 0);
  mat_mom_res.replace(datum::inf, 0);

  return mat_mom_res;
}

// [[Rcpp::export]]
arma::mat estimate_empirical_autocov_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                         const arma::vec h, const arma::vec lag,
                                         const std::string kernel_name = "epanechnikov"){
  int n = t.size();
  int nlag = lag.size();
  // Smooth curve
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
    arma::vec Xhat = estimate_nw_cpp(mat_cur.col(2), mat_cur.col(1), t, h_to_use, kernel_name);

    // Store the data
    arma::mat cur_mat(n, 3);
    cur_mat.col(0) = idx_cur_curve * arma::ones(n);
    cur_mat.col(1) = t;
    cur_mat.col(2) = Xhat;
    mat_res_nw(span(i * n, (i + 1) * n - 1), span(0, 2)) = cur_mat;
  }

  // Declare a matrix to store the computed autocovariance for all lags
  arma::mat mat_autocov_res(nlag * n, 3);
  for(int l = 0; l < nlag; ++l){
    // Declare a matrix to store the computed autocovariance for one lag
    arma::mat mat_cur_lag(n, 3);
    for(int idx_t = 0; idx_t < n; ++idx_t){
      // Extract all smooth data for each t
      arma::uvec cur_idx = arma::find(mat_res_nw.col(1) == t(idx_t));
      arma::mat mat_res_nw_by_t = mat_res_nw.rows(cur_idx);
      arma::vec x = mat_res_nw_by_t.col(2);
      double N = cur_idx.size();

      // Estimate the mean
      // arma::uvec finiteIndices_x = arma::find_finite(x);
      arma::vec x_for_mean = x.elem(arma::find_finite(x));
      double muhat = arma::mean(x_for_mean);
      arma::vec xk = x.subvec(0, N - 1 - lag(l));
      arma::vec xk_plus_lag = x.subvec(lag(l), N - 1);
      arma::vec autocov_vec = (xk - muhat) % (xk_plus_lag - muhat);

      // Remove NaN values and compute auto-covariance
      arma::uvec finiteIndices = arma::find_finite(autocov_vec);
      double autocov = 0;
      if (finiteIndices.is_empty()) {
        // If autocov_vec contains only NaN values
        autocov_vec.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
        autocov = arma::mean(autocov_vec);
      } else {
        autocov = arma::mean(autocov_vec.elem(finiteIndices));
      }

      // Store the result for one lag
      mat_cur_lag.row(idx_t) = {t(idx_t), lag(l), autocov};
    }
    // Append the result
    mat_autocov_res(span(l * n, (l + 1) * n - 1), span(0, 2)) = mat_cur_lag;
  }

  // Replace non-finite value by zero
  mat_autocov_res.replace(datum::nan, 0);
  mat_autocov_res.replace(datum::inf, 0);

  return mat_autocov_res;
}

// Function to estimate empirical autocovariance for curves X(s)X(t)
// [[Rcpp::export]]
arma::mat estimate_empirical_XsXt_autocov_cpp(const Rcpp::DataFrame data,
                                              const arma::vec s,
                                              const arma::vec t,
                                              const int cross_lag,
                                              const arma::vec lag,
                                              const arma::vec h,
                                              const std::string kernel_name,
                                              const bool center) {
  int n = s.size();
  int nlag = lag.size();
  double clag = cross_lag;

  // Smooth curves
  // Implement curve smoothing here
  arma::mat data_mat(data.nrows(), 3);
  data_mat.col(0) = as<arma::vec>(data["id_curve"]);
  data_mat.col(1) = as<arma::vec>(data["tobs"]);
  data_mat.col(2) = as<arma::vec>(data["X"]);
  arma::vec unique_id_curve = arma::unique(data_mat.col(0));

  // Initialize a matrix for output
  double n_curve = unique_id_curve.n_elem;
  arma::mat mat_res_nw(n_curve * n, 5);
  for(int i = 0 ; i < n_curve; ++i){
    // Extrat the current curve index data
    double idx_cur_curve = unique_id_curve(i);
    arma::uvec indices_cur = arma::find(data_mat.col(0) == idx_cur_curve);
    arma::mat mat_cur = data_mat.rows(indices_cur);

    // Smooth using Nadaraya-Watson estimator
    arma::vec h_to_use(1, fill::value(h(idx_cur_curve - 1)));
    arma::vec Xhat_s = estimate_nw_cpp(mat_cur.col(2), mat_cur.col(1), s, h_to_use, kernel_name);
    arma::vec Xhat_t = estimate_nw_cpp(mat_cur.col(2), mat_cur.col(1), t, h_to_use, kernel_name);

    // Store the data
    arma::mat cur_mat(n, 5);
    cur_mat.col(0) = idx_cur_curve * arma::ones(n);
    cur_mat.col(1) = s;
    cur_mat.col(2) = Xhat_s;
    cur_mat.col(3) = t;
    cur_mat.col(4) = Xhat_t;
    mat_res_nw(span(i * n, (i + 1) * n - 1), span(0, 4)) = cur_mat;
  }

  // Implement the \gamma_{cross_lag}(s,t) estimation
  // Estimate autocovariance
  // Declare a matrix to store the computed autocovariance for all lags
  arma::mat mat_autocov_XsXt_res(nlag * n, 6);
  for(int j = 0; j < n; ++j){
    // Extract all smooth data for each t
    arma::uvec cur_idx = arma::find((mat_res_nw.col(1) == s(j)) % (mat_res_nw.col(3) == t(j)));
    arma::mat mat_res_nw_by_st = mat_res_nw.rows(cur_idx);
    arma::vec Xhat_s = mat_res_nw_by_st.col(2);
    arma::vec Xhat_t = mat_res_nw_by_st.col(4);
    double N = cur_idx.size();

    // Estimate the mean
    double muhat_s = arma::mean(Xhat_s.elem(arma::find_finite(Xhat_s)));
    double muhat_t = arma::mean(Xhat_t.elem(arma::find_finite(Xhat_t)));
    // Center data if necessary
    if (center) {
      Xhat_s -= muhat_s;
      Xhat_t -= muhat_t;
    }
    arma::vec Xhat_s_n = Xhat_s.subvec(0, N - 1 - cross_lag);
    arma::vec Xhat_t_n_plus_cross_lag = Xhat_t.subvec(cross_lag, N - 1);
    arma::vec XsXt_cross_lag = Xhat_s_n % Xhat_t_n_plus_cross_lag;

    // Remove NaN values and compute auto-covariance
    arma::uvec finite_XsXt = arma::find_finite(XsXt_cross_lag);
    double gamma_cross_lag = arma::mean(XsXt_cross_lag.elem(finite_XsXt));

    // Compute the X(s)X_l(t) autocovariance for lag = 0, 1, ...
    arma::mat mat_autocov_lag(nlag, 6);
    for(int l = 0; l < nlag; ++l){
      // For the argument s
      arma::vec Xhat_s_i = Xhat_s.subvec(0, N - 1 - cross_lag - lag(l));
      arma::vec Xhat_s_i_plus_lag = Xhat_s.subvec(lag(l), N - 1 - cross_lag);

      // For the argument t
      arma::vec Xhat_t_i_plus_cross_lag = Xhat_t.subvec(cross_lag, N - 1 - lag(l));
      arma::vec Xhat_t_i_plus_cross_lag_plus_lag = Xhat_t.subvec(cross_lag + lag(l), N - 1);

      arma::vec autocov_vec_lag = ((Xhat_s_i % Xhat_t_i_plus_cross_lag) - gamma_cross_lag) %
        ((Xhat_s_i_plus_lag % Xhat_t_i_plus_cross_lag_plus_lag) - gamma_cross_lag);

      // Remove NaN values and compute auto-covariance
      arma::uvec finite_autocov_lag = arma::find_finite(autocov_vec_lag);
      double XsXt_autocov = 0;
      if (finite_autocov_lag.is_empty()) {
        // If autocov_vec contains only NaN values
        autocov_vec_lag.replace(arma::datum::nan, 0).replace(arma::datum::inf, 0).replace(-arma::datum::inf, 0);
        XsXt_autocov = arma::mean(autocov_vec_lag);
      } else {
        XsXt_autocov = arma::mean(autocov_vec_lag.elem(finite_autocov_lag));
      }


      // Store the result for the X(s)X_l(t) autocovariance for lag = 0, 1, ...
      mat_autocov_lag.row(l) = {s(j), t(j), clag, lag(l), gamma_cross_lag, XsXt_autocov};
    }

    mat_autocov_XsXt_res(span(j * nlag, (j + 1) * nlag - 1), span(0, 5)) = mat_autocov_lag;
  }

  // Return results
  // Replace non-finite value by zero and sort
  mat_autocov_XsXt_res.replace(datum::nan, 0);
  mat_autocov_XsXt_res.replace(datum::inf, 0);
  mat_autocov_XsXt_res = mat_autocov_XsXt_res.rows(arma::sort_index(mat_autocov_XsXt_res.col(0)));

  return mat_autocov_XsXt_res ;
}

// :::::::::::::: Estimate \mathbb{D}(t, h_t) numerator ::::::::::::
// Estimate \mathbb{D}(t) for a single t
// [[Rcpp::export]]
double estimate_numerator_DD_single_t(const arma::vec& Znvec_raw, const int max_lag) {
  arma::vec Znvec = Znvec_raw;
  Znvec.replace(datum::nan, 0.0);  // Replace NaN with zero

  int N = Znvec.n_elem;
  double D_val_1 = 0.0;
  double D_val_2 = 0.0;

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
        D_val_1 += 24.0 * std::abs(sum_k / max_k);
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
      D_val_2 += 6.0 * std::abs(sum_k / max_k);
    }
  }
  double D_val = D_val_1 + D_val_2;
  return D_val;
}


// Function to estimate empirical autocovariance for curves X(s)X(t)
// [[Rcpp::export]]
arma::mat estimate_numerator_dependence_term_DD_cpp(const Rcpp::DataFrame data,
                                                    const arma::vec& t,
                                                    const arma::vec& bw_vec,
                                                    const arma::vec& h,
                                                    const int max_lag,
                                                    const std::string kernel_name,
                                                    const bool center) {
  int n = t.size();
  int n_bw = bw_vec.size();

  // Prepare data
  arma::mat data_mat(data.nrows(), 3);
  data_mat.col(0) = as<arma::vec>(data["id_curve"]);
  data_mat.col(1) = as<arma::vec>(data["tobs"]);
  data_mat.col(2) = as<arma::vec>(data["X"]);
  arma::vec unique_id_curve = arma::unique(data_mat.col(0));
  int n_curve = unique_id_curve.n_elem;

  // Smooth curves once per curve
  arma::mat mat_res_nw(n_curve * n, 4);
  for (int i = 0; i < n_curve; ++i) {
    double idx_cur_curve = unique_id_curve(i);
    arma::uvec indices_cur = arma::find(data_mat.col(0) == idx_cur_curve);
    arma::mat mat_cur = data_mat.rows(indices_cur);

    arma::vec h_to_use(1, arma::fill::value(h(idx_cur_curve - 1)));
    arma::vec Xhat_t = estimate_nw_cpp(mat_cur.col(2), mat_cur.col(1), t, h_to_use, kernel_name);

    arma::mat cur_mat(n, 4);
    cur_mat.col(0).fill(idx_cur_curve);
    cur_mat.col(1) = t;
    cur_mat.col(3) = Xhat_t;

    mat_res_nw.rows(i * n, (i + 1) * n - 1) = cur_mat;
  }

  // Prepare result matrix
  arma::mat mat_DD_res(n * n_bw, 4);

  for (int b = 0; b < n_bw; ++b) {
    double bw = bw_vec(b);

    for (int j = 0; j < n; ++j) {
      double t_j = t(j);
      arma::uvec cur_idx = arma::find(mat_res_nw.col(1) == t_j);
      arma::mat mat_res_nw_by_t = mat_res_nw.rows(cur_idx);

      // Compute pn_vec for current bandwidth
      arma::vec pn_vec(mat_res_nw_by_t.n_rows);

      for (int i = 0; i < mat_res_nw_by_t.n_rows; ++i) {
        double idx_cur_curve = mat_res_nw_by_t(i, 0);
        arma::uvec indices_cur = arma::find(data_mat.col(0) == idx_cur_curve);
        arma::mat mat_cur = data_mat.rows(indices_cur);

        arma::vec Tn_t_diff_over_bw = (mat_cur.col(1) - t_j) / bw;
        pn_vec(i) = arma::find(arma::abs(Tn_t_diff_over_bw) <= 1).is_empty() ? 0 : 1;
      }

      arma::vec Xhat_t = mat_res_nw_by_t.col(3);
      Xhat_t.replace(arma::datum::nan, 0.0);
      Xhat_t.replace(arma::datum::inf, 0.0);

      double muhat_t = arma::mean(Xhat_t.elem(arma::find_finite(Xhat_t)));
      if (center) {
        Xhat_t -= muhat_t;
      }

      arma::vec Zn_vec_pn = pn_vec % Xhat_t;
      double PN = arma::accu(pn_vec);
      double DD_by_t = estimate_numerator_DD_single_t(Zn_vec_pn, max_lag);

      int idx_res = b * n + j;
      mat_DD_res(idx_res, 0) = t_j;
      mat_DD_res(idx_res, 1) = bw;
      mat_DD_res(idx_res, 2) = DD_by_t;
      mat_DD_res(idx_res, 3) = PN;
    }
  }

  return mat_DD_res;
}

