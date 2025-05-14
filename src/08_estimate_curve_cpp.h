#ifndef ESTIMATE_CURVE_H
#define ESTIMATE_CURVE_H

#include <RcppArmadillo.h>

Rcpp::List estimate_curve_cpp(const Rcpp::DataFrame data,
                              const arma::vec t,
                              const Rcpp::Nullable<int> id_curve = R_NilValue,
                              const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                              const bool use_same_bw = false,
                              const bool center = true,
                              const bool correct_diagonal = true,
                              const std::string kernel_name = "epanechnikov");
arma::mat get_nearest_best_autocov_bw(const arma::mat& mat_opt_param, const arma::vec snew, const arma::vec tnew);
arma::mat reshape_matrix(const arma::mat& A,
                         const arma::uword idx_col_s,
                         const arma::uword idx_col_t,
                         const arma::uword idx_col_value);

#endif
