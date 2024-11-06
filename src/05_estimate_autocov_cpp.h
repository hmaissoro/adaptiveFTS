#ifndef ESTIMATE_AUTOCOV_H
#define ESTIMATE_AUTOCOV_H

#include <RcppArmadillo.h>

arma::mat estimate_autocov_risk_cpp(const Rcpp::DataFrame data, const arma::vec s, const arma::vec t,
                                    const int lag, const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                    const bool use_same_bw = false, const bool center = true,
                                    const std::string kernel_name = "epanechnikov");

arma::mat sort_by_columns(const arma::mat& mat, arma::uword first_col_idx, arma::uword second_col_idx);

arma::mat estimate_autocov_cpp(const Rcpp::DataFrame data, const arma::vec s, const arma::vec t,
                               const int lag,
                               const Rcpp::Nullable<arma::vec> optbw_s = R_NilValue,
                               const Rcpp::Nullable<arma::vec> optbw_t = R_NilValue,
                               const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                               const bool use_same_bw = false, const bool center = true,
                               const bool correct_diagonal = true,
                               const std::string kernel_name = "epanechnikov");

#endif
