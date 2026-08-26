#ifndef ESTIMATE_AUTOCOV_H
#define ESTIMATE_AUTOCOV_H

#include <RcppArmadillo.h>

arma::mat estimate_autocov_risk_cpp(const Rcpp::DataFrame data, const arma::vec s, const arma::vec t,
                                    const int lag, const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                    const bool common_bw = false, const bool center = true,
                                    const std::string kernel_name = "epanechnikov",
                                    const Rcpp::Nullable<arma::vec> presmooth_bw = R_NilValue,
                                    const Rcpp::Nullable<double> Delta = R_NilValue,
                                    const Rcpp::Nullable<arma::vec> presmooth_bw_grid = R_NilValue,
                                    const Rcpp::Nullable<int> presmooth_nsubset = R_NilValue);

arma::mat sort_by_columns(const arma::mat& mat, arma::uword first_col_idx, arma::uword second_col_idx);

arma::mat estimate_autocov_cpp(const Rcpp::DataFrame data, const arma::vec s, const arma::vec t,
                               const int lag,
                               const Rcpp::Nullable<arma::vec> bw_s = R_NilValue,
                               const Rcpp::Nullable<arma::vec> bw_t = R_NilValue,
                               const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                               const bool common_bw = false, const bool center = true,
                               const bool correct_diagonal = true,
                               const std::string kernel_name = "epanechnikov",
                               const Rcpp::Nullable<arma::vec> presmooth_bw = R_NilValue,
                               const Rcpp::Nullable<double> Delta = R_NilValue,
                               const Rcpp::Nullable<arma::vec> presmooth_bw_grid = R_NilValue,
                               const Rcpp::Nullable<int> presmooth_nsubset = R_NilValue);

#endif
