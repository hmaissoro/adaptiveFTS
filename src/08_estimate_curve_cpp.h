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

#endif
