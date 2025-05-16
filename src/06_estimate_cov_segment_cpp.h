#ifndef ESTIMATE_COV_SEGMENT_H
#define ESTIMATE_COV_SEGMENT_H

#include <RcppArmadillo.h>

arma::mat estimate_cov_segment_risk_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                        const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                        const bool center = true,
                                        const std::string kernel_name = "epanechnikov");

arma::mat estimate_cov_segment_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                   const Rcpp::Nullable<arma::vec> optbw = R_NilValue,
                                   const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                   const bool center = true,
                                   const std::string kernel_name = "epanechnikov");

#endif
