#ifndef ESTIMATE_COV_SEGMENT_H
#define ESTIMATE_COV_SEGMENT_H

#include <RcppArmadillo.h>

arma::mat estimate_cov_segment_risk_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                        const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                        const bool center = true,
                                        const std::string kernel_name = "epanechnikov",
                                        const Rcpp::Nullable<arma::vec> presmooth_bw = R_NilValue,
                                        const Rcpp::Nullable<double> Delta = R_NilValue,
                                        const Rcpp::Nullable<arma::vec> presmooth_bw_grid = R_NilValue,
                                        const Rcpp::Nullable<int> presmooth_nsubset = R_NilValue);

arma::mat estimate_cov_segment_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                   const Rcpp::Nullable<arma::vec> bw = R_NilValue,
                                   const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                   const bool center = true,
                                   const std::string kernel_name = "epanechnikov",
                                   const Rcpp::Nullable<arma::vec> presmooth_bw = R_NilValue,
                                   const Rcpp::Nullable<double> Delta = R_NilValue,
                                   const Rcpp::Nullable<arma::vec> presmooth_bw_grid = R_NilValue,
                                   const Rcpp::Nullable<int> presmooth_nsubset = R_NilValue);

#endif
