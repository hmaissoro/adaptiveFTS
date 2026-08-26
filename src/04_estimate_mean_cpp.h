#ifndef ESTIMATE_MEAN_H
#define ESTIMATE_MEAN_H

#include <RcppArmadillo.h>

arma::mat estimate_mean_risk_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                 const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                                 const std::string kernel_name = "epanechnikov",
                                 const Rcpp::Nullable<arma::vec> presmooth_bw = R_NilValue,
                                 const Rcpp::Nullable<double> Delta = R_NilValue,
                                 const Rcpp::Nullable<arma::vec> presmooth_bw_grid = R_NilValue,
                                 const Rcpp::Nullable<int> presmooth_nsubset = R_NilValue);
arma::mat estimate_mean_cpp(const Rcpp::DataFrame data, const arma::vec t,
                            const Rcpp::Nullable<arma::vec> bw = R_NilValue,
                            const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue,
                            const std::string kernel_name = "epanechnikov",
                            const Rcpp::Nullable<arma::vec> presmooth_bw = R_NilValue,
                            const Rcpp::Nullable<double> Delta = R_NilValue,
                            const Rcpp::Nullable<arma::vec> presmooth_bw_grid = R_NilValue,
                            const Rcpp::Nullable<int> presmooth_nsubset = R_NilValue);

#endif
