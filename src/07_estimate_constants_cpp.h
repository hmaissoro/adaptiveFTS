#ifndef ESTIMATE_CONSTANT_H
#define ESTIMATE_CONSTANT_H

#include <RcppArmadillo.h>

arma::mat estimate_sigma_cpp(const Rcpp::DataFrame data, const arma::vec t);
arma::mat estimate_empirical_mom_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                     const arma::vec h, const double mom_order,
                                     const double center, const std::string kernel_name = "epanechnikov");
arma::mat estimate_empirical_autocov_cpp(const Rcpp::DataFrame data, const arma::vec t,
                                         const arma::vec h, const arma::vec lag,
                                         const std::string kernel_name = "epanechnikov");
arma::mat estimate_empirical_XsXt_autocov_cpp(const Rcpp::DataFrame data, const arma::vec s,
                                              const arma::vec t, const int cross_lag,
                                              const arma::vec lag, const arma::vec h,
                                              const std::string kernel_name, const bool center);
arma::mat estimate_numerator_dependence_term_DD_cpp(const Rcpp::DataFrame data,
                                                    const arma::vec t,
                                                    const int max_lag,
                                                    const arma::vec h,
                                                    const std::string kernel_name,
                                                    const bool center);

#endif
