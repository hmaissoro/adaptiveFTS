#ifndef SMOOTHING_RCPP_H
#define SMOOTHING_RCPP_H

#include <RcppArmadillo.h>

arma::vec biweight_kernel(const arma::vec u);
arma::vec triweight_kernel(const arma::vec u);
arma::vec tricube_kernel(const arma::vec u);
arma::vec epanechnikov_kernel(const arma::vec u);
arma::vec triangular_kernel(const arma::vec u);
arma::vec uniform_kernel(const arma::vec u);
std::function<arma::vec(const arma::vec)> select_kernel(const std::string kernel_name = "epanechnikov");
arma::mat estimate_nw_cpp(const arma::vec y, const arma::vec t, const arma::vec tnew, const arma::vec h, const std::string kernel_name = "epanechnikov");
double estimate_nw_bw_cpp(const arma::vec y, const arma::vec t, const arma::vec bw_grid, const std::string kernel_name = "epanechnikov");
double get_nw_optimal_bw_cpp(const Rcpp::DataFrame data, const Rcpp::Nullable<arma::vec> bw_grid = R_NilValue, const Rcpp::Nullable<int> nsubset = R_NilValue, const std::string kernel_name = "epanechnikov");

#endif
