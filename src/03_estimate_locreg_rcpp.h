#ifndef ESTIMATE_LOCREG_H
#define ESTIMATE_LOCREG_H

#include <RcppArmadillo.h>

arma::mat estimate_locreg_cpp(const Rcpp::DataFrame data, const arma::vec t,
                              const bool center,
                              const std::string kernel_name = "epanechnikov",
                              const Rcpp::Nullable<arma::vec> h = R_NilValue,
                              const Rcpp::Nullable<double> Delta = R_NilValue);

#endif
