#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List get_eigen_decomp_cpp(const arma::mat A){
  // Perform SVD
  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, A, "std");

  Rcpp::List result;
  result["vectors"] = eigvec;
  result["values"] = eigval;
  return result;
}
