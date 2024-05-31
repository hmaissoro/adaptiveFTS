// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// biweight_kernel
arma::vec biweight_kernel(const arma::vec u);
RcppExport SEXP _adaptiveFTS_biweight_kernel(SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(biweight_kernel(u));
    return rcpp_result_gen;
END_RCPP
}
// triweight_kernel
arma::vec triweight_kernel(const arma::vec u);
RcppExport SEXP _adaptiveFTS_triweight_kernel(SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(triweight_kernel(u));
    return rcpp_result_gen;
END_RCPP
}
// tricube_kernel
arma::vec tricube_kernel(const arma::vec u);
RcppExport SEXP _adaptiveFTS_tricube_kernel(SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(tricube_kernel(u));
    return rcpp_result_gen;
END_RCPP
}
// epanechnikov_kernel
arma::vec epanechnikov_kernel(const arma::vec u);
RcppExport SEXP _adaptiveFTS_epanechnikov_kernel(SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(epanechnikov_kernel(u));
    return rcpp_result_gen;
END_RCPP
}
// triangular_kernel
arma::vec triangular_kernel(const arma::vec u);
RcppExport SEXP _adaptiveFTS_triangular_kernel(SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(triangular_kernel(u));
    return rcpp_result_gen;
END_RCPP
}
// uniform_kernel
arma::vec uniform_kernel(const arma::vec u);
RcppExport SEXP _adaptiveFTS_uniform_kernel(SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(uniform_kernel(u));
    return rcpp_result_gen;
END_RCPP
}
// estimate_nw_cpp
arma::mat estimate_nw_cpp(const arma::vec y, const arma::vec t, const arma::vec tnew, const arma::vec h, const std::string kernel_name);
RcppExport SEXP _adaptiveFTS_estimate_nw_cpp(SEXP ySEXP, SEXP tSEXP, SEXP tnewSEXP, SEXP hSEXP, SEXP kernel_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type tnew(tnewSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_name(kernel_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_nw_cpp(y, t, tnew, h, kernel_name));
    return rcpp_result_gen;
END_RCPP
}
// estimate_nw_bw_cpp
double estimate_nw_bw_cpp(const arma::vec y, const arma::vec t, const arma::vec bw_grid, const std::string kernel_name);
RcppExport SEXP _adaptiveFTS_estimate_nw_bw_cpp(SEXP ySEXP, SEXP tSEXP, SEXP bw_gridSEXP, SEXP kernel_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type bw_grid(bw_gridSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_name(kernel_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_nw_bw_cpp(y, t, bw_grid, kernel_name));
    return rcpp_result_gen;
END_RCPP
}
// get_nw_optimal_bw_cpp
double get_nw_optimal_bw_cpp(const Rcpp::DataFrame data, const Nullable<arma::vec> bw_grid, const Nullable<int> nsubset, const std::string kernel_name);
RcppExport SEXP _adaptiveFTS_get_nw_optimal_bw_cpp(SEXP dataSEXP, SEXP bw_gridSEXP, SEXP nsubsetSEXP, SEXP kernel_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Nullable<arma::vec> >::type bw_grid(bw_gridSEXP);
    Rcpp::traits::input_parameter< const Nullable<int> >::type nsubset(nsubsetSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_name(kernel_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(get_nw_optimal_bw_cpp(data, bw_grid, nsubset, kernel_name));
    return rcpp_result_gen;
END_RCPP
}
// estimate_locreg_cpp
arma::mat estimate_locreg_cpp(const Rcpp::DataFrame data, const arma::vec t, const bool center, const std::string kernel_name, const Rcpp::Nullable<arma::vec> h, const Rcpp::Nullable<double> Delta);
RcppExport SEXP _adaptiveFTS_estimate_locreg_cpp(SEXP dataSEXP, SEXP tSEXP, SEXP centerSEXP, SEXP kernel_nameSEXP, SEXP hSEXP, SEXP DeltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< const bool >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type h(hSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<double> >::type Delta(DeltaSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_locreg_cpp(data, t, center, kernel_name, h, Delta));
    return rcpp_result_gen;
END_RCPP
}
// estimate_mean_risk_cpp
arma::mat estimate_mean_risk_cpp(const Rcpp::DataFrame data, const arma::vec t, const Rcpp::Nullable<arma::vec> bw_grid, const std::string kernel_name);
RcppExport SEXP _adaptiveFTS_estimate_mean_risk_cpp(SEXP dataSEXP, SEXP tSEXP, SEXP bw_gridSEXP, SEXP kernel_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type bw_grid(bw_gridSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_name(kernel_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_mean_risk_cpp(data, t, bw_grid, kernel_name));
    return rcpp_result_gen;
END_RCPP
}
// estimate_mean_cpp
arma::mat estimate_mean_cpp(const Rcpp::DataFrame data, const arma::vec t, const Rcpp::Nullable<arma::vec> optbw, const Rcpp::Nullable<arma::vec> bw_grid, const std::string kernel_name);
RcppExport SEXP _adaptiveFTS_estimate_mean_cpp(SEXP dataSEXP, SEXP tSEXP, SEXP optbwSEXP, SEXP bw_gridSEXP, SEXP kernel_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type optbw(optbwSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type bw_grid(bw_gridSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_name(kernel_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_mean_cpp(data, t, optbw, bw_grid, kernel_name));
    return rcpp_result_gen;
END_RCPP
}
// estimate_autocov_risk_cpp
arma::mat estimate_autocov_risk_cpp(const Rcpp::DataFrame data, const arma::vec s, const arma::vec t, const int lag, const Rcpp::Nullable<arma::vec> bw_grid, const bool use_same_bw, const bool center, const std::string kernel_name);
RcppExport SEXP _adaptiveFTS_estimate_autocov_risk_cpp(SEXP dataSEXP, SEXP sSEXP, SEXP tSEXP, SEXP lagSEXP, SEXP bw_gridSEXP, SEXP use_same_bwSEXP, SEXP centerSEXP, SEXP kernel_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type bw_grid(bw_gridSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_same_bw(use_same_bwSEXP);
    Rcpp::traits::input_parameter< const bool >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_name(kernel_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_autocov_risk_cpp(data, s, t, lag, bw_grid, use_same_bw, center, kernel_name));
    return rcpp_result_gen;
END_RCPP
}
// build_grid
arma::mat build_grid(const arma::vec& u, const arma::vec& v);
RcppExport SEXP _adaptiveFTS_build_grid(SEXP uSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(build_grid(u, v));
    return rcpp_result_gen;
END_RCPP
}
// estimate_autocov_cpp
arma::mat estimate_autocov_cpp(const Rcpp::DataFrame data, const arma::vec s, const arma::vec t, const int lag, const arma::vec param_grid, const Rcpp::Nullable<arma::vec> optbw_s, const Rcpp::Nullable<arma::vec> optbw_t, const Rcpp::Nullable<arma::vec> bw_grid, const bool use_same_bw, const bool center, const std::string kernel_name);
RcppExport SEXP _adaptiveFTS_estimate_autocov_cpp(SEXP dataSEXP, SEXP sSEXP, SEXP tSEXP, SEXP lagSEXP, SEXP param_gridSEXP, SEXP optbw_sSEXP, SEXP optbw_tSEXP, SEXP bw_gridSEXP, SEXP use_same_bwSEXP, SEXP centerSEXP, SEXP kernel_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type param_grid(param_gridSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type optbw_s(optbw_sSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type optbw_t(optbw_tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type bw_grid(bw_gridSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_same_bw(use_same_bwSEXP);
    Rcpp::traits::input_parameter< const bool >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_name(kernel_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_autocov_cpp(data, s, t, lag, param_grid, optbw_s, optbw_t, bw_grid, use_same_bw, center, kernel_name));
    return rcpp_result_gen;
END_RCPP
}
// estimate_sigma_cpp
arma::mat estimate_sigma_cpp(const Rcpp::DataFrame data, const arma::vec t);
RcppExport SEXP _adaptiveFTS_estimate_sigma_cpp(SEXP dataSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_sigma_cpp(data, t));
    return rcpp_result_gen;
END_RCPP
}
// estimate_empirical_mom_cpp
arma::mat estimate_empirical_mom_cpp(const Rcpp::DataFrame data, const arma::vec t, const arma::vec h, const double mom_order, const double center, const std::string kernel_name);
RcppExport SEXP _adaptiveFTS_estimate_empirical_mom_cpp(SEXP dataSEXP, SEXP tSEXP, SEXP hSEXP, SEXP mom_orderSEXP, SEXP centerSEXP, SEXP kernel_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< const double >::type mom_order(mom_orderSEXP);
    Rcpp::traits::input_parameter< const double >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_name(kernel_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_empirical_mom_cpp(data, t, h, mom_order, center, kernel_name));
    return rcpp_result_gen;
END_RCPP
}
// estimate_empirical_autocov_cpp
arma::mat estimate_empirical_autocov_cpp(const Rcpp::DataFrame data, const arma::vec t, const arma::vec h, const arma::vec lag, const std::string kernel_name);
RcppExport SEXP _adaptiveFTS_estimate_empirical_autocov_cpp(SEXP dataSEXP, SEXP tSEXP, SEXP hSEXP, SEXP lagSEXP, SEXP kernel_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_name(kernel_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_empirical_autocov_cpp(data, t, h, lag, kernel_name));
    return rcpp_result_gen;
END_RCPP
}
// estimate_empirical_XsXt_autocov_cpp
arma::mat estimate_empirical_XsXt_autocov_cpp(const Rcpp::DataFrame data, const arma::vec s, const arma::vec t, const int cross_lag, const arma::vec lag, const arma::vec h, const std::string kernel_name, const bool center);
RcppExport SEXP _adaptiveFTS_estimate_empirical_XsXt_autocov_cpp(SEXP dataSEXP, SEXP sSEXP, SEXP tSEXP, SEXP cross_lagSEXP, SEXP lagSEXP, SEXP hSEXP, SEXP kernel_nameSEXP, SEXP centerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int >::type cross_lag(cross_lagSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_name(kernel_nameSEXP);
    Rcpp::traits::input_parameter< const bool >::type center(centerSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_empirical_XsXt_autocov_cpp(data, s, t, cross_lag, lag, h, kernel_name, center));
    return rcpp_result_gen;
END_RCPP
}
// reshape_matrix
arma::mat reshape_matrix(const arma::mat& A, const arma::uword idx_col_s, const arma::uword idx_col_t, const arma::uword idx_col_value);
RcppExport SEXP _adaptiveFTS_reshape_matrix(SEXP ASEXP, SEXP idx_col_sSEXP, SEXP idx_col_tSEXP, SEXP idx_col_valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type idx_col_s(idx_col_sSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type idx_col_t(idx_col_tSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type idx_col_value(idx_col_valueSEXP);
    rcpp_result_gen = Rcpp::wrap(reshape_matrix(A, idx_col_s, idx_col_t, idx_col_value));
    return rcpp_result_gen;
END_RCPP
}
// combine_matrices
arma::mat combine_matrices(const arma::mat& A, const arma::mat& B, const arma::mat& C, const arma::mat& D);
RcppExport SEXP _adaptiveFTS_combine_matrices(SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(combine_matrices(A, B, C, D));
    return rcpp_result_gen;
END_RCPP
}
// estimate_curve
Rcpp::List estimate_curve(const Rcpp::DataFrame data, const arma::vec t, const Rcpp::Nullable<int> id_curve, const Rcpp::Nullable<arma::vec> param_grid, const Rcpp::Nullable<arma::vec> optbw_s, const Rcpp::Nullable<arma::vec> optbw_t, const Rcpp::Nullable<arma::vec> bw_grid, const bool use_same_bw, const bool center, const std::string kernel_name);
RcppExport SEXP _adaptiveFTS_estimate_curve(SEXP dataSEXP, SEXP tSEXP, SEXP id_curveSEXP, SEXP param_gridSEXP, SEXP optbw_sSEXP, SEXP optbw_tSEXP, SEXP bw_gridSEXP, SEXP use_same_bwSEXP, SEXP centerSEXP, SEXP kernel_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<int> >::type id_curve(id_curveSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type param_grid(param_gridSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type optbw_s(optbw_sSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type optbw_t(optbw_tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type bw_grid(bw_gridSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_same_bw(use_same_bwSEXP);
    Rcpp::traits::input_parameter< const bool >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_name(kernel_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_curve(data, t, id_curve, param_grid, optbw_s, optbw_t, bw_grid, use_same_bw, center, kernel_name));
    return rcpp_result_gen;
END_RCPP
}
// get_pn_vec
arma::mat get_pn_vec(const arma::vec svec, const arma::vec Tnvec, const arma::vec optbw_s_to_use);
RcppExport SEXP _adaptiveFTS_get_pn_vec(SEXP svecSEXP, SEXP TnvecSEXP, SEXP optbw_s_to_useSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type svec(svecSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type Tnvec(TnvecSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type optbw_s_to_use(optbw_s_to_useSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pn_vec(svec, Tnvec, optbw_s_to_use));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_adaptiveFTS_biweight_kernel", (DL_FUNC) &_adaptiveFTS_biweight_kernel, 1},
    {"_adaptiveFTS_triweight_kernel", (DL_FUNC) &_adaptiveFTS_triweight_kernel, 1},
    {"_adaptiveFTS_tricube_kernel", (DL_FUNC) &_adaptiveFTS_tricube_kernel, 1},
    {"_adaptiveFTS_epanechnikov_kernel", (DL_FUNC) &_adaptiveFTS_epanechnikov_kernel, 1},
    {"_adaptiveFTS_triangular_kernel", (DL_FUNC) &_adaptiveFTS_triangular_kernel, 1},
    {"_adaptiveFTS_uniform_kernel", (DL_FUNC) &_adaptiveFTS_uniform_kernel, 1},
    {"_adaptiveFTS_estimate_nw_cpp", (DL_FUNC) &_adaptiveFTS_estimate_nw_cpp, 5},
    {"_adaptiveFTS_estimate_nw_bw_cpp", (DL_FUNC) &_adaptiveFTS_estimate_nw_bw_cpp, 4},
    {"_adaptiveFTS_get_nw_optimal_bw_cpp", (DL_FUNC) &_adaptiveFTS_get_nw_optimal_bw_cpp, 4},
    {"_adaptiveFTS_estimate_locreg_cpp", (DL_FUNC) &_adaptiveFTS_estimate_locreg_cpp, 6},
    {"_adaptiveFTS_estimate_mean_risk_cpp", (DL_FUNC) &_adaptiveFTS_estimate_mean_risk_cpp, 4},
    {"_adaptiveFTS_estimate_mean_cpp", (DL_FUNC) &_adaptiveFTS_estimate_mean_cpp, 5},
    {"_adaptiveFTS_estimate_autocov_risk_cpp", (DL_FUNC) &_adaptiveFTS_estimate_autocov_risk_cpp, 8},
    {"_adaptiveFTS_build_grid", (DL_FUNC) &_adaptiveFTS_build_grid, 2},
    {"_adaptiveFTS_estimate_autocov_cpp", (DL_FUNC) &_adaptiveFTS_estimate_autocov_cpp, 11},
    {"_adaptiveFTS_estimate_sigma_cpp", (DL_FUNC) &_adaptiveFTS_estimate_sigma_cpp, 2},
    {"_adaptiveFTS_estimate_empirical_mom_cpp", (DL_FUNC) &_adaptiveFTS_estimate_empirical_mom_cpp, 6},
    {"_adaptiveFTS_estimate_empirical_autocov_cpp", (DL_FUNC) &_adaptiveFTS_estimate_empirical_autocov_cpp, 5},
    {"_adaptiveFTS_estimate_empirical_XsXt_autocov_cpp", (DL_FUNC) &_adaptiveFTS_estimate_empirical_XsXt_autocov_cpp, 8},
    {"_adaptiveFTS_reshape_matrix", (DL_FUNC) &_adaptiveFTS_reshape_matrix, 4},
    {"_adaptiveFTS_combine_matrices", (DL_FUNC) &_adaptiveFTS_combine_matrices, 4},
    {"_adaptiveFTS_estimate_curve", (DL_FUNC) &_adaptiveFTS_estimate_curve, 10},
    {"_adaptiveFTS_get_pn_vec", (DL_FUNC) &_adaptiveFTS_get_pn_vec, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_adaptiveFTS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}