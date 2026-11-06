#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "00_internal.h"
using namespace Rcpp;
using namespace arma;

// Presmoothing bandwidth (one per curve) and neighbourhood length Delta, resolved
// once: h by cross-validation when NULL, Delta from lambda when NULL.
LocregCfg resolve_locreg(const core::Curves& C, const Rcpp::Nullable<arma::vec> h,
                         const Rcpp::Nullable<double> Delta,
                         const Rcpp::Nullable<arma::vec> presmooth_bw_grid,
                         const Rcpp::Nullable<int> presmooth_nsubset, core::Kern kern) {
  int n_curve = C.N;
  LocregCfg cfg;
  cfg.h = arma::vec(n_curve);
  if (h.isNull()) {
    // Get the bandwidth grid and estimate best bandwidth
    double hbest = nw_optimal_bw_impl(C, presmooth_bw_grid, presmooth_nsubset, kern);
    cfg.h.fill(hbest);
  } else {
    arma::vec hvec = Rcpp::as<arma::vec>(h);
    int h_size = hvec.size();
    if (h_size != n_curve) {
      if (h_size != 1) {
        Rcpp::stop("'h' must be a numeric vector with length equal to the number of curves in 'data' or a scalar value.");
      } else {
        cfg.h.fill(hvec(0));
      }
    } else {
      cfg.h = hvec;
    }
  }

  // Control on Delta. NB : the range test must read the *supplied* value; it
  // previously tested a variable still at 0, so a user-supplied Delta went
  // unchecked.
  if (Delta.isNull()) {
    double lambdahat = C.lambdahat();
    cfg.Delta = std::min(exp(- std::pow(log(lambdahat), 1.0 / 3)), 0.2);
  } else {
    cfg.Delta = as<double>(Delta);
    if (!(cfg.Delta > 0) || cfg.Delta > 1) {
      stop("If 'Delta' is not NULL, it must be a scalar value between 0 and 1.");
    }
  }
  return cfg;
}

// Local regularity at each t. Each curve is smoothed once at t1, t2, t3 (the
// windowed Nadaraya-Watson estimator, bandwidth h(id - 1)); for each point, the
// finite triplets trimmed at their 2.5% / 97.5% quantiles give H_t and L_t^2.
// Points sharing the same (t1, t2, t3) are pooled, as the former look-up of the
// long table of smoothed values did.
arma::mat locreg_impl(const core::Curves& C, const arma::vec& t, const bool center,
                      core::Kern kern, const LocregCfg& cfg) {
  if (t.size() == 0) {
    stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.");
  }
  const int n_curve = C.N;
  const double Delta_to_use = cfg.Delta;
  const arma::vec& hvec_to_use = cfg.h;

  // Take into account the cases where t-Delta/2 < 0 and t + Delta/2
  int n = t.size();
  arma::vec t1(n), t2(n), t3(n);

  // Deterministic tie-breaking offset. Previously this was a random draw in
  // [1e-6, 2e-6] via arma::randu(), which (a) is not controlled by R's
  // set.seed(), and (b) advances Armadillo's RNG on every call, so the local
  // regularity estimates were not reproducible (they varied by ~1e-4 between
  // otherwise identical calls). A fixed offset breaks ties just as well while
  // making the estimates fully reproducible across runs and platforms.
  const double jitter = 1.5e-06;
  for (int i = 0; i < n; ++i) {
    double ti = t(i);

    // Add a small offset to avoid ties
    if ((ti - Delta_to_use / 2) <= 0) {
      t1(i) = ti;
      t2(i) = ti + Delta_to_use / 2 + jitter;
      t3(i) = ti + Delta_to_use + jitter;
    } else if ((ti + Delta_to_use / 2) >= 1) {
      t3(i) = ti;
      t2(i) = ti - Delta_to_use / 2 + jitter;
      t1(i) = ti - Delta_to_use + jitter;
    } else {
      t1(i) = ti - Delta_to_use / 2 + jitter;
      t2(i) = ti;
      t3(i) = ti + Delta_to_use / 2 + jitter;
    }
  }

  // Smoothed curves at t1, t2, t3: n_curve x n matrices
  arma::mat X1(n_curve, n), X2(n_curve, n), X3(n_curve, n);
#pragma omp parallel for schedule(static)
  for (int c = 0; c < n_curve; ++c) {
    double hc = hvec_to_use(C.ids(c) - 1);
    for (int j = 0; j < n; ++j) {
      X1(c, j) = core::nw_at(C, c, t1(j), hc, kern);
      X2(c, j) = core::nw_at(C, c, t2(j), hc, kern);
      X3(c, j) = core::nw_at(C, c, t3(j), hc, kern);
    }
  }

  std::map<std::vector<double>, std::vector<arma::uword>> groups;
  for (int j = 0; j < n; ++j) groups[{t1(j), t2(j), t3(j)}].push_back(j);

  // Estimate locreg
  arma::mat mat_locreg(n, 6);
  const double locreg_bw = arma::median(hvec_to_use);
  for (int j = 0; j < n; ++j) {
    // All smoothed values for this t (curve by curve)
    const std::vector<arma::uword>& g = groups[{t1(j), t2(j), t3(j)}];
    const arma::uword len = n_curve * g.size();
    arma::vec xt1(len), xt2(len), xt3(len);
    arma::uword r = 0;
    for (int c = 0; c < n_curve; ++c)
      for (arma::uword jj : g) { xt1(r) = X1(c, jj); xt2(r) = X2(c, jj); xt3(r) = X3(c, jj); ++r; }

    // Remove NaN values
    arma::uvec finiteIndices = arma::find_finite(xt1 + xt2 + xt3);
    arma::vec filtered_xt1 = xt1.elem(finiteIndices);
    arma::vec filtered_xt2 = xt2.elem(finiteIndices);
    arma::vec filtered_xt3 = xt3.elem(finiteIndices);

    // Remove extreme values (outside the 2.5% / 97.5% quantiles; the former code
    // computed the same quantile once per element)
    const arma::vec probs = {0.025, 0.975};
    arma::uvec rxt1, rxt2, rxt3;
    if (finiteIndices.is_empty()) {
      rxt1 = rxt2 = rxt3 = arma::uvec();
    } else {
      arma::vec q1 = arma::quantile(filtered_xt1, probs);
      arma::vec q2 = arma::quantile(filtered_xt2, probs);
      arma::vec q3 = arma::quantile(filtered_xt3, probs);
      rxt1 = (filtered_xt1 >= q1(0)) % (filtered_xt1 <= q1(1));
      rxt2 = (filtered_xt2 >= q2(0)) % (filtered_xt2 <= q2(1));
      rxt3 = (filtered_xt3 >= q3(0)) % (filtered_xt3 <= q3(1));
    }

    arma::vec xt1_final = filtered_xt1.elem(arma::find(rxt1 % rxt2 % rxt3));
    arma::vec xt2_final = filtered_xt2.elem(arma::find(rxt1 % rxt2 % rxt3));
    arma::vec xt3_final = filtered_xt3.elem(arma::find(rxt1 % rxt2 % rxt3));

    // Center data if necessary
    if (center) {
      xt1_final -= arma::mean(xt1_final);
      xt2_final -= arma::mean(xt2_final);
      xt3_final -= arma::mean(xt3_final);
    }

    // Compute local regularity parameters
    double theta_t1_t3 = arma::mean(arma::square(xt1_final - xt3_final));
    double theta_t2_t3 = arma::mean(arma::square(xt2_final - xt3_final));
    double Ht_row = (std::log(theta_t1_t3) - std::log(theta_t2_t3)) / (2 * std::log(2));
    double Ht = (Ht_row <= 0.1) ? 0.1 : (Ht_row >=1) ? 1 : Ht_row;
    double Lt = theta_t2_t3 / std::pow(std::abs(t2(j) - t3(j)), 2 * Ht);

    // Store the result in matrix mat_locreg initialized above
    double Nused = xt1_final.size();
    mat_locreg.row(j) = {t(j), locreg_bw, Delta_to_use, Nused, Ht, Lt};
  }

  return mat_locreg;
}

//' Local Regularity Parameters Estimation
 //'
 //' @param data A DataFrame containing the columns "id_curve", "tobs", and "X". Typically, the output of the function \link{format_data}.
 //' @param t Numeric vector. Observation points at which we want to estimate the local regularity parameters of the underlying process.
 //' @param center Logical. If \code{TRUE}, the curves are centered.
 //' @param kernel_name String. The kernel function of the Nadaraya-Watson estimator. Default is "epanechnikov".
 //' @param h Numeric (positive vector or scalar). The bandwidth of the Nadaraya-Watson estimator for the local regularity estimation.
 //' Default \code{h = NULL} and thus it will be estimated by Cross-Validation on a subset of curves.
 //' If \code{h} is a scalar, then all curves will be smoothed with the same bandwidth.
 //' Otherwise, if \code{h} is a vector, its length must be equal to the number of curves in \code{data}
 //' and each element of the vector must correspond to a curve given in the same order as in \code{data}.
 //' @param Delta Numeric (positive). The length of the neighborhood of \code{t} around which the local regularity is to be estimated.
 //' Default \code{Delta = NULL} and thus it will be estimated from the data.
 //' @param presmooth_bw_grid Numeric vector. Candidate bandwidths of the cross-validation
 //' that selects \code{h} when \code{h = NULL}. Default \code{presmooth_bw_grid = NULL} uses
 //' the default grid of `get_nw_optimal_bw_cpp`. Ignored when \code{h} is supplied.
 //' @param presmooth_nsubset Integer (positive). Number of curves used by that cross-validation.
 //' Default \code{presmooth_nsubset = NULL} uses min(70, floor(N / 2)) curves.
 //' Ignored when \code{h} is supplied.
 //'
 //' @return A \code{matrix} containing the following six columns in order:
 //' \enumerate{
 //'   \item t: The points around which the local regularity parameters are estimated.
 //'   \item locreg_bw: The presmoothing bandwidth.
 //'   \item Delta: The length of the neighborhood of \code{t} around which the local regularity is to be estimated.
 //'   \item Nused: The number of curves that give non-degenerate estimates around \code{t}.
 //'   \item Ht: The local exponent estimates for each \code{t}.
 //'   \item Lt2: The Hölder constant estimates for each \code{t}.
 //' }
 //'
 //' @examples
 //' \dontrun{
 //' data <- data.frame(id_curve = rep(1:5, each = 10), tobs = rep(seq(0, 1, length.out = 10), 5), X = rnorm(50))
 //' estimate_locreg_cpp(data, seq(0, 1, length.out = 20), center = TRUE)
 //' }
 //'
 // [[Rcpp::export]]
 arma::mat estimate_locreg_cpp(const Rcpp::DataFrame data, const arma::vec t,
                               const bool center,
                               const std::string kernel_name = "epanechnikov",
                               const Rcpp::Nullable<arma::vec> h = R_NilValue,
                               const Rcpp::Nullable<double> Delta = R_NilValue,
                               const Rcpp::Nullable<arma::vec> presmooth_bw_grid = R_NilValue,
                               const Rcpp::Nullable<int> presmooth_nsubset = R_NilValue){
   if (t.size() == 0) {
     stop("'t' must be a numeric vector or scalar value(s) between 0 and 1.");
   }
   core::Kern kern = core::parse_kernel(kernel_name);
   core::Curves C = core::make_curves(data);
   LocregCfg cfg = resolve_locreg(C, h, Delta, presmooth_bw_grid, presmooth_nsubset, kern);
   return locreg_impl(C, t, center, kern, cfg);
 }
