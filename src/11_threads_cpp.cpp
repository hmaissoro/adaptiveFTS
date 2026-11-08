#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// Number of OpenMP threads of the compiled estimators. The setting is that of
// this package's OpenMP runtime (other packages, e.g. data.table, keep theirs).

//' Set the number of threads of the compiled estimators (C++ core)
//'
//' @param n Number of threads; `0` or less uses every available processor.
//' @return The previous number of threads.
//' @keywords internal
// [[Rcpp::export]]
int set_threads_cpp(int n) {
#ifdef _OPENMP
  int old = omp_get_max_threads();
  if (n <= 0) n = omp_get_num_procs();
  omp_set_num_threads(n);
  return old;
#else
  (void) n;
  return 1;
#endif
}

//' Get the number of threads of the compiled estimators (C++ core)
//'
//' @return The number of threads used by the next parallel computation.
//' @keywords internal
// [[Rcpp::export]]
int get_threads_cpp() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}
