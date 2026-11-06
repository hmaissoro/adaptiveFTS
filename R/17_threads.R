#' Number of threads of the compiled estimators
#'
#' The C++ estimators are parallelised with OpenMP and, by default, use every
#' available processor. `set_threads()` changes that number for the current R
#' session; `get_threads()` returns it.
#'
#' @details
#' The results do not depend on the number of threads: every parallel loop
#' computes each output on a single thread, with the same summation order.
#'
#' **Monte Carlo studies.** When replications run in parallel R processes (for
#' instance with `parallel::parLapply()` or the `future` package), each process
#' starts as many threads as there are processors, so the machine is
#' oversubscribed. Call `set_threads(1)` in every worker (or set the environment
#' variable `ADAPTIVEFTS_NUM_THREADS=1` before the workers start) and run one
#' worker per processor: the replications then run side by side at full speed.
#'
#' The environment variable `ADAPTIVEFTS_NUM_THREADS`, when set to a positive
#' integer, is applied when the package is loaded. OpenMP's own limits (such as
#' `OMP_THREAD_LIMIT`) still apply.
#'
#' @param n `integer (positive)`. Number of threads. `NULL` (default) uses every
#'   available processor.
#'
#' @return `set_threads()` returns the previous number of threads, invisibly;
#'   `get_threads()` returns the current number.
#'
#' @examples
#' old <- set_threads(1)
#' get_threads()
#' set_threads(old)
#'
#' @export
set_threads <- function(n = NULL) {
  if (is.null(n)) {
    n <- 0L
  } else if (! (is.numeric(n) && length(n) == 1 && is.finite(n) && n >= 1 && n == floor(n))) {
    stop("'n' must be NULL or a single positive integer.")
  }
  return(invisible(set_threads_cpp(as.integer(n))))
}

#' @rdname set_threads
#' @export
get_threads <- function() {
  return(get_threads_cpp())
}
