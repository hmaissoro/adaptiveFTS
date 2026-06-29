## Submission summary

This is the first submission of adaptiveFTS to CRAN. The package implements the
adaptive estimation and prediction procedures for weakly dependent functional
time series of Maissoro, Patilea and Vimond (2024, 2025), with the numerical
core written in C++ via Rcpp/RcppArmadillo.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Hassan Maissoro <maissorohassan@gmail.com>'
  New submission

The remaining note is the standard "New submission" note.

## Test environments

* Local: Windows 11, R 4.2.1
* GitHub Actions (r-lib/actions):
  * ubuntu-latest: R-devel, R-release, R-oldrel-1
  * macOS-latest: R-release, R-oldrel-1
  * windows-latest: R-release, R-oldrel-1

## Downstream dependencies

There are currently no reverse dependencies (new submission).

## Notes for the CRAN team

* The package compiles C++ (Rcpp/RcppArmadillo) and uses OpenMP where available;
  the OpenMP header include is guarded with `#ifdef _OPENMP`, and all results are
  deterministic regardless of the number of threads.
* Examples that are computationally heavier are wrapped in `\donttest{}`.
