% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/07_estimate_constants.R
\name{estimate_empirical_XsXt_autocov}
\alias{estimate_empirical_XsXt_autocov}
\title{Estimate Empirical \eqn{X_0(s)X_{\ell}(t)} Autocovariance Function for \eqn{\ell} = 0, 1, ...}
\usage{
estimate_empirical_XsXt_autocov(
  data,
  idcol = NULL,
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  cross_lag = 1,
  lag = c(0, 1, 2),
  h = NULL,
  center = FALSE,
  kernel_name = "epanechnikov"
)
}
\arguments{
\item{data}{A \code{data.table} (or \code{data.frame}), a \code{list} of \code{data.table} (or \code{data.frame}), or a \code{list} of \code{list}.
\itemize{
\item{If \code{data.table}:}{
It should contain the raw curve observations in at least three columns.
\itemize{
\item{\code{idcol} :}{ The name of the column containing the curve index in the sample.
Each curve index is repeated according to the number of observation points.}
\item{\code{tcol} :}{ The name of the column with observation points associated with each curve index.}
\item{\code{ycol} :}{ The name of the column with observed values at each observation point for each curve index.}
}
}
\item{If \code{list} of \code{data.table}:}{
In this case, each element in the \code{list} represents the observation data of a curve in the form of a \code{data.table} or \code{data.frame}.
Each \code{data.table} contains at least two columns.
\itemize{
\item{\code{tcol} :}{ The name of the column with observation points for the curve.}
\item{\code{ycol} :}{ The name of the column with observed values for the curve.}
}
}
\item{If \code{list} of \code{list}:}{
In this case, \code{data} is a list where each element is the observation data of a curve, given as a \code{list} of two vectors.
\itemize{
\item{\code{tcol} :}{ The vector containing observation points for the curve.}
\item{\code{ycol} :}{ The vector containing observed values for the curve.}
}
}
}}

\item{idcol}{\code{character}. If \code{data} is given as a \code{data.table} or \code{data.frame}, this is the name of the column that holds the curve index.
Each curve index is repeated according to the number of observation points. If \code{data} is a \code{list} of \code{data.table} (or \code{data.frame}) or a \code{list} of \code{list}, set \code{idcol = NULL}.}

\item{tcol}{\code{character}. The name of the column (or vector) containing the observation points for the curves.}

\item{ycol}{\code{character}. The name of the column with observed values for the curves.}

\item{s}{\code{vector (numeric)}. First argument in \eqn{X_0(s)X_{\ell}(t)}, corresponding to observation points \code{s} in the pair (\code{s}, \code{t}).
Must be of the same length as \code{t}.}

\item{t}{\code{vector (numeric)}. Second argument in \eqn{X_0(s)X_{\ell}(t)}, corresponding to observation points \code{t} in the pair (\code{s}, \code{t}).
Must be of the same length as \code{s}.}

\item{cross_lag}{\code{integer (positive integer)}. The lag \eqn{\ell} in \eqn{X_0(s)X_{\ell}(t)}.}

\item{lag}{\code{vector (integer)}. Lag for the autocovariance of \eqn{X_0(s)X_{\ell}(t)}.
If \code{lag = NULL}, only \eqn{\mathbb{E}X_0(s)X_{\ell}(t)} is returned.}

\item{h}{\code{numeric (positive vector or scalar)}. Smoothing bandwidth parameter.
Defaults to \code{NULL}, in which case \code{h} is estimated via Cross-Validation on a subset of curves.
If \code{h} is a scalar, all curves are smoothed with the same bandwidth; if a vector, it should match the number of curves in \code{data}.}

\item{center}{\code{logical}. If \code{TRUE}, the estimated autocovariance is centered: \eqn{\mathbb{E}(X_0(s) - \mu(s))(X_{\ell}(t) - \mu(t))}. Defaults to \code{FALSE}, providing \eqn{\mathbb{E}X_0(s)X_{\ell}(t)}.}

\item{kernel_name}{\code{string}. Kernel function for estimation; defaults to "epanechnikov". Supported kernels are: "epanechnikov", "biweight", "triweight", "tricube", "triangular", and "uniform".}
}
\value{
A \code{data.table} with columns:
\itemize{
\item{s :}{ First argument in \eqn{X_0(s)X_{\ell}(t)}.}
\item{t :}{ Second argument in \eqn{X_0(s)X_{\ell}(t)}.}
\item{cross_lag :}{ Lag \eqn{\ell} in \eqn{X_0(s)X_{\ell}(t)}.}
\item{lag :}{ Lags for autocovariance estimation of \eqn{X_0(s)X_{\ell}(t)}; contains \code{NA} if \code{lag = NULL}.}
\item{EXsXt_cross_lag :}{ Mean of \eqn{X_0(s)X_{\ell}(t)}.}
\item{XsXt_autocov :}{ Autocovariance estimates of \eqn{X_0(s)X_{\ell}(t)} for each \code{lag}; contains \code{NA} if \code{lag = NULL}.}
}
}
\description{
This function estimates the empirical \eqn{X_0(s)X_{\ell}(t)} autocovariance function for \eqn{\ell} = 0, 1, ...,
used in the empirical study of the papers \insertCite{maissoro2024adaptive;textual}{adaptiveFTS} and \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
}
\examples{
\dontrun{
# Load data
data("data_far")

# Example 1: Estimate autocovariance without centering
dt_empirical_cov <- estimate_empirical_XsXt_autocov(
  data = data_far,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  cross_lag = 1,
  lag = c(0, 1, 2),
  h = 0.1,
  center = FALSE,
  kernel_name = "epanechnikov"
)
dt_empirical_cov

# Example 2: Estimate autocovariance with centering
dt_empirical_cov_centered <- estimate_empirical_XsXt_autocov(
  data = data_far,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  s = c(1/5, 2/5, 4/5),
  t = c(1/4, 1/2, 3/4),
  cross_lag = 1,
  lag = c(0, 1, 2),
  h = 0.1,
  center = TRUE,
  kernel_name = "epanechnikov"
)
dt_empirical_cov_centered
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link[=get_nw_optimal_bw]{get_nw_optimal_bw()}}.
}
