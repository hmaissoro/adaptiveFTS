% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/06_estimate_covariance_segment.R
\name{estimate_cov_segment_risk}
\alias{estimate_cov_segment_risk}
\title{Estimate the Risk of the Covariance Segment Function}
\usage{
estimate_cov_segment_risk(
  data,
  idcol = "id_curve",
  tcol = "tobs",
  ycol = "X",
  t = c(1/4, 1/2, 3/4),
  bw_grid = NULL,
  center = TRUE,
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

\item{t}{A numeric vector. Observation points where the mean function of the underlying process is estimated.}

\item{bw_grid}{A numeric vector. A bandwidth grid from which the best smoothing parameter is selected for each \code{t}.
Default is \code{NULL}, in which case it is defined as an exponential grid of \eqn{N \times \lambda}.}

\item{center}{Logical. If \code{TRUE}, centers the data before estimation. Default is \code{TRUE}.}

\item{kernel_name}{Character string. Specifies the kernel function for estimation; default is \code{"epanechnikov"}.
Supported kernels include: \code{"epanechnikov"}, \code{"biweight"}, \code{"triweight"}, \code{"tricube"},
\code{"triangular"}, and \code{"uniform"}.}
}
\value{
A \link[data.table]{data.table} with columns:
\itemize{
\item \code{t}: The points at which the risk function is estimated.
\item \code{h}: The candidate bandwidth.
\item \code{PN}: The number of curves used to estimate the mean at \code{t}, i.e., \eqn{P_N(t;h)}.
\item \code{locreg_bw}: The bandwidth used to estimate the local regularity parameters.
\item \code{Ht}: The estimates of the local exponent \eqn{H_t}.
\item \code{Lt}: The estimates of the Hölder constant \eqn{L_t^2}.
\item \code{bias_term}: The bias term of the risk function.
\item \code{variance_term}: The variance term of the risk function.
\item \code{dependence_term}: The dependence term of the risk function.
\item \code{cov_segment_risk}: The estimated risk of the covariance segment function.
}
}
\description{
Estimates the risk \eqn{R_{\Gamma_0}(t; h)} associated with the covariance segment line estimation
proposed by \insertCite{maissoro2024pred;textual}{adaptiveFTS}.
}
\details{
The local regularity parameters are estimated within the function using \link{estimate_locreg_cpp}.

The dependence term includes contributions from both a term based on \eqn{\mathbb{D}(t; h_t)} derived from fourth-moment tensors,
and an empirical autocovariance term computed using \link{estimate_empirical_XsXt_autocov_cpp}.
}
\examples{
# Example coming soon

}
\references{
\insertAllCited{}
}
\seealso{
\link{estimate_mean}, \link{estimate_locreg}, \link{estimate_sigma},
\link{estimate_nw}, \link{estimate_empirical_autocov}
}
