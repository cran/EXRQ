% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/EXRQ.R
\name{EVI.CFG.func}
\alias{EVI.CFG.func}
\title{Hill Estimator of the Extreme Value Index}
\usage{
EVI.CFG.func(x, tol = 1e-04, min.prop = 0.3, taus)
}
\arguments{
\item{x}{the estimated quantiles at intermediate quantile levels}

\item{tol}{the tolerance level used for checking quantile crossing}

\item{min.prop}{the minimum proportion of quantiles that are estimated higher than the adjacent lower quantiles}

\item{taus}{the corresponding quantile levels}
}
\value{
The estimated extreme value index is returned. If the proportion of cases with quantile crossing is too high, an NA is returned.
}
\description{
Hill Estimator of the Extreme Value Index
}
\details{
The function estimates the extreme value index using Hill estimator based on the estimated
intermediate quantiles.
}
\references{
Chernozhukov, C., Fernandez-Val, I., and Galichon, A. (2010).
Quantile and probability curves without crossing. Econometrica, 78, 1093-1125.
}

