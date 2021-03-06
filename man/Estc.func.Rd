% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/EXRQ.R
\name{Estc.func}
\alias{Estc.func}
\title{Estimation of the C vector}
\usage{
Estc.func(y, x, tau = 0.99, M = 2)
}
\arguments{
\item{y}{a vector of n untransformed responses}

\item{x}{a n x p matrix of n observations and p predictors}

\item{tau}{an upper quantile level close to one}

\item{M}{a constant larger than one that is used for estimating the c vector and thus K(x) function. The default is two}
}
\value{
A p-dimensional vector is returned.
}
\description{
This function estimates the C vector involved in the function test.EVI for testing the constancy of EVI
}

