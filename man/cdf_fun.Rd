% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood_copula_models.R
\name{cdf_fun}
\alias{cdf_fun}
\title{Function factory for distribution functions}
\usage{
cdf_fun(para, family)
}
\arguments{
\item{para}{Parameter vector.}

\item{family}{Distributional family, one of the following:
\itemize{
\item \code{"normal"}: normal distribution where \code{para[1]} is the mean and \code{para[2]} is
the standard deviation.
\item \code{"logistic"}: logistic distribution as parameterized in \code{stats::plogis()} where
\code{para[1]} and \code{para[2]} correspond to \code{location} and \code{scale}, respectively.
\item \code{"t"}: t distribution as parameterized in \code{stats::pt()} where \code{para[1]} and
\code{para[2]} correspond to \code{ncp} and \code{df}, respectively.
}}
}
\value{
A distribution function that has a single argument. This is the
vector of values in which the distribution function is evaluated.
}
\description{
Function factory for distribution functions
}
