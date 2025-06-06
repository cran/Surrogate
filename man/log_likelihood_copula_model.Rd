% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood_copula_models.R
\name{log_likelihood_copula_model}
\alias{log_likelihood_copula_model}
\title{Computes loglikelihood for a given copula model}
\usage{
log_likelihood_copula_model(
  theta,
  X,
  Y,
  d1,
  d2,
  copula_family,
  cdf_X,
  cdf_Y,
  pdf_X,
  pdf_Y,
  return_sum = TRUE
)
}
\arguments{
\item{theta}{Copula parameter}

\item{X}{Numeric vector corresponding to first outcome variable.}

\item{Y}{Numeric vector corresponding to second outcome variable.}

\item{d1}{An integer vector. Indicates whether first variable is observed or
right-censored,
\itemize{
\item \code{d1[i] = 1} if \code{u[i]} corresponds to non-censored value
\item \code{d1[i] = 0} if \code{u[i]} corresponds to right-censored value
\item \code{d1[i] = -1} if \code{u[i]} corresponds to left-censored value
}}

\item{d2}{An integer vector. Indicates whether first variable is observed or
right-censored,
\itemize{
\item \code{d2[i] = 1} if \code{v[i]} corresponds to non-censored value
\item \code{d2[i] = 0} if \code{v[i]} corresponds to right-censored value
\item \code{d2[i] = -1} if \code{v[i]} corresponds to left-censored value
}}

\item{copula_family}{Copula family, one of the following:
\itemize{
\item \code{"clayton"}
\item \code{"frank"}
\item \code{"gumbel"}
\item \code{"gaussian"}
}}

\item{cdf_X}{Distribution function for the first outcome variable.}

\item{cdf_Y}{Distribution function for the second outcome variable.}

\item{pdf_X}{Density function for the first outcome variable.}

\item{pdf_Y}{Density function for the second outcome variable.}

\item{return_sum}{Return the sum of the individual loglikelihoods? If \code{FALSE},
a vector with the individual loglikelihood contributions is returned.}
}
\value{
Loglikelihood of the bivariate copula model evaluated in the observed
data.
}
\description{
\code{log_likelihood_copula_model()} computes the loglikelihood for a given
bivariate copula model and data set while allowin for right-censoring of both
outcome variables.
}
