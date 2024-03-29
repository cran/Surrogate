% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_model_BinCont_copula.R
\name{twostep_BinCont}
\alias{twostep_BinCont}
\title{Fit binary-continuous copula submodel with two-step estimator}
\usage{
twostep_BinCont(
  X,
  Y,
  copula_family,
  marginal_surrogate,
  marginal_surrogate_estimator = NULL,
  method = "BFGS"
)
}
\arguments{
\item{X}{(numeric) Continuous surrogate variable}

\item{Y}{(integer) Binary true endpoint variable (\eqn{T_k \, \in \, \{0, 1\}})}

\item{copula_family}{Copula family, one of the following:
\itemize{
\item \code{"clayton"}
\item \code{"frank"}
\item \code{"gumbel"}
\item \code{"gaussian"}
}}

\item{marginal_surrogate}{Marginal distribution for the surrogate. For all
available options, see \code{?Surrogate::cdf_fun}.}

\item{marginal_surrogate_estimator}{Not yet implemented}

\item{method}{Optimization algorithm for maximizing the objective function.
For all options, see \code{?maxLik::maxLik}. Defaults to \code{"BFGRS"}.}
}
\value{
A list with three elements:
\itemize{
\item ml_fit: object of class \code{maxLik::maxLik} that contains the estimated copula
model.
\item marginal_S_dist: object of class \code{fitdistrplus::fitdist} that represents the
marginal surrogate distribution.
\item copula_family: string that indicates the copula family
}
}
\description{
The \code{twostep_BinCont()} function fits the copula (sub)model fir a continuous
surrogate and binary true endpoint with a two-step estimator. In the first
step, the marginal distribution parameters are estimated through maximum
likelihood. In the second step, the copula parameter is estimate while
holding the marginal distribution parameters fixed.
}
