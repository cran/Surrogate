% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ITCI-Dvine-inference-utils.R
\name{ICA_given_model_constructor}
\alias{ICA_given_model_constructor}
\title{Constructor for the function that returns that ICA as a function of the
identifiable parameters}
\usage{
ICA_given_model_constructor(
  fitted_model,
  copula_par_unid,
  copula_family2,
  rotation_par_unid,
  n_prec,
  measure = "ICA",
  mutinfo_estimator = NULL,
  ICA_estimator = NULL,
  seed,
  composite = NULL,
  restr_time = +Inf
)
}
\arguments{
\item{fitted_model}{Returned value from \code{\link[=fit_copula_OrdOrd]{fit_copula_OrdOrd()}},
\code{\link[=fit_copula_OrdCont]{fit_copula_OrdCont()}}, or \code{\link[=fit_copula_ContCont]{fit_copula_ContCont()}}. This object
contains the estimated identifiable part of the joint distribution for the
potential outcomes.}

\item{copula_par_unid}{Parameter vector for the sequence of \emph{unidentifiable}
bivariate copulas that define the D-vine copula. The elements of
\code{copula_par} correspond to \eqn{(c_{23}, c_{13;2}, c_{24;3}, c_{14;23})}.}

\item{copula_family2}{Copula family of the other bivariate copulas. For the
possible options, see \code{loglik_copula_scale()}. The elements of
\code{copula_family2} correspond to \eqn{(c_{23}, c_{13;2}, c_{24;3}, c_{14;23})}.}

\item{rotation_par_unid}{Vector of rotation parameters for the sequence of
\emph{unidentifiable} bivariate copulas that define the D-vine copula. The elements of
\code{rotation_par} correspond to \eqn{(c_{23}, c_{13;2},
  c_{24;3}, c_{14;23})}.}

\item{n_prec}{Number of Monte Carlo samples for the computation of the mutual
information.}

\item{measure}{Compute intervals for which measure of surrogacy? Defaults to
\code{"ICA"}. See first column names of \code{sens_results} for other possibilities.}

\item{mutinfo_estimator}{Function that estimates the mutual information
between the first two arguments which are numeric vectors. Defaults to
\code{FNN::mutinfo()} with default arguments in the survival-survival setting. This
argument is not used for non-survival-survival settings.}

\item{ICA_estimator}{Function that estimates the ICA between the first two
arguments which are numeric vectors. Defaults to \code{NULL} which corresponds
to using \code{\link[=estimate_ICA_ContCont]{estimate_ICA_ContCont()}}, \code{\link[=estimate_ICA_OrdCont]{estimate_ICA_OrdCont()}}, or
\code{\link[=estimate_ICA_OrdOrd]{estimate_ICA_OrdOrd()}} (depending on the endpoint types). This argument is
not used in the survival-survival setting.}

\item{seed}{Seed for Monte Carlo sampling. This seed does not affect the global
environment.}

\item{composite}{(boolean)  If \code{composite} is \code{TRUE}, then the surrogate
endpoint is a composite of both a "pure" surrogate endpoint and the true
endpoint, e.g., progression-free survival is the minimum of time-to-progression
and time-to-death.}

\item{restr_time}{Restriction time for the potential outcomes. Defaults to
\code{+Inf} which means no restriction. Otherwise, the sampled potential outcomes
are replace by \code{pmin(S0, restr_time)} (and similarly for the other potential
outcomes).}
}
\value{
A function that computes the ICA as a function of the identifiable
parameters. In this computation, the unidentifiable parameters are fixed at
the values supplied as arguments to \code{\link[=ICA_given_model_constructor_SurvSurv]{ICA_given_model_constructor_SurvSurv()}} or
\code{\link[=ICA_given_model_constructor]{ICA_given_model_constructor()}}.
}
\description{
\code{\link[=ICA_given_model_constructor]{ICA_given_model_constructor()}} returns a function fixes the unidentifiable
parameters at user-specified values and takes the identifiable parameters as
argument.
}
