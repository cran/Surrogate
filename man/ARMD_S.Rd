\name{ARMD_S}
\alias{ARMD_S}
\docType{data}
\title{
Data of the Age-Related Macular Degeneration Study (Strobl workshop)
}
\description{
These are the data of a clinical trial involving patients suffering from age-related macular degeneration (ARMD), a condition that involves a progressive loss of vision. A total of \eqn{190} patients from \eqn{42} centers participated in the trial. Patients' visual acuity was assessed using standardized vision charts. There were two treatment conditions (placebo and interferon-\eqn{\alpha}). The potential surrogate endpoint is the change in the visual acuity at \eqn{24} weeks (\eqn{6} months) after starting treatment. The true endpoint is the change in the visual acuity at \eqn{52} weeks.}
\usage{data(ARMD_S)}
\format{
A \code{data.frame} with \eqn{190} observations on \eqn{7} variables.
  \describe{
    \item{\code{Center}}{The center in which the patient was treated.}
    \item{\code{Id}}{The Patient ID.}
    \item{\code{Diff24}}{The change in the visual acuity at \eqn{24} weeks after starting treatment. This endpoint is a potential surrogate for \code{Diff52}.}
    \item{\code{Diff52}}{The change in the visual acuity at \eqn{52} weeks after starting treatment. This outcome serves as the true endpoint.}
    \item{\code{Bin.Diff24}}{The dichotomized change in the visual acuity at \eqn{24} weeks after starting treatment (with \eqn{\leq 0} dichotomized as 0 and \eqn{>0} as 1).}
    \item{\code{Bin.Diff52}}{The dichotomized change in the visual acuity at \eqn{52} weeks after starting treatment (with \eqn{\leq 0} dichotomized as 0 and \eqn{>0} as 1).}
    \item{\code{Treat}}{The treatment indicator, coded as \eqn{-1} = placebo and \eqn{1} = interferon-\eqn{\alpha}.}

  }
}
\keyword{ARMD_S}
