% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survcat.R
\name{summary.survcat}
\alias{summary.survcat}
\title{Provides a summary of the surrogacy measures for an object fitted with the 'survcat()' function.}
\usage{
\method{summary}{survcat}(object, ...)
}
\arguments{
\item{object}{An object of class 'survcat' fitted with the 'survcat()' function.}

\item{...}{...}
}
\value{
The surrogacy measures with their 95\% confidence intervals.
}
\description{
Provides a summary of the surrogacy measures for an object fitted with the 'survcat()' function.
}
\examples{
\dontrun{
data("colorectal4")
fit <- survcat(data = colorectal4, true = truend, trueind = trueind, surrog = surrogend,
               trt = treatn, center = center, trial = trialend, patientid = patid)
summary(fit)
}
}
