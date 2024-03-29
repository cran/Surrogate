% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survbin.R
\name{summary.survbin}
\alias{summary.survbin}
\title{Provides a summary of the surrogacy measures for an object fitted with the 'survbin()' function.}
\usage{
\method{summary}{survbin}(object, ...)
}
\arguments{
\item{object}{An object of class 'survbin' fitted with the 'survbin()' function.}

\item{...}{...}
}
\value{
The surrogacy measures with their 95\% confidence intervals.
}
\description{
Provides a summary of the surrogacy measures for an object fitted with the 'survbin()' function.
}
\examples{
\dontrun{
data("colorectal")
fit_bin <- survbin(data = colorectal, true = surv, trueind = SURVIND, surrog = responder,
                   trt = TREAT, center = CENTER, trial = TRIAL, patientid = patientid)
summary(fit_bin)
}
}
