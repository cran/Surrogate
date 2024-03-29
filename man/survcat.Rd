% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survcat.R
\name{survcat}
\alias{survcat}
\title{Compute surrogacy measures for a categorical (ordinal) surrogate and a time-to-event true endpoint in the meta-analytic multiple-trial setting.}
\usage{
survcat(data, true, trueind, surrog, trt, center, trial, patientid)
}
\arguments{
\item{data}{A data frame with the correct columns (See details).}

\item{true}{Observed time-to-event (true endpoint).}

\item{trueind}{Time-to-event indicator.}

\item{surrog}{Ordinal surrogate endpoint.}

\item{trt}{Treatment indicator.}

\item{center}{Center indicator (equal to trial if there are no different centers).}

\item{trial}{Trial indicator.}

\item{patientid}{Patient indicator.}
}
\value{
Returns an object of class "survcat" that can be used to evaluate surrogacy and contains the following elements:
\itemize{
\item Indiv.GlobalOdds: a data frame that contains the Global Odds and 95\% confidence interval to evaluate surrogacy at the individual level.
\item Trial.R2: a data frame that contains the correlation coefficient and 95\% confidence interval to evaluate surrogacy at the trial level.
\item EstTreatEffects: a data frame that contains the estimated treatment effects and sample size for each trial.
}
}
\description{
The function 'survcat()' fits the model for a categorical (ordinal) surrogate and time-to-event true endpoint developed by Burzykowski et al. (2004) in the meta-analytic multiple-trial setting.
}
\section{Model}{
In the model developed by Burzykowski et al. (2004), a copula-based model is used for the true endpoint and a latent continuous variable, underlying the surrogate endpoint.
More specifically, the Plackett copula is used. The marginal model for the surrogate endpoint is a proportional odds model. For the true endpoint, the proportional hazards model is used.
The quality of the surrogate at the individual level can be evaluated by using the copula parameter Theta, which takes the form of a global odds ratio.
The quality of the surrogate at the trial level can be evaluated by considering the correlation coefficient between the estimated treatment effects, while adjusting for the estimation error.
}

\section{Data Format}{
The data frame must contains the following columns:
\itemize{
\item a column with the observed time-to-event (true endpoint)
\item a column with the time-to-event indicator: 1 if true event is observed, 0 otherwise
\item a column with the ordinal surrogate endpoint: 1 2 3 ... K
\item a column with the treatment indicator: 0 or 1
\item a column with the trial indicator
\item a column with the center indicator. If there are no different centers within each trial, the center indicator is equal to the trial indicator
\item a column with the patient indicator
}
}

\examples{
\dontrun{
data("colorectal4")
fit <- survcat(data = colorectal4, true = truend, trueind = trueind, surrog = surrogend,
               trt = treatn, center = center, trial = trialend, patientid = patid)
print(fit)
summary(fit)
plot(fit)
}
}
\references{
Burzykowski, T., Molenberghs, G., & Buyse, M. (2004). The validation of surrogate end points by using data from randomized clinical trials: a case-study in advanced colorectal cancer. Journal of the Royal Statistical Society Series A: Statistics in Society, 167(1), 103-124.
}
\author{
Dries De Witte
}
