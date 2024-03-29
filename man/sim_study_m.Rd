% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim-utils.R
\name{sim_study_m}
\alias{sim_study_m}
\title{sim_study_m}
\usage{
sim_study_m(mus, nt, nc, r)
}
\arguments{
\item{nt}{number indicating the number of participants for the treated group.}

\item{nc}{number indicating the number of participants for the control group. If \code{NULL} the value assigned to \code{nt} will be used.}

\item{r}{number indicating the correlation between sampling errors.}

\item{thetas}{a numeric vector with population level effect sizes.}
}
\value{
a dataframe
}
\description{
Simulate a single study with multiple outcomes. The function simulate a compound symmetry structure where the correlation between sampling errors is a single value and the sampling variance is fixed to 1.
}
\examples{
sim_study_m(mus = c(0.1, 0.5, 1), nt = 30, nc = 30, r = 0.7)
}
