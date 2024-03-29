% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-spOccupancy.R
\name{diagnoseGofSpOccu}
\alias{diagnoseGofSpOccu}
\title{Conduct goodness-of-fit test for spOccupancy models}
\usage{
diagnoseGofSpOccu(object, post_sims, fit_stat = "chi-squared", group = 1)
}
\arguments{
\item{object}{an spOccupancy fit}

\item{post_sims}{A list of posterior simulations obtained with simDetSpOccu().}

\item{fit_stat}{Goodness of fit statistic to compute. Currently only Chi-squared
is supported.}

\item{group}{Whether the GOF statistic should be computed for sites (1) or
for visits (2). Currently, only grouping by site is supported}
}
\value{
A list GOF statistics computed for observed data and for simulated
data. One statistic for is computed for each MCMC iteration.
The function will create a list with several objects. The most important are:
\itemize{
\item fit.y: For each MCMC sample, chi-squared statistic for total number of detections
for all sites in data wrt expected model expectation
\item fit.y.rep: For each MCMC sample, chi-squared statistic for total number of detections
for all sites in simulated data wrt expected model expectation
\item fit.y.group.quants: For each site, chi-squared statistic for total number of detections
in data wrt expected model expectation (posterior distribution quantiles)
\item fit.y.rep.group.quants: For each site, chi-squared statistic for total number of detections
in simulated data wrt expected model expectation (posterior distribution quantiles)
\item y.summ.per.site: Data frame with number of detections per site obtained from
data simulated from the posterior distribution. Also the observed number of detections
in the data.
}
}
\description{
This function comes largely from the ppcOcc.R from the
\href{https://github.com/doserjef/spOccupancy}{spOccupancy} package
}
