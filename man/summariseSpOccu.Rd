% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-spOccupancy.R
\name{summariseSpOccu}
\alias{summariseSpOccu}
\title{Summarise spOccupancy occupancy estimates}
\usage{
summariseSpOccu(pred_p, pred_psi, quants)
}
\arguments{
\item{pred_p}{Detection probabilities estimated from a spOccupancy model.
It must be a matrix (or mcmc object) with each row corresponding to an MCMC sample and each column
corresponding to a visit. This object must contain a "pentad" attribute indicating
which pentad each column correspond to, an attribute "year" indicating what
year the probabilities correspond to, and an attribute "obs" indicating whether
the species was detected in the visit or not. Outputs from \code{\link{predictSpOccu}},
should be readily appropriate.}

\item{pred_psi}{Occupancy probabilities estimated from an a spOccupancy model.
It must be a matrix (or mcmc object) with each row corresponding to an MCMC sample and each column
corresponding to a pentad. This object must contain a "pentad" attribute indicating
which pentad each column correspond to, and an attribute "year" indicating what
year the probabilities correspond to. Outputs from \code{\link{predictSpOccu}},
should be readily appropriate.}

\item{quants}{Quantiles to summarise predictions distribution passed as c("lower", "med", "upper").}
}
\value{
A tibble with estimates and/or quantiles for each pentad in site_data:
\itemize{
\item psi: occupancy probability,
\item p: detection probability,
\item occu: realized occupancy (occupancy conditional on data).
}
}
\description{
spOccupancy produces estimates of detection probability for each
visit and occupancy probabilities for each site. This functions takes these
predictions and creates detection probabilities, occupancy probabilities and
realized occupancy probabilities (occupancy probabilities conditional on
observed data) for each site.
}
