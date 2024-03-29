% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDetections.R
\name{plotDetections}
\alias{plotDetections}
\title{Plot detection from occupancy data}
\usage{
plotDetections(site_data, visit_data)
}
\arguments{
\item{site_data}{An sf object with the sites visited in a given period.}

\item{visit_data}{A data frame with the visits that have occurred in that a
given period. visit_data must have the same name as site_data for the sites
and the variable must be called 'site'. Detections must be represented by a
variable named 'obs', which should be 1 when the species was detected and 0
otherwise.}
}
\description{
Plot detection from occupancy data
}
