% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-spOccupancy.R
\name{prepSpOccuData_multi}
\alias{prepSpOccuData_multi}
\title{Prepare spOccupancy data for multi-season model fitting}
\usage{
prepSpOccuData_multi(sp_code, year, config, spatial = FALSE, ...)
}
\arguments{
\item{sp_code}{SAFRING code of the species to run the pipeline for}

\item{year}{Year to run to the pipeline for}

\item{config}{A list with pipeline configuration parameters.
See \code{\link{configPipeline}}}

\item{spatial}{Logical, indicating whether spatial random effects should be
included in the model. Defaults to FALSE.}

\item{...}{Other arguments passed on to other functions}
}
\description{
Prepare spOccupancy data for multi-season model fitting
}
