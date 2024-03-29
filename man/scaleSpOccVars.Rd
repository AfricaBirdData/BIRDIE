% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-spOccupancy.R
\name{scaleSpOccVars}
\alias{scaleSpOccVars}
\title{Scale covariates in spOccupancy-type data}
\usage{
scaleSpOccVars(spOcc_data, var_type, scale_vars)
}
\arguments{
\item{spOcc_data}{an spOccupancy data list.}

\item{var_type}{Type of variables we want to scale. Currently, one of "occ"
occupancy covariates, "det" detection covariates.}

\item{scale_vars}{A vector with the names of the covariates that we want to
scale.}
}
\value{
An spOccupancy data list with the scaled covariates, substituting the
original, unscaled covariates. A single factor is used to center and scale all
data (across all dimensions) of the covariate. That means, for example, that
all seasons will be scaled by the same amount. The factors used
for centering and scaling are stored as attributes of each covariate.
}
\description{
Scale covariates in spOccupancy-type data
}
