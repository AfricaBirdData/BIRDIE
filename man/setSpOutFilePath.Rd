% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{setSpOutFilePath}
\alias{setSpOutFilePath}
\title{Set a standard path for a pipeline output file associated with a species}
\usage{
setSpOutFilePath(prefix, config, years, sp_code, ext)
}
\arguments{
\item{prefix}{A character string with the prefix that will appear at the
beginning of the name. This is what distinguishes files within the species
directory.}

\item{config}{A list with pipeline configuration parameters.
See \code{\link{configPipeline}}.}

\item{years}{Character string with years to include in the name.}

\item{sp_code}{SAFRING reference number of the species we want to analyze.}

\item{ext}{The extension of the output file. Note that we must add the
trailing '.' to the extension.}
}
\value{
A character string with the path to/for a pipeline output file
}
\description{
The datapipeline produces many output files. This function
creates an standard output file path for those files that are associated with
a particular species. There might be instances where we need name that
deviates from the standard. We need to handle these exceptions
case by case. Those files that are not associated with a particular species
follow different standards.
}
\examples{
\dontrun{
config <- configAbu(2017, server = FALSE)
sp_code <- 4

# Set a file path for a state-space model output
setSpOutFilePath("ssm_fit", config, config$years_ch, sp_code, ".rds)

# Set a file path for model-ready data
setSpOutFilePath("model_data", config, sp_code, ".csv)
}
}
