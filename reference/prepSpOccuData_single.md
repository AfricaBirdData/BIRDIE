# Prepare spOccupancy data for single-season model fitting

Prepare spOccupancy data for single-season model fitting

## Usage

``` r
prepSpOccuData_single(
  site_data,
  visit_data,
  config,
  spatial = FALSE,
  sp_sites = NULL
)
```

## Arguments

- site_data:

  A data frame containing information about covariates associated with
  ABAP pentads for a given year.

- visit_data:

  A data frame with information associated with sampling visits to ABAP
  pentads in a given year. Detection/non-detection data for the species
  of interest must also be included in this data frame.

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

- spatial:

  Logical, indicating whether spatial random effects should be included
  in the model. Defaults to FALSE.

- sp_sites:

  Spatial object with the sites to be used for fitting spatial models
