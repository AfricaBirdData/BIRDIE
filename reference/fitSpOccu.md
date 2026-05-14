# Fit spOccupancy model

Fit spOccupancy model

## Usage

``` r
fitSpOccu(
  site_data_year,
  visit_data_year,
  config,
  sp_code,
  spatial = FALSE,
  sp_sites,
  ...
)
```

## Arguments

- site_data_year:

  Occupancy site data for a single year and species (see
  [`ppl_create_site_visit`](https://africabirddata.github.io/BIRDIE/reference/ppl_create_site_visit.md))

- visit_data_year:

  Occupancy visit data for a single year and species (see
  [`ppl_create_site_visit`](https://africabirddata.github.io/BIRDIE/reference/ppl_create_site_visit.md))

- config:

  A list with pipeline configuration parameters (see
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)).

- sp_code:

  SAFRING code of the species the pipeline is running for

- spatial:

  Logical, indicating whether spatial random effects should be included
  in the model (TRUE) or not (FALSE, default)

- sp_sites:

  Spatial object containing the pentads in `site_data_year`.

- ...:

  Other arguments that might be needed (e.g. for messages)

## Value

Either a spOccupancy model fit or the integer 3, indicating that model
fit failed.
