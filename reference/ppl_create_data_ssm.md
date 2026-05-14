# Prepare covariate data for abundance pipeline

Prepare covariate data for abundance pipeline

## Usage

``` r
ppl_create_data_ssm(
  sp_code,
  year,
  catchment,
  config,
  force_gee = TRUE,
  force_gee_upload = TRUE,
  steps = c("subset", "missing", "gee", "model"),
  ...
)
```

## Arguments

- sp_code:

  SAFRING reference number of the species we want to analyse.

- year:

  Year for which the SSM data should be prepared.

- catchment:

  A sf object with polygons corresponding to the catchment, or reference
  area considered for the covariates associated with CWAC sites.

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

- force_gee:

  Logical. If TRUE (default), then count data will be annotated with
  environmental information from GEE. If FALSE, then we assume data have
  already been annotated and only some manipulation to join catchment
  environmental data with count data is required.

- force_gee_upload:

  Logical. If TRUE (default), then the catchment polygons will be
  uploaded to GEE under the name 'quin_catchm'. If FALSE, then we assume
  these polygons are already \#' present in GEE server and we don't need
  to upload them again.

- steps:

  A character vector expressing the processing steps for the CWAC data.
  It can be one or more of: "missing" - add missing counts as missing
  data, "gee" - annotate data with covariates from Google Earth Engine,
  "subset" - subset data to those sites with presence of the species and
  a coverage of at least 10 years from 1993 to 2021, and "model" prepare
  data for model fitting by adding some aux variables and ordering the
  data. It defaults to all steps.

- ...:

  Other arguments passed on to
  [prepGEECatchmData](https://africabirddata.github.io/BIRDIE/reference/prepGEECatchmData.md)
