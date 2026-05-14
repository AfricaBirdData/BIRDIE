# Configure pipeline parameters

Set basic variables to run the BIRDIE pipeline locally or remotely.

## Usage

``` r
configPipeline(
  year,
  dur,
  region = c("southafrica", "kenya"),
  module = c("dst", "abu"),
  occ_mod = NULL,
  det_mod = NULL,
  fixed_vars = NULL,
  mod_file = NULL,
  server = FALSE,
  data_dir = NULL,
  mod_dir = NULL,
  out_dir = NULL,
  package = NULL
)
```

## Arguments

- year:

  Year of interest.

- dur:

  Temporal coverage of the analysis in years. `year` will be the last
  year covered by the analysis.

- region:

  A character string with the region we want to run the pipeline for.
  Currently only "South Africa" and "Kenya" are supported.

- module:

  A character string defining the module the pipeline should run. At the
  moment it can be one of `c("dst", "abu")` for distributions and
  abundance respectively.

- occ_mod:

  A character vector with the names of the variables to include in the
  occupancy process in the occupancy model. Random effects and
  interactions are specified as in
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html). Note that only
  second order interactions are accepted at the moment (i.e.,
  interactions of two variables).

- det_mod:

  A character vector with the names of the variables to include in the
  detection process in the occupancy model. Random effects and
  interactions are specified as in
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html). Note that only
  second order interactions are accepted at the moment. (i.e.,
  interactions of two variables).

- fixed_vars:

  A character vector with the names of the variables included in the
  occupancy model that don't change over time.

- mod_file:

  Name of file containing model, with out path to directory. Directory
  is specified in `mod_dir`. This is typically used for JAGS or Stan
  where models are written on an external file.

- server:

  Logical. If TRUE the preamble is prepared to run remotely, otherwise
  it is prepared to run locally.

- data_dir:

  Path to data directory. There are a few inputs to the pipeline that it
  doesn't generate itself, such as some environmental layers that are
  not on Google Earth Engine. Those would be stored here.

- mod_dir:

  Path to directory where models are saved.

- out_dir:

  Path to output directory. Pipeline outputs will be stored here,
  including intermediate outputs, so most what we need is here.

- package:

  A character string with the name of the package that should be used
  for fitting occupancy or state-space models.

## Value

A list of parameters that will be passed on to other functions in the
pipeline.

## Examples

``` r
config <- configPipeline(
    year = 2021,
    dur = 29,
    mod_file = "cwac_ssm_two_season_mean_rev.R",
    package = "jagsUI",
    data_dir = NULL,
    out_dir = NULL,
    server = FALSE
    )
```
