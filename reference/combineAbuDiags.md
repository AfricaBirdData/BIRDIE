# Combine species abundance diagnostics into a single object

The data pipeline produces one abundance model diagnostic file for each
species and year. This function reads diagnostics files and combines
them into a single object

## Usage

``` r
combineAbuDiags(config, sp_codes, year)
```

## Arguments

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md).

- sp_codes:

  SAFRING reference numbers of the species we want diagnostics for.

- year:

  The year(s) for which diagnostics are required. Years should be
  formated as the two last digits of the beginning and last years. See

## Value

A dataframe with diagnostic paramters

## Examples

``` r
if (FALSE) { # \dontrun{
config <- configPipeline(
    year = 2021,
    dur = 29,
    mod_file = "cwac_ssm_two_season_mean_rev.R",
    package = "jagsUI",
    data_dir = NULL,
    out_dir = NULL,
    server = FALSE
)
sp_codes <- config$species

combineAbuDiags(config, sp_codes, "93_21")
} # }
```
