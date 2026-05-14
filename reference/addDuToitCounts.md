# Add DuToit's Pan data from Doug

This function is used to attach data from DuToit's Pan that Doug
Harebottle provided and that are not on the CWAC database.

## Usage

``` r
addDuToitCounts(counts, config)
```

## Arguments

- counts:

  A dataframe with CWAC data coming out of
  [`getCwacSppCounts`](https://rdrr.io/pkg/CWAC/man/getCwacSppCounts.html)
  SAFRING reference number of the species we want to add data to.

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md).

## Value

A dataframe with DuToit's counts added to `counts`

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

site_code <- 28462448

counts <- CWAC::getCwacSiteCounts(site_code)

complete_counts <- addDuToitCounts(counts, config)


} # }
```
