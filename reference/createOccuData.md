# Create data for fitting an occupancy model

This function takes generic occupancy site data and, optionally, visit
data without species observations, adds detection data for a species in
a given year, runs some checks and formats the output.

## Usage

``` r
createOccuData(config, sp_code, years, site_data, visit_data = NULL)
```

## Arguments

- config:

  A list of configuration parameters see
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

- sp_code:

  SAFRING_No of the species of interest as extracted from ABAP. Ignored
  if download set to FALSE.

- years:

  A numeric vector with elements corresponding to the years we want data
  for. Ignored if download set to FALSE.

- site_data:

  A dataframe with occupancy site data.

- visit_data:

  Optional. A dataframe with occupancy visit data.

## Value

A list containing two data frames: one with site data and one with visit
data, if visit data is provided. The second element will be NULL if
visit data not provided. Note that if `visit_data` is provided then all
the sites without visits will be removed - i.e., both site and visit
data frames will have the same set of sites.
