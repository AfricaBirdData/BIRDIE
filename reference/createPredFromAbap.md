# Create prediction object directly from ABAP data

This function is used when it is not possible to run an occupancy model
for an species and year. Then, the raw data is presented but they must
have the same structure as the prediction object to be able to integrate
them in the database.

## Usage

``` r
createPredFromAbap(sp_code, year_sel, config)
```

## Arguments

- sp_code:

  SAFRING code of the species to run the pipeline for

- year_sel:

  Year in data to run model for.

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)
