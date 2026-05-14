# Prepare spOccupancy data for multi-season model fitting

Prepare spOccupancy data for multi-season model fitting

## Usage

``` r
prepSpOccuData_multi(sp_code, year, config, spatial = FALSE, ...)
```

## Arguments

- sp_code:

  SAFRING code of the species to run the pipeline for

- year:

  Year to run to the pipeline for

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

- spatial:

  Logical, indicating whether spatial random effects should be included
  in the model. Defaults to FALSE.

- ...:

  Other arguments passed on to other functions
