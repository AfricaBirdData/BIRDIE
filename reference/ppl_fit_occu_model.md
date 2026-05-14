# Fit occupancy model

Fit occupancy model

## Usage

``` r
ppl_fit_occu_model(sp_code, year_sel, config, spatial = FALSE, ...)
```

## Arguments

- sp_code:

  SAFRING code of the species to run the pipeline for

- year_sel:

  Year in data to run model for.

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

- spatial:

  Whether a spatial model should be fit. Defaults to FALSE.

- ...:

  Other arguments passed on to other functions

## Value

It can return either a number correspoding to the status of the fitting
process: 1 - No detections, 2 - too few detections, 3 - model fitting
failed (inherited from
[`fitSpOccu`](https://africabirddata.github.io/BIRDIE/reference/fitSpOccu.md))
