# Predict from spOccupancy model fit

Predict from spOccupancy model fit

## Usage

``` r
predictSpOccu(fit, sp_code, year_sel, config, ...)
```

## Arguments

- fit:

  An occupancy model fit to summarise occupancy and detection
  probabilities from.

- sp_code:

  SAFRING code of the species to run the pipeline for

- year_sel:

  Year in data to run model for.

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

- ...:

  Other arguments passed on to other functions

## Value

A list with two elements: 1) posterior occupancy probability samples for
the South African pentads, and 2) posterior detection probability
samples for each visit in the ABAP data for the year `year_sel`.
