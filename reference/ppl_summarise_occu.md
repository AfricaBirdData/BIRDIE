# Summarise predictions from occupancy model

Summarise predictions from occupancy model

## Usage

``` r
ppl_summarise_occu(fit, sp_code, year_sel, config, ...)
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
