# Diagnose convergence for spOccupancy model

This function runs basic Rhat checks for an spOccupancy model

## Usage

``` r
diagnoseRhatSpOccu(fit, sp_code, year)
```

## Arguments

- fit:

  A `spOccupancy` model fit to ABAP data.

- sp_code:

  SAFRING reference number of the species we want to analyse.

- year:

  The year the pipeline is run for

## Value

A data frame with Rhat values for the different parameters estimated by
the model.
