# Diagnose convergence for JAGS state-space model

This function runs basic Rhat checks for a JAGS SSM

## Usage

``` r
diagnoseRhatJagsSsm(fit, sp_code, config)
```

## Arguments

- fit:

  A JAGS state-space model fitted to CWAC data

- sp_code:

  SAFRING reference number of the species we want to analyse.

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

## Value

A data frame with Rhat values for the different parameters estimated by
the model.
