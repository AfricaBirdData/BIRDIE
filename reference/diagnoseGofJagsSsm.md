# Diagnose goodness-of-fit for JAGS state-space model

This function runs basic posterior predictive checks for a JAGS SSM

## Usage

``` r
diagnoseGofJagsSsm(fit, counts, linear = TRUE)
```

## Arguments

- fit:

  A JAGS state-space model fitted to CWAC data

- counts:

  A dataframe with CWAC counts ready for model fit. See
  [`ppl_create_data_ssm`](https://africabirddata.github.io/BIRDIE/reference/ppl_create_data_ssm.md)

- linear:

  Logical. If TRUE (default) then posterior predictive checks are made
  on the data scale, otherwise they are done in the model (log) scale.

## Value

A data frame with Rhat values for the different parameters estimated by
the model.
