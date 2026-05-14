# Run diagnostics for occupancy model fit

This function performs basic checks for an occupancy model fit, such as
convergence of parameters and posterior predictive checks

## Usage

``` r
ppl_diagnose_occu(fit, ppc = NULL, data = NULL, sp_code, year)
```

## Arguments

- fit:

  A `spOccupancy` model fit to ABAP data.

- ppc:

  A list with the outputs from
  [diagnoseGofSpOccu](https://africabirddata.github.io/BIRDIE/reference/diagnoseGofSpOccu.md).

- data:

  A dataset with detection/non-detection data. Not used at the moment
  because spOccupancy (currently used) fit objects contain the data used
  to fit the model.

- sp_code:

  SAFRING code of the species to run the pipeline for

- year:

  The year the pipeline is run for

## Value

A data frame with Rhat values for the different parameters estimated by
the model is returned. If `ppc` is provided it will save some time
because posterior simulations can be time consuming. Sometimes the
posterior GOF test and posterior simulations might be available from
previous pipeline runs.
