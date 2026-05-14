# Run diagnostics for a state-space model fit

This function performs basic checks for an SSM fit, such as convergence
of parameters and posterior predictive checks

## Usage

``` r
ppl_diagnose_ssm(fit, counts, sp_code, config)
```

## Arguments

- fit:

  A JAGS state-space model fitted to CWAC data

- counts:

  A dataframe with CWAC counts ready for model fit. See
  [`ppl_create_data_ssm`](https://africabirddata.github.io/BIRDIE/reference/ppl_create_data_ssm.md)

- sp_code:

  SAFRING reference number of the species we want to analyse.

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

## Value

A list with Rhat values and posterior check statistics is returned. At
them moment, we obtain Rhat values for all monitored parameters and
three posterior check statistics: "Tmean" proportion of posterior
simulations with mean greater than that observed in the data (we would
like values close to 0.5), "Tsd" proportion of posterior simulations
with sd greater than that observed in the data (we would like values
close to 0.5), "Tdiff" mean difference between observed data and
posterior simulations (we would like values close to 0).
