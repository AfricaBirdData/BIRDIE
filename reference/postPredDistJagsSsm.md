# Prepare posterior predictive distribution for JAGS state-space model

This function prepares a data frame with the observed response and the
corresponding estimates obtained from a Bayesian state-space model
fitted with JAGS

## Usage

``` r
postPredDistJagsSsm(fit, data, obs_error = TRUE, nsamples)
```

## Arguments

- fit:

  A state-space model fitted with JAGS

- data:

  The data used to fit the model `fit`. This must be counts from CWAC
  data. See
  [`ppl_fit_ssm_model`](https://africabirddata.github.io/BIRDIE/reference/ppl_fit_ssm_model.md)

- obs_error:

  Logical. Indicates whether observation error should be included in the
  posterior simulations. Defaults to TRUE.

- nsamples:

  Number of posterior samples that should be used to build the data
  frame. Defaults to 500.

## Value

A data frame with the observed response and the corresponding estimates
obtained from a Bayesian state-space model fitted with JAGS. The
variable `obs_sim` correspond contains the real and simulated data. The
variable `iter` identifies the iteration each observation corresponds
to, and those observation with `iter = 0` and/or `obs = 1` correspond to
observed data.
