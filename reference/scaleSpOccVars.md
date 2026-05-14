# Scale covariates in spOccupancy-type data

Scale covariates in spOccupancy-type data

## Usage

``` r
scaleSpOccVars(spOcc_data, var_type, scale_vars)
```

## Arguments

- spOcc_data:

  an spOccupancy data list.

- var_type:

  Type of variables we want to scale. Currently, one of "occ" occupancy
  covariates, "det" detection covariates.

- scale_vars:

  A vector with the names of the covariates that we want to scale.

## Value

An spOccupancy data list with the scaled covariates, substituting the
original, unscaled covariates. A single factor is used to center and
scale all data (across all dimensions) of the covariate. That means, for
example, that all seasons will be scaled by the same amount. The factors
used for centering and scaling are stored as attributes of each
covariate.
