# Plot occupancy covariates

Plots selected covariates against detection/non-detection data

## Usage

``` r
plotOccuVars(occu_data, vars, type = "ind")
```

## Arguments

- occu_data:

  An list with occupancy site and visit data. See
  [`ppl_create_site_visit`](https://africabirddata.github.io/BIRDIE/reference/ppl_create_site_visit.md)

- vars:

  Character vector with the names of the variables to plot. They must
  match those names in occu_data.

- type:

  The type of plot we want to display. Either "ind" which shows
  individual plots for each variable, or "cross", which shows a
  scatterplot with two variables. We need to pass on exactly two
  variables.

## Value

A plot with detection/non-detection against selected variables
