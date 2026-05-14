# Process JAGS model outputs

Process JAGS model outputs

## Usage

``` r
processJAGSoutput(fit, DIC, params.omit, verbose = TRUE)
```

## Arguments

- fit:

  A JAGS state-space model fitted to CWAC data

- DIC:

  Logical stating whether DIC should be computed

- params.omit:

  A character vector with the name of the parameters that should not be
  processed.

- verbose:

  Logical. If TRUE (default), then several messages are displayed during
  processing.

## Value

A list with two elements: i) plot: a plot with summer and winter fitted
states, as well as the long-term trend, ii) data: the data used to
create the individual plots. This is useful for extracting the data used
by ggplot to render the plots (e.g. for exporting to the dashboard)
