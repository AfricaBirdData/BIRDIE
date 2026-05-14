# Plot state series of a JAGS state-space model with two seasons

Plot state series of a JAGS state-space model with two seasons

## Usage

``` r
plotJagsSsm2ss(
  fit,
  ssm_counts,
  linear = TRUE,
  plot_options = list(colors = NULL, pers_theme = NULL)
)
```

## Arguments

- fit:

  A JAGS state-space model fitted to CWAC data

- ssm_counts:

  A data frame with the count data use to fit the state-space model

- linear:

  If TRUE (default) abundance estimates and data are transformed back to
  its original scale.

- plot_options:

  A list with two elements: colors - the colours of the points that will
  appear in the plot (two values), and pers_theme - A personalized
  ggplot theme.

## Value

A list with two elements: i) plot: a plot with summer and winter fitted
states, as well as the long-term trend, ii) data: the data used to
create the individual plots. This is useful for extracting the data used
by ggplot to render the plots (e.g. for exporting to the dashboard)

## Examples

``` r
if (FALSE) { # \dontrun{
counts <- barberspan
ssmcounts <- prepSsmData(counts, species = NULL)
fit <- fitCwacSsm(ssmcounts, mod_file = "mymodel.jags",
param = c("beta", "lambda", "sig.zeta",
"sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))
plotSsm2ss(fit = fit, ssm_counts = ssmcounts)
} # }
```
