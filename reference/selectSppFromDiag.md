# Select species based on model diagnostics

The datapipeline produces one occupancy diagnostic file for each species
and year. This function combines all these files and returns those
species with convergence or goodnes-of-fit issues.

## Usage

``` r
selectSppFromDiag(config, sp_codes, year, module = c("dst", "abu"))
```

## Arguments

- config:

  A list with pipeline configuration parameters. See
  [configPipeline](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md).

- sp_codes:

  SAFRING reference numbers of the species we want diagnostics for.

- year:

  The year for which diagnostics are required.

- module:

  Either "dst" if we want occupancy diagnostics or "abu" if we want
  state-space model diagnostics.

## Value

A list with two elements: \$no_converge, contains codes of all species
with convergence issues, \$bad_fit, contains codes of species with
goodness-of-fit issues (small Bayesian p-value).

## Examples

``` r
if (FALSE) { # \dontrun{
config <- configPipeline(year = 2010,
                         dur = 3,
                         occ_mod = c("log_dist_coast", "elev", "log_hum.km2", "wetcon",
                                     "watrec", "watext", "log_watext", "watext:watrec",
                                     "ndvi", "prcp", "tdiff"),
                         det_mod = c("(1|obs_id)", "log_hours", "prcp", "tdiff", "cwac"),
                         fixed_vars = c("Pentad", "lon", "lat", "watocc_ever", "wetext_2018","wetcon_2018",
                                        "dist_coast", "elev"),
                         package = "spOccupancy",
                         data_dir = "analysis/hpc/imports",
                         out_dir = "analysis/hpc/imports",
                         server = TRUE)
sp_codes <- config$species

selectSppFromDiag(config, sp_codes, 2008)
} # }
```
