# Create pipeline event log

This function will create a customized log entry for a log file that
records the activity of the pipeline. More precisely, it will record the
time, species and outcome of data preparation, model fitting,
diagnostics and preparation.

## Usage

``` r
createLog(
  config,
  log_file = NULL,
  date_time = NULL,
  species = NULL,
  model = NULL,
  year = NA,
  data = NA,
  fit = NA,
  diagnose = NA,
  summary = NA,
  package = NA,
  notes = "",
  full_log = NULL
)
```

## Arguments

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)
  and
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md).

- log_file:

  Optional. A character string with the path to the file that should be
  updated. If NULL (default) a new log file will be created.

- date_time:

  Optional. A character string with the log date. If NULL (default) the
  time given by [`Sys.time`](https://rdrr.io/r/base/Sys.time.html),
  formatted as SAST will be used.

- species:

  SAFRING code of the species the log corresponds to.

- model:

  Type of model being processed. Either "occ" or "ssm".

- year:

  Year the model is fit for.

- data:

  Data preparation status. Defaults to NA.

- fit:

  model fitting status. Defaults to NA.

- diagnose:

  Model diagnostics status. Defaults to NA.

- summary:

  Model summary status. Defaults to NA.

- package:

  A character string with the name of the package that should be used to
  fit models. Currently, only `spOccupancy` for occupancy and `jagsUI`
  for state-space modelling are supported. maintained.

- notes:

  Any additional comments. Defaults to NA.

- full_log:

  Optionally, we can pass a full vector log, with all the above
  elements. In this case, all other arguments but `config` and `logfile`
  will be ignored.

## Value

A .csv will be saved on the reports directory. The location of this
directory is configured in the configuration object passed as `config`
at `config$output_dir`.
