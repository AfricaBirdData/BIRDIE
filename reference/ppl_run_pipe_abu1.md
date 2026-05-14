# Run abundance pipeline ABU1

Run abundance pipeline ABU1

## Usage

``` r
ppl_run_pipe_abu1(
  sp_code,
  config,
  steps = c("data", "fit", "diagnose", "summary"),
  prep_data_steps,
  summary_scale = c("linear", "model"),
  ...
)
```

## Arguments

- sp_code:

  SAFRING reference number of the species we want to analyse.

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

- steps:

  Pipeline steps to run. It can be one or more of: c("data", "fit",
  "diagnose", "summary").

- prep_data_steps:

  Data preparation steps to pass on to
  [`ppl_create_data_ssm`](https://africabirddata.github.io/BIRDIE/reference/ppl_create_data_ssm.md)

- summary_scale:

  Either "linear", summaries are given the linear scale or "model",
  summaries are given in the modelling scale (typically log scale)

- ...:

  Other parameters to pass on to
  [`prepGEECatchmData`](https://africabirddata.github.io/BIRDIE/reference/prepGEECatchmData.md)

## Value

This function will run the whole abundance module of the BIRDIE pipeline
