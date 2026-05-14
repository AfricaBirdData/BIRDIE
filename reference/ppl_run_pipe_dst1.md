# Run distribution indicators pipeline module 1

Module 1 of the distribution indicators pipeline estimates occupancy
probabilities in South Africa for a selected species. These occupancy
probabilities form the basis for building more elaborated indicators in
other pipeline modules.

## Usage

``` r
ppl_run_pipe_dst1(
  sp_code,
  year,
  config,
  steps = c("data", "fit", "diagnose", "summary"),
  ...
)
```

## Arguments

- sp_code:

  SAFRING code of the species to run the pipeline for

- year:

  Year to run to the pipeline for

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

- steps:

  A character vector containing the steps of the pipeline to run. Can
  contain: "data", "fit", "diagnose", "summary". Defaults to all of
  them.

- ...:

  Other arguments passed on to other functions
