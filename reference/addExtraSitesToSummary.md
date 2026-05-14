# Add raw counts to CWAC SSM predictions

We only run state-space models for those CWAC sites that have enough
data. See
[`ppl_create_data_ssm`](https://africabirddata.github.io/BIRDIE/reference/ppl_create_data_ssm.md).
For those sites that don't have enough data, we present the raw data.
This function is used to bind the output of
[`ppl_summarise_ssm`](https://africabirddata.github.io/BIRDIE/reference/ppl_summarise_ssm.md)
with the raw data from those sites that didn't have enough data to run
an analysis.

## Usage

``` r
addExtraSitesToSummary(counts, preds)
```

## Arguments

- counts:

  A dataframe with CWAC counts. It is preferable that the dataframe
  containd missing counts as well. See `addMissingCwacCounts`

- preds:

  A dataframe with estimates from a state-space model fitted to CWAC
  data. See
  [`ppl_summarise_ssm`](https://africabirddata.github.io/BIRDIE/reference/ppl_summarise_ssm.md)

## Value

A dataframe with predictions for sites with good data, together with raw
counts for sites with not enough data to run an analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
sp_code <- 87
config <- configPipeline(
    year = 2022,
    dur = 30,
    module = "abu",
    mod_file = "cwac_ssm_two_season_mean_rev_jump.R",
    package = "jagsUI",
    data_dir = NULL,     # this might have to be adapted?
    out_dir = NULL,     # this might have to be adapted?
    server = FALSE
)
counts <- read.csv(setSpOutFilePath("cwac_data_w_miss", config, config$years_ch, sp_code, ".csv"))
preds <- setSpOutFilePath("ssm_pred", config, config$years_ch, sp_code, "_all.csv")
preds_w_raw <- addExtraSitesToSummary(counts, preds)
} # }
```
