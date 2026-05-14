# Create site and visit occupancy files

This function prepares site and visit occupancy data to fit an occupancy
model from ABAP data. It has two parts: the first part downloads ABAP
data and annotates them with covariates from Google Earth Engine using
the functions
[`prepGEESiteData`](https://africabirddata.github.io/BIRDIE/reference/prepGEESiteData.md)
and
[`prepGEEVisitData`](https://africabirddata.github.io/BIRDIE/reference/prepGEEVisitData.md),
the second part uses the function
[`createOccuData`](https://africabirddata.github.io/BIRDIE/reference/createOccuData.md)
to format the data.

## Usage

``` r
ppl_create_site_visit(
  config,
  sp_code,
  force_gee_dwld = FALSE,
  force_site_visit = FALSE,
  force_abap_dwld = FALSE,
  monitor_gee = TRUE
)
```

## Arguments

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

- sp_code:

  SAFRING code of the species to run the pipeline for

- force_gee_dwld:

  Whether covariates from Google Earth Engine should be downloaded, even
  if a file with covariates is already present on disk. Defaults to
  FALSE.

- force_site_visit:

  Whether site and visit data should be prepared even if visit and site
  data files are already on disk. Defaults to FALSE

- force_abap_dwld:

  Indicates whether ABAP data must be downloaded for the species and
  years indicated by 'sp_code' and 'years'. If TRUE, data will be
  downloaded from ABAP once per session and cached in a temp file. After
  this the cached file will be used, unless download is set to FALSE, in
  which case data will be downloaded regardless of the cached file.

- monitor_gee:

  if TRUE (default) periodic messages of the state of the downloads from
  GEE will be printed on screen.

## Value

The first part of the function creates two data frames (in .csv format)
that will be saved to disk: GEE annotated ABAP site data and GEE
annotated ABAP visit data. The second part of the functions creates
three data frames that will be saved to disk: site, visit and species
detection data frames, all in .csv format.
