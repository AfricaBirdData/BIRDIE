# Annotate site data with Google Earth Engine with a time limit

Annotate site data with Google Earth Engine with a time limit

## Usage

``` r
addSiteVarEETimeLimit(
  ee_feats,
  ee_collection,
  band,
  last_year,
  reducer,
  unmask,
  monitor,
  config
)
```

## Arguments

- ee_feats:

  A feature collection with the sites we want to annotate. We should
  have uploaded an sf object with the sites to GEE, previously. See
  `uploadFeaturesToEE`

- ee_collection:

  A GEE collection produced with `ee$ImageCollection()`. See [GEE
  catalog](https://developers.google.com/earth-engine/datasets/catalog).

- band:

  The band of the collection we are going to use to annotate our data

- last_year:

  The last year provided by the GEE image collection

- reducer:

  The reducer we will use to summarize the images of the image
  collection

- unmask:

  GEE masks missing values, which means they are not used for computing
  means, counts, etc. Sometimes we might want to avoid this behaviour
  and use 0 instead of NA. If so, set unmask to TRUE.

- monitor:

  Logical. If TRUE (default) monitoring messages produced by `rgee` will
  displayed. If FALSE, only high-level messages will be displayed.

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

## Details

The function annotates data with the corresponding year if it is present
in the GEE image collection. Years after the last year present in
'config\$years' are annotated with the last year present in the image
collection.
