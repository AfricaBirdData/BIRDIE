# Prepare Google Earth Engine visit data

Prepare Google Earth Engine visit data

## Usage

``` r
prepGEEVisitData(config, visits, asset_id, upload_asset = TRUE, monitor = TRUE)
```

## Arguments

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

- visits:

  An sf spatial object with the visit data that we want to annotate with
  environmental covariates. Note that we need these data to be in a
  spatial format, so we should probably join them with pentad data. See
  [`getRegionPentads`](https://rdrr.io/pkg/ABAP/man/getRegionPentads.html)

- asset_id:

  Character string with the name given to the object created in Google
  Earth Engine (asset) that contains the sites in `visits`.

- upload_asset:

  If TRUE (default), the object `visits` will be uploaded to Google
  Earth Engine and an asset under the name of `asset_id` will be
  created. If FALSE, it will be assumed that an asset named after
  `asset_id` is already present in GEE and `visits` will not be
  uploaded.

- monitor:

  Logical. If TRUE (default) monitoring printed messages produced by
  `rgee` will displayed. If FALSE, only high-level messages will be
  displayed.

## Details

Note that some GEE layers don't have information past a certain date. At
the time of writing surface water layers only have information up until
2021 and human population density up until 2020. We have set up the code
in such a way that visits past the last date of the layer get annotated
with the latest available information. Take this into consideration for
the analyses. Code should be updated as more information becomes
available.
