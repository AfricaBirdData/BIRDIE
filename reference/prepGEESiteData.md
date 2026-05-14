# Prepare Google Earth Engine site data

Prepare Google Earth Engine site data

## Usage

``` r
prepGEESiteData(config, pentads, asset_id, upload_asset = TRUE, monitor = TRUE)
```

## Arguments

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

- pentads:

  An sf spatial object with the set of pentads that we want to annotate
  with environmental covariates. See
  [`getRegionPentads`](https://rdrr.io/pkg/ABAP/man/getRegionPentads.html)

- asset_id:

  Character string with the name given to the object created in Google
  Earth Engine (asset) that contains the sites in `pentads`.

- upload_asset:

  If TRUE (default), the object `pentads` will be uploaded to Google
  Earth Engine and an asset under the name of `asset_id` will be
  created. If FALSE, it will be assumed that an asset named after
  `asset_id` is already present in GEE and `pentads` will not be
  uploaded.

- monitor:

  Logical. If TRUE (default) monitoring printed messages produced by
  `rgee` will displayed. If FALSE, only high-level messages will be
  displayed.

## Details

Note that some GEE layers don't have information past a certain date. At
the time of writing surface water layers only have information up until
2021 and human population density up until 2020. We have set up the code
in such a way that data past the last date of the layer get annotated
with the latest available information. Take this into consideration for
the analyses. Code should be updated as more information becomes
available.
