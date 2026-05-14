# Prepare Google Earth Engine data for catchments

Prepare Google Earth Engine data for catchments

## Usage

``` r
prepGEECatchmData(sp_code, catchment, config, monitor = TRUE)
```

## Arguments

- sp_code:

  SAFRING reference number of the species we want to analyze.

- catchment:

  An sf object with the polygons defining the catchments to be
  annotated.

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

- monitor:

  Logical. If TRUE (default) monitoring printed messages produced by
  `rgee` will displayed. If FALSE, only high-level messages will be
  displayed.

## Details

It is assumed that there is an asset on
[`rgee::ee_get_assethome()`](https://r-spatial.github.io/rgee/reference/ee_get_assethome.html)
named 'quin_catchm' that has the polygons defining the quinary
catchments. Note that we add an extra year before the start of the
series. This is because summer waterbird populations should be affected
by the conditions in the previous year rather than by conditions in the
following year.
