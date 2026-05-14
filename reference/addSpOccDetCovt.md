# Add detection covariate to spOccupancy data list

Add detection covariate to spOccupancy data list

## Usage

``` r
addSpOccDetCovt(spOcc_data, covt_data, seasons = NULL)
```

## Arguments

- spOcc_data:

  An spOccupancy data list

- covt_data:

  A data frame with the ABAP visits used to create spOcc_data. The data
  frame must contain the columns 'Pentad', 'StartDate' and the covariate
  that we want to add to spOcc_data. If the data is multi-season, then
  we also need the variable that identifies the season. Only one
  covariate at a time is allowed at the moment.

- seasons:

  The name of the variable used to identify the seasons in multi-season
  data sets. Defaults to NULL.

## Value

An spOccupancy data list with the additional detection covariate.
