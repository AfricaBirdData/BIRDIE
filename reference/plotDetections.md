# Plot detection from occupancy data

Plot detection from occupancy data

## Usage

``` r
plotDetections(site_data, visit_data)
```

## Arguments

- site_data:

  An sf object with the sites visited in a given period.

- visit_data:

  A data frame with the visits that have occurred in that a given
  period. visit_data must have the same name as site_data for the sites
  and the variable must be called 'site'. Detections must be represented
  by a variable named 'obs', which should be 1 when the species was
  detected and 0 otherwise.
