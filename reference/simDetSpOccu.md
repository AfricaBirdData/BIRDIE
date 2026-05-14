# Simulate detections from spOccupancy fit

This function comes largely from the fitted.PGOcc.R from the
[spOccupancy](https://github.com/doserjef/spOccupancy) package

## Usage

``` r
simDetSpOccu(object)
```

## Arguments

- object:

  an spOccupancy fit

## Value

A list with posterior detection probabilities samples and posterior
detection predictions samples. The results are given in a long format
and as an attribute the indices of the sites the samples correspond to
