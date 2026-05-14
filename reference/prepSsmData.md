# Prepare CWAC data to fit a state-space model

Prepare CWAC data to fit a state-space model

## Usage

``` r
prepSsmData(counts, spp_sel = NULL, keep = NULL)
```

## Arguments

- counts:

  A data frame with raw CWAC count data

- spp_sel:

  An optional vector of species codes. Only records for these species
  are prepared for fitting an SSM, others are discarded.

- keep:

  A vector of variables to keep after the processing. Useful, if there
  are covariates of interest. If NULL, only year, season, start date,
  count and species name are returned.

## Value

A tibble with clean and prepared data for fitting an SSM (e.g. filled
with missing years)

## Examples

``` r
counts <- barberspan
prepSsmData(counts)
#> # A tibble: 34 × 5
#>     year Season StartDate  count spp  
#>    <dbl> <chr>  <date>     <int> <chr>
#>  1  1993 S      1993-01-23 11637 multi
#>  2  1993 W      1993-07-10  4113 multi
#>  3  1994 S      1994-01-15  2141 multi
#>  4  1994 W      1994-07-16  7191 multi
#>  5  1995 S      1995-01-28 16567 multi
#>  6  1996 S      1996-01-01  2440 multi
#>  7  1999 S      1999-01-16 11932 multi
#>  8  1999 W      1999-09-12 22112 multi
#>  9  2000 S      2000-01-29  1674 multi
#> 10  2000 W      2000-07-15  8649 multi
#> # ℹ 24 more rows
prepSsmData(counts, spp_sel = 212)
#> # A tibble: 34 × 5
#>     year Season StartDate  count spp             
#>    <dbl> <chr>  <date>     <int> <chr>           
#>  1  1993 S      1993-01-23  6996 Red-knobbed Coot
#>  2  1993 W      1993-07-10  1451 Red-knobbed Coot
#>  3  1994 S      1994-01-15   536 Red-knobbed Coot
#>  4  1994 W      1994-07-16  2883 Red-knobbed Coot
#>  5  1995 S      1995-01-28 14826 Red-knobbed Coot
#>  6  1996 S      1996-01-01   559 Red-knobbed Coot
#>  7  1999 S      1999-01-16  8981 Red-knobbed Coot
#>  8  1999 W      1999-09-12 15628 Red-knobbed Coot
#>  9  2000 S      2000-01-29   462 Red-knobbed Coot
#> 10  2000 W      2000-07-15  4940 Red-knobbed Coot
#> # ℹ 24 more rows
prepSsmData(counts, spp_sel = c(212, 50))
#> # A tibble: 34 × 5
#>     year Season StartDate  count spp  
#>    <dbl> <chr>  <date>     <int> <chr>
#>  1  1993 S      1993-01-23  7002 multi
#>  2  1993 W      1993-07-10  1451 multi
#>  3  1994 S      1994-01-15   560 multi
#>  4  1994 W      1994-07-16  2921 multi
#>  5  1995 S      1995-01-28 14860 multi
#>  6  1996 S      1996-01-01   680 multi
#>  7  1999 S      1999-01-16  9321 multi
#>  8  1999 W      1999-09-12 15774 multi
#>  9  2000 S      2000-01-29   655 multi
#> 10  2000 W      2000-07-15  5571 multi
#> # ℹ 24 more rows
```
