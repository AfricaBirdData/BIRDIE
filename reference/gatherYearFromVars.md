# Gather years from variables with year

Gather years from variables with year

## Usage

``` r
gatherYearFromVars(x, vars, sep)
```

## Arguments

- x:

  A dataframe with columns that represent a variable on a given year.
  Note that there must be a unique key that identifies each case
  uniquely.

- vars:

  Variables to separate in variable and year (in different columns).

- sep:

  A character that identifies the separator between the name of the
  variable and the year.

## Value

A dataframe with separate columns for variables and years. There will be
a column that represents the year and several columns that represent the
variables measured across years.

## Examples

``` r
df <- data.frame(id = 1:10,
                 v1_2010 = rnorm(10, 0, 1),
                 v1_2011 = rnorm(10, 100, 10),
                 v2_2010 = runif(10, 0, 1),
                 v2_2011 = runif(10, 100, 110))

gatherYearFromVars(x = df, vars = names(df)[-1], sep = "_")
#> # A tibble: 20 × 4
#>       id  year        v1       v2
#>    <int> <int>     <dbl>    <dbl>
#>  1     1  2010  -1.40      0.680 
#>  2     1  2011  94.5     105.    
#>  3     2  2010   0.255     0.499 
#>  4     2  2011 106.      103.    
#>  5     3  2010  -2.44      0.642 
#>  6     3  2011 121.      102.    
#>  7     4  2010  -0.00557   0.660 
#>  8     4  2011  83.7     105.    
#>  9     5  2010   0.622     0.0960
#> 10     5  2011 105.      105.    
#> 11     6  2010   1.15      0.766 
#> 12     6  2011  81.4     108.    
#> 13     7  2010  -1.82      0.770 
#> 14     7  2011  94.8     102.    
#> 15     8  2010  -0.247     0.991 
#> 16     8  2011  99.5     107.    
#> 17     9  2010  -0.244     0.971 
#> 18     9  2011 105.      101.    
#> 19    10  2010  -0.283     0.389 
#> 20    10  2011  90.9     104.    
```
