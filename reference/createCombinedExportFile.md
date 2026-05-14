# Create a combined file for export to data mart

Create a combined file for export to data mart

## Usage

``` r
createCombinedExportFile(config, type = c("abu", "dst"))
```

## Arguments

- config:

  A list with pipeline configuration parameters. See
  [configPipeline](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md)

- type:

  A character string with three options: "abu", for abundance estimates
  files, or "dst" for occupancy estimates files.

## Value

It will create a file in `config$out_dir/exports` combining all files of
the selected type of all species and years in `config$years`.
