# Set a standard path for a pipeline output file associated with a species

The datapipeline produces many output files. This function creates an
standard output file path for those files that are associated with a
particular species. There might be instances where we need name that
deviates from the standard. We need to handle these exceptions case by
case. Those files that are not associated with a particular species
follow different standards.

## Usage

``` r
setSpOutFilePath(prefix, config, years, sp_code, ext)
```

## Arguments

- prefix:

  A character string with the prefix that will appear at the beginning
  of the name. This is what distinguishes files within the species
  directory.

- config:

  A list with pipeline configuration parameters. See
  [`configPipeline`](https://africabirddata.github.io/BIRDIE/reference/configPipeline.md).

- years:

  Character string with years to include in the name.

- sp_code:

  SAFRING reference number of the species we want to analyze.

- ext:

  The extension of the output file. Note that we must add the trailing
  '.' to the extension.

## Value

A character string with the path to/for a pipeline output file

## Examples
