# BIRDIE
This repository contains code and documentation related to the South Africa Biodiversity Data Pipeline for Wetlands and Waterbirds (BIRDIE) project: https://jrsbiodiversity.org/grants/sanbi-2020/

## Basic structure
There are two basic elements in the BIRDIE repo: i) a package that contains functions and data that is installed as an extension to base R and ii) an analysis folder that contains analysis scripts and data that runs in R but that is not integrated in R. In other words, the package adds functionality and the analysis uses this functionality.

The folders *R*, *data* and *man* contain functions, data and documentation for the BIRDIE package respectively. See the file *main_dev.R* and https://r-pkgs.org/, to get and idea of how to include new functions, data and modify existing ones. The folder *data_prep* stores scripts that are used to create data that is contained in the package, so that we can keep track of how data were produced.

The folder *analysis* contains data, scripts and outputs from analyses related to the BIRDIE project. This folder is not part of the package, so nothing in there is meant to be installed.

Data that is meant to be included in and exported with the package should go to the *data* folder, and it should be properly documented (see *main_dev.R* and https://r-pkgs.org/). These data is usually meant for running examples or it may also be data that is used very frequently. Data that is only meant for running anaysis and not to exported with the package, should go to the *analysis/data* folder. It is not necessary to document these data.
