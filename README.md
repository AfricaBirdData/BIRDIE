# BIRDIE
This repository contains code and documentation related to the South Africa Biodiversity Data Pipeline for Wetlands and Waterbirds (BIRDIE) project: https://jrsbiodiversity.org/grants/sanbi-2020/

## Basic structure
The folders *R*, *data* and *man* contain functions, data and documentation for the BIRDIE package respectively. See the file *main_dev.R* and https://r-pkgs.org/, to get and idea of how to include new functions, data and modify existing ones. The folder *data_prep* stores scripts that are used to create data that is contained in the package, so that we can keep track of how data were produced.

The folder *analysis* contains data, scripts and outputs from analyses related to the BIRDIE project. This folder is not part of the package, so nothing in there is meant to be installed.

Data that should be included in and exported with the package should go to the *data* folder and it should be properly documented (see *main_dev.R* and https://r-pkgs.org/). Data is meant for running anaysis should go to the *analysis/data* folder. It is not necessary to document these data.
