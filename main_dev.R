# 12-5-2021

# Script to keep track of package development

# DO NOT RUN EVERYTHING AGAIN

# If you want to add a function copy and paste a function block
# at the end of the script (before the install part) and run
# that section that you just pasted.

rm(list = ls())

library(devtools)
library(testthat)


# Create a package structure
# create_package("D:/Documentos/mis_trabajos/Academic/BIRDIE/BIRDIE")

# Add license
# use_mit_license("BIRDIE Development Team")

# Remember to edit the DESCRIPTION file



# Imports -----------------------------------------------------------------

# Import pipe
use_pipe()

# Import packages
use_package("dplyr")
use_package("tidyr")
use_package("lubridate")
use_package("jagsUI")
use_package("ggplot2")
use_package("gridExtra")
use_package("CWAC")
use_package("sf")
use_package("lwgeom")
use_package("raster")
use_package("rgeos")
use_package("utils")
use_package("grDevices")
use_package("mgcv")
use_package("stocc")


# Function fitCwacSsm2ss --------------------------------------------------

# Add function
use_r("fitCwacSsm2ss")

# test locally
load_all()

fitCwacSsm2ss(counts, mod_file = "analysis/models/cwac_ssm_2ss_fxd.jags",
              param = c("beta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))

# Add documentation
# Add ROxygen skeleton manually
document()

check()

# Add tests

use_testthat()

use_test()

test()


# Function prepSsmData --------------------------------------------------

# Add function
use_r("prepSsmData")

# test locally
load_all()

counts <- getCwacSiteCounts(26352535)
prepSsmData(counts)

# Add documentation
# Add ROxygen skeleton manually
document()

check()

# Add tests

use_testthat()

use_test()

test()


# Function plotSsm --------------------------------------------------

# Add function
use_r("plotSsm")

# test locally
load_all()

counts <- getCwacSiteCounts(26352535)
prepSsmData(counts)

# Add documentation
# Add ROxygen skeleton manually
document()

check()

# Add tests

use_testthat()

use_test()

test()


# Data barberspan -------------------------------------------------------------

# Download Barberspan data and save as package data
# source("data_prep/barberspan_prep.R")

# Create an ROxygen2 file and document
document()


# Function getRegionPentads --------------------------------------------------

# Add function
use_r("getRegionPentads")

# test locally
load_all()

getRegionPentads(country = "South Africa", province = "North West")
getRegionPentads(country = "South Africa", province = "North West", path = "analysis/data")

# Add documentation
# Add ROxygen skeleton manually
document()

check()

# Add tests

use_testthat()

use_test()

test()


# Data pentads_sabap2 -------------------------------------------------------------

# Prepare pentads
# source("data_prep/barberspan_prep.R")

# Create an ROxygen2 file and document
document()


# Function loopSsmAllSpp --------------------------------------------------

# Add function
use_r("loopSsmAllSpp")

# test locally
load_all()

counts <- barberspan %>%
    dplyr::filter(spp %in% c(6, 41, 54))

loopSsmAllSpp(counts,
              data_outdir = "analysis/output/dashboard_out/",
              plot_outdir = "analysis/output/dashboard_out/",
              param = c("beta", "lambda", "sig.zeta", "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"),
              jags_control = list(ncores = 3))

# Add documentation
# Add ROxygen skeleton manually
document()

check()

# Add tests

use_testthat()

use_test()

test()


# Function writeJagsModelFile --------------------------------------------------

# Add function
use_r("writeJagsModelFile")

# test locally
load_all()

modpath <- writeJagsModelFile()

# Add documentation
# Add ROxygen skeleton manually
document()

check()

# Add tests

use_testthat()

use_test()

test()


# Function findNextIndex --------------------------------------------------

# Add function
use_r("findNextIndex")

# test locally
load_all()

x <- rbinom(20, 1, 0.2)
findNextIndex(x, 1)

# Add documentation
# Add ROxygen skeleton manually
document()

check()

# Add tests

use_testthat()

use_test()

test()


# Function prepCtSsmData --------------------------------------------------

# Add function
use_r("prepCtSsmData")

# test locally
load_all()

counts <- barberspan
prepCtSsmData(counts)
prepCtSsmData(counts, species = 212)

# Add documentation
# Add ROxygen skeleton manually
document()

check()

# Add tests

use_testthat()

use_test()

test()


# Function prepOccuData --------------------------------------------------

# Add function
use_r("prepOccuData")

# test locally
load_all()



# Add documentation
# Add ROxygen skeleton manually
document()

check()

# Add tests

use_testthat()

use_test()

test()

# Install -----------------------------------------------------------------

devtools::install()


