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


# Install -----------------------------------------------------------------

install()


