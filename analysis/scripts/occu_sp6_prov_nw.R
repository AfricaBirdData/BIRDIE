
# In this script we prepare data for fitting a spatial occupancy model
# to a pilot species in the North West province of South Africa

library(BIRDIE)
library(sf)
library(tidyverse)

rm(list = ls())


# Prepare data occupancy visit data ----------------------------------------

# Prepare SABAP data and pentads for species 6 in the North West province.
pa_dat <- getOccVisitData(region_type = "province", region = "North West",
                          species = 6, path = "analysis/data")

saveRDS(pa_dat, "analysis/out_nosync/pa_dat_6_nw.rds")

pa_dat <- readRDS("analysis/out_nosync/pa_dat_6_nw.rds")

# Annotate occupancy visit data with climatic covariates
covts = c("prcp", "tmax", "tmin", "aet", "pet")

future::plan("multisession", workers = 6) # Set up multicore if needed

for(i in seq_along(covts)){
    pa_dat$visits <- addOccVisitCovt(visits = pa_dat$visits,
                                     sites = pa_dat$sites,
                                     covt = covts[i],
                                     covts_dir = "analysis/out_nosync/",
                                     file_fix = c("terraClim_", "_03_19_nw"))
}

future::plan("sequential")

saveRDS(pa_dat, "analysis/out_nosync/pa_dat_6_wcovts_nw.rds")


# Prepare site covariates -------------------------------------------------

# Load detection/non-detection data
pa_dat <- readRDS("analysis/out_nosync/pa_dat_6_wcovts_nw.rds")

# Define covariates
covts <- c("prcp", "tmax", "tmin", "aet", "pet")

future::plan("multisession", workers = 6)

for(i in seq_along(covts)){
    pa_dat$sites <- addOccSiteCovt(pa_dat$sites,
                                   covt = covts[i],
                                   years = 2007:2010,
                                   covts_dir = "analysis/out_nosync/",
                                   file_fix = c("terraClim_", "_03_19_nw"))
}

future::plan("sequential")

saveRDS(pa_dat, "analysis/out_nosync/pa_dat_6_wcovts_nw.rds")

