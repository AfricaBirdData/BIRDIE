
# In this script we prepare data for fitting a spatial occupancy model
# to a pilot species in the North West province of South Africa

library(BIRDIE)
library(sf)
library(tidyverse)

rm(list = ls())


# Select a species and a region -------------------------------------------

sp_sel <- 6

region <- "North West"


# Prepare occupancy visit data -------------------------------------------

# Prepare SABAP data and pentads for species 6 in the North West province.
pa_dat <- getOccVisitData(region_type = "province", region = region,
                          species = sp_sel, path = "analysis/data")

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


# Prepare occupancy site data -------------------------------------------

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


# Prepare non-linear terms ------------------------------------------------

# Load detection/non-detection data
pa_dat <- readRDS("analysis/out_nosync/pa_dat_6_wcovts_nw.rds")

# Add cyclic basis values for time of year to visit data
spl_bs <- mgcv::cSplineDes(x = pa_dat$visits$month,
                           knots = seq(1, 12, length.out = 5),
                           ord = 4, derivs = 0)
colnames(spl_bs) <- paste0("toy.", 1:ncol(spl_bs))
pa_dat$visits <- cbind(pa_dat$visits, spl_bs)

saveRDS(pa_dat, "analysis/data/pa_dat_6_wcovts_nw.rds")
