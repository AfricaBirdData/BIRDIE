
# In this script we prepare VISIT data for fitting a spatial occupancy model
# to a pilot species in South Africa

library(BIRDIE)
library(sf)
library(raster)
library(tidyverse)

rm(list = ls())


# Select a species and a region -------------------------------------------

sp_sel <- 6

# region <- "North West"
region <- "South Africa"

# Load site info
sites <- readRDS("analysis/data/site_dat_sa_wcovts.rds")


# Prepare occupancy visit data -------------------------------------------

# Prepare SABAP data and pentads for species 6. (There should be an option to supply pentads object)
pa_dat <- getOccVisitData(region_type = "country", region = region, pentads = sites,
                          species = sp_sel, path = "analysis/downloads")

# Subset to years of interest
pa_dat <- pa_dat %>%
    dplyr::filter(year > 2007, year < 2012)

saveRDS(pa_dat, paste0("analysis/data/visit_dat_", sp_sel, ".rds"))

pa_dat <- readRDS(paste0("analysis/data/visit_dat_", sp_sel, ".rds"))

# Annotate occupancy visit data with climatic covariates
covts = c("prcp", "tmax", "tmin", "aet", "pet")

future::plan("multisession", workers = 6) # Set up multicore if needed

for(i in seq_along(covts)){
    pa_dat <- addOccVisitCovt(visits = pa_dat,
                                     sites = sites,
                                     covt = covts[i],
                                     covts_dir = "analysis/downloads/",
                                     file_fix = c("terraClim_", "_03_19"))
}

future::plan("sequential")

saveRDS(pa_dat, paste0("analysis/data/visit_dat_", sp_sel, "_wcovts.rds"))


# Prepare non-linear terms ------------------------------------------------

# Load detection/non-detection data
pa_dat <- readRDS(paste0("analysis/data/visit_dat_", sp_sel, "_wcovts.rds"))

# Add cyclic basis values for time of year to visit data
spl_bs <- mgcv::cSplineDes(x = pa_dat$month,
                           knots = seq(1, 12, length.out = 5),
                           ord = 4, derivs = 0)
colnames(spl_bs) <- paste0("toy.", 1:ncol(spl_bs))
pa_dat <- cbind(pa_dat, spl_bs)

saveRDS(pa_dat, paste0("analysis/data/visit_dat_", sp_sel, "_wcovts.rds"))


