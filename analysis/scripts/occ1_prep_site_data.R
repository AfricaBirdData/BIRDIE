
# In this script we prepare SITE data for fitting a spatial occupancy model
# to a pilot species in South Africa

library(BIRDIE)
library(sf)
library(raster)
library(tidyverse)

rm(list = ls())


# Select a region --------------------------------------------------------

# region <- "North West"
region <- "South Africa"



# Extract region pentads --------------------------------------------------

# Run only if pentads_sa is not already available
# sf::sf_use_s2(FALSE) # s2 intersection takes very long
# pentads_sel <- getRegionPentads(.country = region, .path = "analysis/downloads")
#
# plot(st_geometry(pentads_sel), col = "blue", lwd = 0.1)
#
# saveRDS(pentads_sel, "analysis/data/pentads_sa.rds")


# Prepare site covariate data --------------------------------------------

pentads_sel <- readRDS("analysis/data/pentads_sa.rds")

# Define climatic covariates
covts <- c("prcp", "tmax", "tmin", "aet", "pet")

future::plan("multisession", workers = 6)

for(i in seq_along(covts)){
    pentads_sel <- addOccSiteCovt(pentads_sel,
                                  covt = covts[i],
                                  years = 2008:2011,
                                  covts_dir = "analysis/downloads/",
                                  file_fix = c("terraClim_", "_03_19"))
}

future::plan("sequential")

saveRDS(pentads_sel, "analysis/data/site_dat_sa_wcovts.rds")


# Prepare surface water ---------------------------------------------------

# Load data
pentads_sel <- readRDS("analysis/data/site_dat_sa_wcovts.rds")

# Extract pentads
pentads_sel$id <- 1:nrow(pentads_sel)

# Define coordinates of surface water maps
cc <- list(c(10, 20),
           c(20, 20),
           c(30, 20),
           c(10, 30),
           c(20, 30),
           c(30, 30))

out <- vector("list",length = length(cc))

future::plan("multisession", workers = 6)

for(i in seq_along(cc)){

    # Load water raster
    water <- raster(paste0("analysis/downloads/surf_water_", cc[[i]][1], "E_", cc[[i]][2], "S.tif"))

    water <- raster::crop(water, pentads_sel)

    # Extract values
    pentads_water <- BIRDIE::exactExtractParll(water, pentads_sel,
                                               ncores = future::nbrOfWorkers(),
                                               fun = "mean")

    pentads_water <- pentads_water %>%
        rename(water = vals)

    out[[i]] <- pentads_water %>%
        filter(!is.na(water))

    pentads_sel <- pentads_water %>%
        filter(is.na(water)) %>%
        dplyr::select(-water)

}

future::plan("sequential")

pentads_water <- do.call("rbind", out)
pentads_water <- arrange(pentads_water, id)

plot(pentads_water["water"])

# Save
pentads_sel <- pentads_water

saveRDS(pentads_sel, "analysis/data/site_dat_sa_wcovts.rds")


# Add centroids -----------------------------------------------------------

pentads_sel <- readRDS("analysis/data/site_dat_sa_wcovts.rds")

# Add centroid coordinates
sf::sf_use_s2(TRUE)
pentads_sel <- pentads_sel %>%
    dplyr::mutate(lon = sf::st_coordinates(sf::st_centroid(.))[,1],
                  lat = sf::st_coordinates(sf::st_centroid(.))[,2])

saveRDS(pentads_sel, "analysis/data/site_dat_sa_wcovts.rds")
