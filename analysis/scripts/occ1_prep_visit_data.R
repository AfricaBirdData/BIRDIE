
# In this script we prepare VISIT data for fitting a spatial occupancy model
# We will use Google Earth Engine to annotate the pentads with covariates

# Note that visit data is the same for all species except for the detection
# column. Therefore, covariates will be the same for any species.

library(rgee)
library(SABAP)
library(dplyr)
library(sf)

# Initialize Earth Engine
ee_check()
ee_Initialize(drive = TRUE)

rm(list = ls())

years <- 2009:2019

for(i in seq_along(years)){

    # Download SABAP data for any species -------------------------------------

    year <- years[i]

    visit <- SABAP::getSabapData(.spp_code = 6,
                                 .region_type = "country",
                                 .region = "South Africa",
                                 .years = year)


    # Annotate data with NDVI values ------------------------------------------

    # We will use the NDVI values closer to the date of the visit

    # Make spatial object and select relevant columns
    visit <- visit %>%
        left_join(SABAP::pentads_sabap , by = c("Pentad" = "Name")) %>%
        st_sf() %>%
        filter(!st_is_empty(.)) %>%     # Remove rows without geometry
        mutate(Date = as.character(StartDate)) %>%   # GEE doesn't like dates
        dplyr::select(CardNo, Date, Pentad, TotalHours)

    # Upload to GEE
    ee_visit <- visit %>%
        dplyr::select(-TotalHours) %>%
        sf_as_ee(via = "getInfo")

    # Annotate with GEE TerraClimate
    visit_new <- addVarEEclosestImage(ee_pentads = ee_visit,
                                      collection = "MODIS/006/MOD13A2",
                                      reducer = "mean",                          # We only need spatial reducer
                                      maxdiff = 15,                              # This is the maximum time difference that GEE checks
                                      bands = c("NDVI"))

    visit <- visit %>%
        left_join(visit_new %>%
                      st_drop_geometry() %>%
                      dplyr::select(CardNo, val) %>%
                      rename(ndvi = val),
                  by = c("CardNo"))

    # Update
    if(file.exists(paste0("analysis/data/visit_dat_gee_08_19.rds"))){
        visit_old <- readRDS(paste0("analysis/data/visit_dat_gee_08_19.rds"))
        visit <- rbind(visit_old, visit)
    }

    # Save
    saveRDS(visit, paste0("analysis/data/visit_dat_gee_08_19.rds"))

}
