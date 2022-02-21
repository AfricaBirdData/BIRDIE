
# In this script we prepare VISIT data for fitting a spatial occupancy model
# We will use Google Earth Engine to annotate the pentads with covariates

# Note that visit data is the same for all species except for the detection
# column. Therefore, covariates will be the same for any species.

library(rgee)
library(ABAP)
library(dplyr)
library(sf)

# Initialize Earth Engine
ee_check()
ee_Initialize(drive = TRUE)

rm(list = ls())


config <- BIRDIE::configPreambOccuR(year = 2009, server = TRUE) # The year doesn't matter

years <- 2008:2019

pentads_sa <- ABAP::getRegionPentads(.region_type = "country", .region = "South Africa")

for(i in seq_along(years)){
#for(i in 7:12){

    # Download SABAP data for any species -------------------------------------

    year <- years[i]

    visit <- ABAP::getAbapData(.spp_code = 6,
                               .region_type = "country",
                               .region = "South Africa",
                               .years = year)


    # Annotate data with NDVI values ------------------------------------------

    # We will use the NDVI values closer to the date of the visit

    # Make spatial object and select relevant columns
    visit <- visit %>%
        left_join(pentads_sa,
                  by = c("Pentad" = "Name")) %>%
        st_sf() %>%
        filter(!st_is_empty(.)) %>%     # Remove rows without geometry
        mutate(Date = as.character(StartDate)) %>%   # GEE doesn't like dates
        dplyr::select(CardNo, StartDate, Date, Pentad, TotalHours)

    # Upload to GEE
    ee_visit <- visit %>%
        dplyr::select(-c(StartDate, TotalHours)) %>%
        sf_as_ee(via = "getInfo")

    # Annotate with GEE TerraClimate
    visit_new <- addVarEEclosestImage(ee_pentads = ee_visit,
                                      collection = "MODIS/006/MOD13A2",
                                      reducer = "mean",                          # We only need spatial reducer
                                      maxdiff = 15,                              # This is the maximum time difference that GEE checks
                                      bands = c("NDVI"))

    visit <- visit %>%
        st_drop_geometry() %>%
        left_join(visit_new %>%
                      st_drop_geometry() %>%
                      dplyr::select(CardNo, val) %>%
                      rename(ndvi = val),
                  by = c("CardNo"))

    # Update
    if(file.exists(file.path(config$data_dir, "visit_dat_sa_gee_08_19.rds"))){
        visit_old <- readRDS(file.path(config$data_dir, "visit_dat_sa_gee_08_19.rds"))
        visit <- rbind(visit_old, visit)
    }

    # Save
    saveRDS(visit, file.path(config$data_dir, "visit_dat_sa_gee_08_19.rds"))

}


# Prepare variables for fitting -------------------------------------------

visit <- readRDS(file.path(config$data_dir, "visit_dat_sa_gee_08_19.rds"))

visit <- visit %>%
    mutate(Date = lubridate::date(Date),
           year = lubridate::year(Date),
           month = lubridate::month(Date),
           ndvi = ndvi/1e4)

saveRDS(visit, file.path(config$data_dir, "visit_dat_sa_gee_08_19.rds"))
