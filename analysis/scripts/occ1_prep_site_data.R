
# In this script we prepare SITE data for fitting a spatial occupancy model

# We will use Google Earth Engine to annotate the pentads with covariates

library(rgee)
library(SABAP)
library(dplyr)
library(sf)
# library(BIRDIE)

# Initialize Earth Engine
ee_check()
ee_Initialize(drive = TRUE)

rm(list = ls())


# Load pentads to GEE -----------------------------------------------------

# Load SABAP pentads
pentads <- readRDS("analysis/data/pentads_sa.rds")

# Load to EE (if not done already)
assetId <- sprintf("%s/%s", ee_get_assethome(), 'pentads')
# pentads %>%
#     sf_as_ee(assetId = assetId,
#              via = "getInfo_to_asset")

ee_pentads <- ee$FeatureCollection(assetId)


# Annotate with TerraClimate ----------------------------------------------

# Define bands
bands <- c("pr", "tmmn", "tmmx")

# Define function
f <- function(band){

    stackCollection <- EEcollectionToMultiband(collection = "IDAHO_EPSCOR/TERRACLIMATE",
                                               dates = c("2008-01-01", "2020-01-01"),
                                               band = band,
                                               group_type = "year",
                                               groups = 2008:2019,
                                               reducer = "mean")

    # Extract values from all bands of the image
    out <- addVarEEimage(ee_pentads, stackCollection, "mean")

    out <- out %>%
        rename_with(~gsub("X", band, .x), .cols = starts_with("X"))

    return(out)

}

pentad_vars <- purrr::map(bands, ~f(.x))

out <- Reduce("cbind", pentad_vars)

out <- out %>%
    dplyr::select(!contains("."))

saveRDS(out, "analysis/data/site_dat_gee_08_19.rds")


# Annotate with yearly surface water occurrence --------------------------------

# Find mean water presence for each pixel and year
band <- "waterClass"
stackCollection <- EEcollectionToMultiband(collection = "JRC/GSW1_3/YearlyHistory",
                                           dates = c("2008-01-01", "2020-01-01"),
                                           band = band,
                                           group_type = "year",
                                           groups = 2008:2019,
                                           unmask = TRUE)

# Find mean (mean) water presence for each pentad and year
out <- addVarEEimage(ee_pentads, stackCollection, "mean")

out <- out %>%
    rename_with(~gsub("X", band, .x), .cols = starts_with("X"))

plot(out['waterClass2011'], lwd = 0.01)

readRDS("analysis/data/site_dat_gee_08_19.rds") %>%
    left_join(out %>%
                  st_drop_geometry() %>%
                  dplyr::select(-id),
              by = "Name") %>%
    saveRDS("analysis/data/site_dat_gee_08_19.rds")


# Annotate with overall surface water occurrence --------------------------

out <- addVarEEimage(ee_pentads = ee_pentads,
                     image = "JRC/GSW1_3/GlobalSurfaceWater",
                     reducer = "mean",
                     bands = "occurrence",
                     unmask = TRUE)

out <- out %>%
    rename(water_occur = occurrence)

readRDS("analysis/data/site_dat_gee_08_19.rds") %>%
    left_join(out %>%
                  st_drop_geometry() %>%
                  dplyr::select(-id),
              by = "Name") %>%
    saveRDS("analysis/data/site_dat_gee_08_19.rds")


# Annotate with yearly NDVI -----------------------------------------------

# Find mean water presence for each pixel and year
band <- "NDVI"
stackCollection <- EEcollectionToMultiband(collection = "MODIS/006/MOD13A2",
                                           dates = c("2008-01-01", "2020-01-01"),
                                           band = band,
                                           group_type = "year",
                                           groups = 2008:2019,
                                           reducer = "mean",
                                           unmask = FALSE)

# Find mean (mean) water presence for each pentad and year
out <- addVarEEimage(ee_pentads, stackCollection, "mean")

out <- out %>%
    rename_with(~gsub("X", band, .x), .cols = starts_with("X"))

plot(out['NDVI2011'], lwd = 0.01)

readRDS("analysis/data/site_dat_gee_08_19.rds") %>%
    left_join(out %>%
                  st_drop_geometry() %>%
                  dplyr::select(-id),
              by = "Name") %>%
    saveRDS("analysis/data/site_dat_gee_08_19.rds")
