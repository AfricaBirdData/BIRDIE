
# In this script we prepare SITE data for fitting a spatial occupancy model
# We will use Google Earth Engine to annotate the pentads with covariates

library(rgee)
library(ABAP)
library(dplyr)
library(sf)

# Initialize Earth Engine
ee_check()
ee_Initialize(drive = TRUE)

rm(list = ls())


config <- BIRDIE::configPreambOccuR(year = 2009, server = TRUE) # The year doesn't matter


# Load pentads to GEE -----------------------------------------------------

# Load ABAP pentads
pentads_sa <- getRegionPentads(.region_type = "country", .region = "South Africa")

# Set a name for the asset
eeid <- sprintf("%s/%s", ee_get_assethome(), 'pentads')

# Upload to EE (if not done already)
uploadPentadsToEE(pentads = dplyr::select(pentads_sa, Name),
                  asset_id = eeid,
                  load = FALSE)

# Load pentads from GEE
ee_pentads <- ee$FeatureCollection(eeid)


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
                                               reducer = "mean",
                                               unmask = FALSE)

    # Extract values from all bands of the image
    out <- addVarEEimage(ee_pentads, stackCollection, "mean")

    out <- out %>%
        rename_with(~gsub("X", band, .x), .cols = starts_with("X"))

    return(out)

}

pentad_vars <- purrr::map(bands, ~f(.x))

out <- Reduce("cbind", pentad_vars)

# Fix covariates
sitedata <- sitedata %>%
    dplyr::select(!contains(".")) %>%
    dplyr::select(Name, everything()) %>%
    dplyr::select(-id) %>%
    rename_with(.fn =  ~gsub("pr", "prcp_", .x),
                .cols = starts_with("pr")) %>%
    rename_with(.fn =  ~gsub("tmmn", "tmmn_", .x),
                .cols = starts_with("tmmn")) %>%
    rename_with(.fn =  ~gsub("tmmx", "tmmx_", .x),
                .cols = starts_with("tmmx")) %>%
    mutate(across(.cols = starts_with("tmmn"),
                  .fns = ~.x/10),
           across(.cols = starts_with("tmmx"),
                  .fns = ~.x/10)) %>%
    st_drop_geometry()

saveRDS(out, file.path(config$data_dir, "site_dat_sa_gee_08_19.rds"))


# Annotate with yearly surface water occurrence --------------------------------

# Number of pixels with water each year
band <- "waterClass"
stackCollection <- EEcollectionToMultiband(collection = "JRC/GSW1_3/YearlyHistory",
                                           dates = c("2008-01-01", "2020-01-01"),
                                           band = band,
                                           group_type = "year",
                                           groups = 2008:2019,
                                           unmask = FALSE)

out <- addVarEEimage(ee_pentads, stackCollection, "count")

# Fix covariates
out <- out %>%
    dplyr::select(Name, everything()) %>%
    dplyr::select(-id) %>%
    rename_with(~gsub("X", "watext_", .x), .cols = starts_with("X")) %>%
    st_drop_geometry()

readRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds")) %>%
    left_join(out, by = "Name") %>%
    saveRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds"))

# Recurrence of pixels with water each year
out <- addVarEEimage(ee_pentads, stackCollection, "mean")

# Fix covariates
out <- out %>%
    dplyr::select(Name, everything()) %>%
    dplyr::select(-id) %>%
    rename_with(~gsub("X", "watrec_", .x), .cols = starts_with("X"))
mutate(across(.fns = ~tidyr::replace_na(.x, 0))) %>%
    st_drop_geometry()

readRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds")) %>%
    left_join(out, by = "Name") %>%
    saveRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds"))


# Annotate with overall surface water occurrence --------------------------

out <- addVarEEimage(ee_pentads = ee_pentads,
                     image = "JRC/GSW1_3/GlobalSurfaceWater",
                     reducer = "mean",
                     bands = "occurrence",
                     unmask = TRUE)

# Fix covariates
out <- out %>%
    dplyr::select(Name, everything()) %>%
    dplyr::select(-id) %>%
    rename(water_occur = occurrence) %>%
    rename_with(.fn =  ~gsub("water_occur", "watocc_ever", .x),
                .cols = starts_with("water_occur")) %>%
    st_drop_geometry()

readRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds")) %>%
    left_join(out, by = "Name") %>%
    saveRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds"))


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

# Fix covariates
out <- out %>%
    dplyr::select(Name, everything()) %>%
    dplyr::select(-id) %>%
    rename_with(~gsub("X", "ndvi_", .x), .cols = starts_with("X")) %>%
    mutate(across(.cols = starts_with("ndvi"),
                  .fns = ~.x/10000)) %>%
    st_drop_geometry()

readRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds")) %>%
    left_join(out, by = "Name") %>%
    saveRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds"))


# Add geometries -----------------------------------------------------------

# Load data
sitedata <- readRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds"))

# Add geometries
sitedata <- sitedata %>%
    left_join(dplyr::select(pentads_sa, Name))

# Add centroids
cc <- st_centroid(sitedata) %>%
    st_coordinates() %>%
    as.data.frame()

sitedata <- sitedata %>%
    mutate(lon = cc$X,
           lat = cc$Y)

saveRDS(sitedata, file.path(config$data_dir, "site_dat_sa_gee_08_19.rds"))


# Annotate with distance to coast -----------------------------------------

library(rnaturalearth)

pp <- readRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds"))

# Create bounding box
sabox <- st_bbox(pp) + c(-1, -1, 3, 1)

# Download coastline for the World and crop
coast <- ne_coastline(scale = 10, returnclass = "sf") %>%
    st_crop(sabox)

# Check that makes sense
plot(st_geometry(coast[1,]))
plot(st_geometry(pp), add = TRUE)

# Find distances in kilometers
pp$dist_coast <- as.numeric(st_distance(pp, coast[1,]))/1000 # in kilometers

plot(pp["dist_coast"])

# Save
saveRDS(pp, file.path(config$data_dir, "site_dat_sa_gee_08_19.rds"))




# Annotate with elevation -------------------------------------------------

out <- addVarEEimage(ee_pentads = ee_pentads,
                     image = "MERIT/DEM/v1_0_3",
                     reducer = "mean",
                     bands = "dem",
                     unmask = FALSE)

# Fix covariates
out <- out %>%
    dplyr::select(Name, elev = dem) %>%
    st_drop_geometry()

readRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds")) %>%
    left_join(out, by = "Name") %>%
    saveRDS(file.path(config$data_dir, "site_dat_sa_gee_08_19.rds"))
