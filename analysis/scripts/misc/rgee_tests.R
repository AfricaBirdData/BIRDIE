library(terra)
library(sf)
library(dplyr)

rm(list = ls())

tmin <- readRDS("analysis/downloads/terraClim_tmin_03_19.rds")
prcp <- readRDS("analysis/downloads/terraClim_prcp_03_19.rds")

# subset 2010
tmin10 <- terra::subset(tmin, grep("2015", names(tmin)))
prcp10 <- terra::subset(prcp, grep("2015", names(prcp)))

# Find the mean for 2010
tmin_mean <- mean(tmin10)
prcp_mean <- mean(prcp10)


# Test GEE collection annotation ------------------------------------------

# Annotate only with tmmn 2010
library(rgee)
ee_check()
ee_Initialize(drive = TRUE)

# Load to EE (if not done already)
assetId <- sprintf("%s/%s", ee_get_assethome(), 'pentads')
ee_pentads <- ee$FeatureCollection(assetId)

# Test minimum temperature
rgee_tmmn <- SABAP::addVarEEcollection(ee_pentads = ee_pentads,
                                       collection = "IDAHO_EPSCOR/TERRACLIMATE",
                                       dates = c("2010-01-01", "2011-01-01"),
                                       temp_reducer = "mean",
                                       spt_reducer = "mean",
                                       bands = "tmmn")

# Annotate with mean tmin
r_tmin <- rgee_tmmn %>%
    dplyr::select(Name, geometry)

r_tmin$tmin <- exactextractr::exact_extract(tmin_mean, r_tmin, fun = 'mean', weights = 'area')

# Compare
plot(rgee_tmmn$tmmn/10, r_tmin$tmin)


# Test precipitation
rgee_pr <- SABAP::addVarEEcollection(ee_pentads = ee_pentads,
                                       collection = "IDAHO_EPSCOR/TERRACLIMATE",
                                       dates = c("2010-01-01", "2011-01-01"),
                                       temp_reducer = "mean",
                                       spt_reducer = "mean",
                                       bands = "pr")

# Annotate with mean prcp
r_prcp <- rgee_pr %>%
    dplyr::select(Name, geometry)

r_prcp$prcp <- exactextractr::exact_extract(prcp_mean, r_prcp, fun = 'mean', weights = 'area')

# Compare
plot(rgee_pr$pr, r_prcp$prcp)

# Rain seems to come out slightly different from GEE, but not too worrying.


# Test multiple GEE annotation --------------------------------------------

# load data annotated with GEE
gee_data <- readRDS("analysis/data/site_dat_gee_08_19.rds")

# Annotate with mean tmin
r_tmin <- gee_data %>%
    dplyr::select(Name, geometry)

r_tmin$tmin <- exactextractr::exact_extract(tmin_mean, r_tmin, fun = 'mean', weights = 'area')

# Compare
plot(gee_data$tmmn2015/10, r_tmin$tmin)
# YES!

# Annotate with mean prcp
r_prcp <- gee_data %>%
    dplyr::select(Name, geometry)

r_prcp$prcp <- exactextractr::exact_extract(prcp_mean, r_prcp, fun = 'mean', weights = 'area')

# Compare
plot(gee_data$pr2015, r_prcp$prcp)
# YES!
