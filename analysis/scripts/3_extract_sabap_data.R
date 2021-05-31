# 24-05-2021

# In this script we download SABAP data for the Barberspan region (North West province)

library(SABAP)
library(tidyverse)
library(sf)
library(BIRDIE)
library(raster)

rm(list = ls())


# Prepare background data -------------------------------------------------

# Download/load geographic data
region <- raster::getData("GADM", download = FALSE, country = "South Africa", level = 1,
                          path = "analysis/data") %>% st_as_sf()

# Load pentads
pentads <- read.csv("analysis/data/pentads_sabap.csv") %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = F)

# Cut to study area
region <- region %>%
    filter(NAME_1 == "North West")

pentads <- st_crop(pentads, extent(region))

# pentads <- st_intersection(pentads, region) %>%
    # dplyr::select(Name, lon, lat)

plot(st_geometry(pentads))
plot(st_geometry(region), add = T)


# Save North West pentads in raster format --------------------------------

pentads_m <- pentads %>%
    transmute(x = lon,
              y = lat,
              pentad_id = row_number()) %>%
    st_drop_geometry() %>%
    as.matrix()

pixels <- sp::SpatialPixelsDataFrame(pentads_m, tolerance = 0.1, as.data.frame(pentads_m))

pentads_r <- raster(pixels, values = FALSE)
values(pentads_r) <- pentads_m[,"pentad_id"]

plot(pentads_r)

attr(pentads_r, "pentads") <- pentads$Name

saveRDS(pentads_r, "analysis/data/pentads_nw.rds")


# Get species list --------------------------------------------------------

# Load CWAC data for Barberspan
cwacdata <- readRDS("analysis/data/cwacdata.rds")

# Seems like CWAC and SABAP species codes are the same
cwacdata %>%
    dplyr::select(spp, taxon.Common_group, taxon.Common_species)

SABAP::searchSabapSpecies(6)

# Extract species list
spplist <- tibble(spp = unique(cwacdata$spp))

# Add conservation status
spplist <- left_join(spplist %>%
                         mutate(spp = as.character(spp)),
                     CWAC::searchCwacTerm("SpeciesList")$options,
                     by = c("spp" = "ref"))

# Select species with unfavorable conservation status
spp_sel <- spplist %>%
    filter(iucn_status != "")


# Download SABAP data ------------------------------------------------------

# Downloas SABAP data for a test species
sabapdata <- SABAP::getSabapData(spp_sel$spp[1], region_type = "province", region = "North West")


# Prepare occupancy data --------------------------------------------------

occudata <- sabapdata %>%
    mutate(year = lubridate::year(StartDate),
           detc = if_else(Spp == "-", 0, 1)) %>%
    dplyr::select(year, Pentad, TotalHours, detc) %>%
    left_join(st_drop_geometry(pentads), by = c("Pentad" = "Name"))

# THERE ARE A NUMBER OF PENTADS THAT SABAP CONSIDERS IN NORTH WEST THAT ARE NOT
# PRESENT IN THE REGION OBJECT. FIX!!!

occudata <- occudata %>%
    filter(!is.na(lon))

occudata %>%
    group_by(Pentad) %>%
    summarize(n = sum(detc),
              lon = first(lon),
              lat = first(lat)) %>%
    mutate(detc = if_else(n == 0, 0, 1)) %>%
    ungroup() %>%
    ggplot() +
    geom_tile(aes(x = lon, y = lat, fill = factor(detc))) +
    geom_tile(data = pentads, aes(x = lon, y = lat), fill = NA, col = "black") +
    scale_fill_viridis_d() +
    coord_equal()

# Filter out pentads where the species has never been observed
# pentads_sel <- occudata %>%
#     filter(detc == 1) %>%
#     pull(Pentad) %>%
#     unique()
#
# occudata <- occudata %>%
#     filter(Pentad %in% pentads_sel)


# Save data ---------------------------------------------------------------

saveRDS(occudata, "analysis/data/occudata.rds")
