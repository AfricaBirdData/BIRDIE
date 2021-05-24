# 24-05-2021

# In this script we download SABAP data for the Barberspan region (North West province)

library(SABAP)
library(tidyverse)
library(sf)
devtools::load_all() # To load BIRDIE package functions

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

pentads <- st_intersection(pentads, region) %>%
    dplyr::select(Name, lon, lat)

plot(st_geometry(pentads))
plot(st_geometry(region), add = T)


# Get species list --------------------------------------------------------

# Load CWAC data for Barberspan
cwacdata <- readRDS("analysis/data/barberspan.rds")

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


