rm(list = ls())

library(BIRDIE)
library(sf)
library(tidyverse)


# Prepare data ------------------------------------------------------------

# Find pentads in the North West province
pentads_nw <- readRDS("analysis/data/pentads_nw.rds")

# OR

# sf_use_s2(FALSE) # s2 intersection takes very long
# pentads_nw <- getRegionPentads(country = "South Africa",
#                                province = "North West",
#                                path = "analysis/data")
#
# plot(st_geometry(pentads_nw))
#
# # Add centroid coordinates
# sf_use_s2(TRUE)
# pentads_nw <- pentads_nw %>%
#     mutate(lon = st_coordinates(st_centroid(.))[,1],
#            lat = st_coordinates(st_centroid(.))[,2])
#
# saveRDS(pentads_nw, "analysis/data/pentads_nw.rds")

# Download SABAP data for the region and Baberspan species
spp <- barberspan %>%
    dplyr::pull(spp) %>%
    unique()

# For all species run the following. It takes quite long. Consider furrr?
# occu_nw <- purrr::map_dfr(spp, ~SABAP::getSabapData(.x,
#                                                     region_type = "province",
#                                                     region = "North West"))

# For a single species
occu_nw <- SABAP::getSabapData(spp[1], region_type = "province",
                               region = "North West")

occu_nw <- occu_nw %>%
    mutate(year = lubridate::year(StartDate))

# Filter out those pentads that are not in NW province
occu_nw <- occu_nw %>%
    dplyr::filter(Pentad %in% unique(pentads_nw$Name))

# Join occu data and geographic data
occu_nw <- occu_nw %>%
    left_join(pentads_nw %>%
                  st_drop_geometry() %>%
                  dplyr::select(Name, lon, lat),
              by = c("Pentad" = "Name")) %>%
    mutate(detc = if_else(Spp == "-", 0, 1),
           year = lubridate::year(StartDate)) %>%
    dplyr::select(CardNo, lon, lat, year, StartDate, EndDate, Pentad, TotalHours, detc)

# There are no covariates for 2020 so...
occu_nw <- occu_nw %>%
    filter(year < 2020)

saveRDS(occu_nw, "analysis/data/occudata.rds")


# Plot detection map
det_nw <- occu_nw %>%
    filter(detc == 1) %>%
    dplyr::select(Pentad, detc) %>%
    distinct(Pentad, .keep_all = TRUE)

pentads_nw <- pentads_nw %>%
    left_join(det_nw, by = c("Name" = "Pentad"))

plot(pentads_nw["detc"])


# Prepare site covariates -------------------------------------------------

library(BIRDIE)
library(dplyr)
library(raster)
library(sf)

rm(list = ls())
gc()

# Load pentads in the North West province
pentads_nw <- readRDS("analysis/data/pentads_nw.rds")

# Load occupancy data
occu_nw <- readRDS("analysis/data/occudata.rds")

# Load covariates
prcp <- raster::stack("analysis/out_nosync/prcp_03_19_nw_prov.gri")

# Define sites and years
sites <- occu_nw %>%
    distinct(Pentad, year)

years <- sort(unique(sites$year))

prcp <- subset(prcp, paste0("X", years))

pent_prcp <- raster::extract(prcp, as(pentads_nw, "Spatial"), fun = mean)

pent_prcp <- cbind(pentads_nw, pent_prcp)

saveRDS(pent_prcp, "analysis/data/pentd_nw_prcp.rds")


# Extract covariates (this could be a function and run in multicore)

f <- function(yr, covt, sites, pentads){

    rr <- covt[[grep(yr, names(covt))]]

    sites_yr <- sites %>%
        filter(year == yr)

    pent_yr <- pentads %>%
        filter(Name %in% unique(sites_yr$Pentad)) %>%
        mutate(avg_prcp = raster::extract(rr, as(.,"Spatial"), fun = mean),
               year = yr)

    out <- pent_yr %>%
        st_drop_geometry() %>%
        dplyr::select(Pentad = Name, lon, lat, year, avg_prcp)

    return(out)

}

future::plan("multisession", workers = 5)
pent_covt <- furrr::future_map2_dfr(years, as.list(prcp), ~f(.x, .y, sites, pentads_nw))
future::plan("sequential")

saveRDS(pent_covt, "analysis/data/pentd_nw_prcp.rds")
