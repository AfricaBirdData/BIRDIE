rm(list = ls())

library(BIRDIE)
library(sf)
library(tidyverse)

sf_use_s2(FALSE) # s2 intersection takes very long

# Find pentads in the North West province
pentads_nw <- getRegionPentads(country = "South Africa",
                               province = "North West",
                               path = "analysis/data")

plot(st_geometry(pentads_nw))

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

# Any pentads in occu_nw that not in pentads_nw?
any(unique(occu_nw$Pentad) %in% unique(pentads_nw$Name)) # YES, WHY? Check with Michael

# Filter out those pentads that are not in NW province
occu_nw <- occu_nw %>%
    dplyr::filter(Pentad %in% unique(pentads_nw$Name))

# Plot
det_nw <- occu_nw %>%
    filter(Spp != "-") %>%
    dplyr::select(Pentad, Spp) %>%
    distinct(Pentad, .keep_all = TRUE)

pentads_nw <- pentads_nw %>%
    left_join(det_nw, by = c("Name" = "Pentad"))

plot(pentads_nw["Spp"])






# # Load wetland layer
# wetlands <- readRDS("../resources/wetland_nw.rds")
#
# plot(st_geometry(wetlands))
#
# # Make a single feature
# wetlands_one <- st_union(wetlands)
#
# st_area(wetlands_one)
#
# wet_int <- st_intersection(st_transform(pentads_nw, st_crs(wetlands_one)),
#                            wetlands_one)

