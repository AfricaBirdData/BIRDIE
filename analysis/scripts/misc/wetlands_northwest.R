library(sf)
library(tidyverse)


wetland <- st_read("../resources/NBA2018_National Wetland Map 5 and confidence map/NBA2018_NWM5_AEA.shp")

# Download/load geographic data
region <- raster::getData("GADM", download = FALSE, country = "South Africa", level = 1,
                          path = "analysis/data") %>% st_as_sf()

# Cut to study area
region <- region %>%
    filter(NAME_1 == "North West")

region <- st_transform(region, st_crs(wetland))

wetland_nw <- st_intersection(wetland, region)

saveRDS(wetland_nw, "../resources/wetland_nw.rds")

plot(st_geometry(pentads))
plot(st_geometry(region), add = T)
