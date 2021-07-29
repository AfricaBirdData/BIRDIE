# In this script we download surface water data from the JRC website.

# For now they don't have the time series for the surface water only single
# map summaries of different types.

# We can try rgee Google Earth Engine for R to download the time series if we
# want

rm(list = ls())

# Set URL to download the watersurface occurrence map
url <- "https://storage.googleapis.com/global-surface-water/downloads2020/occurrence/occurrence_20E_20Sv1_3_2020.tif"

# Download
download.file(url, destfile = "analysis/data/surf_water_20E_20S.tif")

# Visualize
library(raster)

water <- raster("analysis/data/surf_water_20E_20S.tif")
plot(water)


# Subset to North West province -------------------------------------------
library(sf)

pentads_nw <- BIRDIE::getRegionPentads("South Africa", "North West", .path = "analysis/data/")

plot(st_geometry(pentads_nw))

water <- raster::crop(water, pentads_nw)

plot(water)
plot(st_geometry(pentads_nw), fill = NA, add = TRUE)

future::plan("multisession", workers = 6)
pentads_water <- BIRDIE::exactExtractParll(water, pentads_nw,
                                           ncores = future::nbrOfWorkers(),
                                           fun = "sum")
future::plan("sequential")

