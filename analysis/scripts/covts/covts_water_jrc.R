# In this script we download surface water data from the JRC website.

# For now they don't have the time series for the surface water only single
# map summaries of different types.

# We can try rgee Google Earth Engine for R to download the time series if we
# want

library(raster)

rm(list = ls())

# Set URL to download the watersurface occurrence map
cc <- list(c(10, 20),
           c(20, 20),
           c(30, 20),
           c(10, 30),
           c(20, 30),
           c(30, 30))

# url <- "https://storage.googleapis.com/global-surface-water/downloads2020/occurrence/occurrence_20E_20Sv1_3_2020.tif"

for(i in seq_along(cc)){

    url <- paste0("https://storage.googleapis.com/global-surface-water/downloads2020/occurrence/occurrence_", cc[[i]][1], "E_", cc[[i]][2], "Sv1_3_2020.tif")

    # Download
    download.file(url, destfile = paste0("analysis/downloads/surf_water_", cc[[i]][1], "E_", cc[[i]][2], "S.tif"))

}

# Visualize
water <- raster("analysis/downloads/surf_water_20E_20S.tif")
plot(water)

# Resolution
res(water)


# Merge all water layers --------------------------------------------------

ff <- list.files("analysis/downloads", pattern = "surf_water")

rr <- vector("list", length = length(ff))

for(i in seq_along(ff)){
    rr[[i]] <- raster(paste0("analysis/downloads/", ff[i]))
}

r <- do.call(merge, rr)


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

