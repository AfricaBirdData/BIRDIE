library(rgee)
library(dplyr)
library(stars)
library(raster)

ee_check()
ee_Initialize()

rm(list = ls())

# Upload National Wetland Map ---------------------------------------------

# Read in NWM
nwm <- raster("analysis/data/NWM5_WETCON_100ddw_reclass/NWM5_WETCON_100ddw_reclass.tif")
# nwm <- read_stars("analysis/data/NWM5_WETCON_100ddw_reclass/NWM5_WETCON_100ddw_reclass.tif")
nwm[nwm < 0] <- NA
plot(nwm)
writeRaster(nwm, "analysis/data/NWM5_WETCON_100ddw_reclass/NWM5_WETCON_100ddw_na.tif")

nwm <- st_as_stars(nwm)

# Set a name for the asset
eenwm_id <- file.path(rgee::ee_get_assethome(), 'wetland_map_sa')

# Uploading to earth egnine
ee_upload(x = nwm,
          filename = eenwm_id)
