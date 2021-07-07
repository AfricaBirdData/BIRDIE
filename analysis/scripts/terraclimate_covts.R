# 7-7-2021

# In this script we download and prepare TerraClimate covariates.
# For now we focus on Barberspan and the North West province

library(climateR)
library(raster)
library(sf)

rm(list = ls())


# Define region of interest -----------------------------------------------

region <- raster::getData("GADM", download = TRUE, country = "South Africa",
                          level = 1, path = "analysis/data") %>%
    sf::st_as_sf() %>%
    dplyr::filter(NAME_1 == "North West")


# Download precipitation data ---------------------------------------------

# Download data from TerraClimate server
dat <- getTerraClim(region,
                  param = "prcp",
                  startDate = "1992-01-01",
                  endDate = "2019-12-31")

dat <- dat[[1]]


# Calculate annual means for the 27 years ---------------------------------

ss <- vector("list", length = 27)

yy <- 1992
for(i in 1:27){
    ry <- raster::subset(dat, grep(yy, names(dat)))
    ry <- calc(ry, fun = mean)
    names(ry) <- yy
    ss[[i]] <- ry
    yy <- yy + 1
}

ss <- stack(ss)

plot(ss)

# Save
writeRaster(ss, "analysis/out_nosync/prcp_92_19_nw_prov.grd", format = "raster")
