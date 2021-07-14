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


# Download climatic data ---------------------------------------------

vars <- c("prcp", "tmax", "tmin", "aet", "pet")

# Download data from TerraClimate server
for(i in seq_along(vars)){
    dat <- getTerraClim(region,
                        param = vars[i],
                        startDate = "2003-01-01",
                        endDate = "2019-12-31")

    dat <- dat[[1]]

    fileout <- paste0("analysis/out_nosync/terraClim_", vars[i], "_03_19_nw.rds")
    saveRDS(dat, fileout)

}


# Calculate annual means for the 17 years ---------------------------------

ss <- vector("list", length = 17)

yy <- 2003
for(i in 1:17){
    ry <- raster::subset(dat, grep(yy, names(dat)))
    ry <- calc(ry, fun = mean)
    names(ry) <- yy
    ss[[i]] <- ry
    yy <- yy + 1
}

ss <- stack(ss)

plot(ss)

# Save
writeRaster(ss, "analysis/out_nosync/prcp_03_19_nw_prov.grd", format = "raster")
