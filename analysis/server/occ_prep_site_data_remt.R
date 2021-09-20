
# In this script we prepare SITE data for fitting a spatial occupancy model
# to a pilot species in South Africa

library(BIRDIE)

birdie_dir <- "/home/birdie/analysis/"

rm(list = ls())

future::plan("multisession", workers = 6)

sites <- readRDS("analysis/data/pentads_sa.rds")

sites <- prepOccSiteData(#region = "South Africa",
                         sites = sites,
                         years = 2008:2019,
                         clim_covts = c("prcp", "tmax", "tmin", "aet", "pet"),
                         covts_dir = paste0(birdie_dir, "downloads/"),
                         file_fix = c("terraClim_", "_03_19"),
                         savedir = paste0(birdie_dir, "data/pentads_sa.rds"))

future::plan("sequential")

saveRDS(sites, paste0(birdie_dir, "data/site_dat_sa_wcovts_16_19.rds"))
