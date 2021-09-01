
# In this script we prepare VISIT data for fitting a spatial occupancy model
# to a pilot species in South Africa

library(BIRDIE)

rm(list = ls())


future::plan("multisession", workers = 6)

sp_sel <- 6

sites <- readRDS("analysis/data/site_dat_sa_wcovts.rds")

pa_dat <- prepOccVisitData(region = "South Africa",
                           sites = sites,
                           species = sp_sel,
                           years = 2008:2011,
                           clim_covts = c("prcp", "tmax", "tmin", "aet", "pet"),
                           covts_dir = "analysis/downloads/",
                           file_fix = c("terraClim_", "_03_19"),
                           savedir = "analysis/data/pentads_sa.rds")

future::plan("sequential")


saveRDS(pa_dat, paste0("analysis/data/visit_dat_", sp_sel, "_wcovts.rds"))


