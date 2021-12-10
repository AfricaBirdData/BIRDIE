
# In this script we prepare VISIT data for fitting a spatial occupancy model
# to a pilot species in South Africa

library(BIRDIE)
library(dplyr)

rm(list = ls())


# For now, we want to select species present at Barberspan
bbpan <- BIRDIE::barberspan %>%
    pull(spp) %>%
    unique()


# Select a species and a region -------------------------------------------

# We are most interested in the Maccoa Duck VU (103) and the Cape Cormorant EN (48)
# i <- 1
# sp_sel <- bbpan[i]
sp_sel <- 103

sp_name <- BIRDIE::barberspan %>%
    dplyr::filter(spp == sp_sel) %>%
    mutate(name = paste(taxon.Common_species, taxon.Common_group)) %>%
    pull(name) %>%
    unique()


# Load occupancy site data -------------------------------------------------

# Load site data
sitedata <- readRDS("analysis/data/site_dat_sa_wcovts_16_19.rds")


# Prepare occupancy visit data --------------------------------------------

future::plan("multisession", workers = 6)

visitdata <- prepOccVisitData(region = "South Africa",
                              sites = sitedata,
                              species = sp_sel,
                              years = 2016:2019,
                              clim_covts = c("prcp", "tmax", "tmin", "aet", "pet"),
                              covts_dir = "analysis/downloads/",
                              file_fix = c("terraClim_", "_03_19"),
                              savedir = "analysis/data/pentads_sa.rds")

future::plan("sequential")


saveRDS(visitdata, paste0("analysis/data/visit_dat_", sp_sel, "_wcovts_16_19.rds"))


